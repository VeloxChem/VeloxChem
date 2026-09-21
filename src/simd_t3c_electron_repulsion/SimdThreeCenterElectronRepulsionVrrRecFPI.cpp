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


#include "SimdThreeCenterElectronRepulsionVrrRecFPI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
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

    const auto *dpi0_0 = buffer.data(dpi0 + 0);
    const auto *dpi0_3 = buffer.data(dpi0 + 3);
    const auto *dpi0_5 = buffer.data(dpi0 + 5);
    const auto *dpi0_6 = buffer.data(dpi0 + 6);
    const auto *dpi0_9 = buffer.data(dpi0 + 9);
    const auto *dpi0_10 = buffer.data(dpi0 + 10);
    const auto *dpi0_14 = buffer.data(dpi0 + 14);
    const auto *dpi0_20 = buffer.data(dpi0 + 20);
    const auto *dpi0_27 = buffer.data(dpi0 + 27);

    const auto *dph_0 = buffer.data(dph + 0);
    const auto *dph_1 = buffer.data(dph + 1);
    const auto *dph_3 = buffer.data(dph + 3);
    const auto *dph_5 = buffer.data(dph + 5);
    const auto *dph_6 = buffer.data(dph + 6);
    const auto *dph_9 = buffer.data(dph + 9);
    const auto *dph_15 = buffer.data(dph + 15);
    const auto *dph_17 = buffer.data(dph + 17);
    const auto *dph_18 = buffer.data(dph + 18);
    const auto *dph_20 = buffer.data(dph + 20);
    const auto *dph_21 = buffer.data(dph + 21);
    const auto *dph_26 = buffer.data(dph + 26);
    const auto *dph_30 = buffer.data(dph + 30);
    const auto *dph_36 = buffer.data(dph + 36);
    const auto *dph_38 = buffer.data(dph + 38);
    const auto *dph_39 = buffer.data(dph + 39);
    const auto *dph_41 = buffer.data(dph + 41);
    const auto *dph_57 = buffer.data(dph + 57);
    const auto *dph_59 = buffer.data(dph + 59);
    const auto *dph_60 = buffer.data(dph + 60);
    const auto *dph_62 = buffer.data(dph + 62);
    const auto *dph_78 = buffer.data(dph + 78);
    const auto *dph_80 = buffer.data(dph + 80);
    const auto *dph_81 = buffer.data(dph + 81);
    const auto *dph_82 = buffer.data(dph + 82);
    const auto *dph_84 = buffer.data(dph + 84);
    const auto *dph_87 = buffer.data(dph + 87);
    const auto *dph_90 = buffer.data(dph + 90);
    const auto *dph_94 = buffer.data(dph + 94);

    const auto *dpi1_0 = buffer.data(dpi1 + 0);
    const auto *dpi1_3 = buffer.data(dpi1 + 3);
    const auto *dpi1_5 = buffer.data(dpi1 + 5);
    const auto *dpi1_6 = buffer.data(dpi1 + 6);
    const auto *dpi1_9 = buffer.data(dpi1 + 9);
    const auto *dpi1_10 = buffer.data(dpi1 + 10);
    const auto *dpi1_14 = buffer.data(dpi1 + 14);
    const auto *dpi1_20 = buffer.data(dpi1 + 20);
    const auto *dpi1_27 = buffer.data(dpi1 + 27);

    const auto *fsi0_0 = buffer.data(fsi0 + 0);
    const auto *fsi0_3 = buffer.data(fsi0 + 3);
    const auto *fsi0_5 = buffer.data(fsi0 + 5);
    const auto *fsi0_6 = buffer.data(fsi0 + 6);
    const auto *fsi0_9 = buffer.data(fsi0 + 9);
    const auto *fsi0_10 = buffer.data(fsi0 + 10);
    const auto *fsi0_12 = buffer.data(fsi0 + 12);
    const auto *fsi0_14 = buffer.data(fsi0 + 14);
    const auto *fsi0_21 = buffer.data(fsi0 + 21);
    const auto *fsi0_23 = buffer.data(fsi0 + 23);
    const auto *fsi0_24 = buffer.data(fsi0 + 24);
    const auto *fsi0_25 = buffer.data(fsi0 + 25);
    const auto *fsi0_27 = buffer.data(fsi0 + 27);

    const auto *fsh_0 = buffer.data(fsh + 0);
    const auto *fsh_1 = buffer.data(fsh + 1);
    const auto *fsh_2 = buffer.data(fsh + 2);
    const auto *fsh_3 = buffer.data(fsh + 3);
    const auto *fsh_5 = buffer.data(fsh + 5);
    const auto *fsh_6 = buffer.data(fsh + 6);
    const auto *fsh_8 = buffer.data(fsh + 8);
    const auto *fsh_9 = buffer.data(fsh + 9);
    const auto *fsh_10 = buffer.data(fsh + 10);
    const auto *fsh_14 = buffer.data(fsh + 14);
    const auto *fsh_15 = buffer.data(fsh + 15);
    const auto *fsh_17 = buffer.data(fsh + 17);
    const auto *fsh_18 = buffer.data(fsh + 18);
    const auto *fsh_19 = buffer.data(fsh + 19);
    const auto *fsh_20 = buffer.data(fsh + 20);
    const auto *fsh_21 = buffer.data(fsh + 21);
    const auto *fsh_26 = buffer.data(fsh + 26);
    const auto *fsh_30 = buffer.data(fsh + 30);
    const auto *fsh_36 = buffer.data(fsh + 36);
    const auto *fsh_38 = buffer.data(fsh + 38);
    const auto *fsh_39 = buffer.data(fsh + 39);
    const auto *fsh_40 = buffer.data(fsh + 40);

    const auto *fsi1_0 = buffer.data(fsi1 + 0);
    const auto *fsi1_3 = buffer.data(fsi1 + 3);
    const auto *fsi1_5 = buffer.data(fsi1 + 5);
    const auto *fsi1_6 = buffer.data(fsi1 + 6);
    const auto *fsi1_9 = buffer.data(fsi1 + 9);
    const auto *fsi1_10 = buffer.data(fsi1 + 10);
    const auto *fsi1_12 = buffer.data(fsi1 + 12);
    const auto *fsi1_14 = buffer.data(fsi1 + 14);
    const auto *fsi1_21 = buffer.data(fsi1 + 21);
    const auto *fsi1_23 = buffer.data(fsi1 + 23);
    const auto *fsi1_24 = buffer.data(fsi1 + 24);
    const auto *fsi1_25 = buffer.data(fsi1 + 25);
    const auto *fsi1_27 = buffer.data(fsi1 + 27);

    const auto *fpg0_0 = buffer.data(fpg0 + 0);
    const auto *fpg0_1 = buffer.data(fpg0 + 1);
    const auto *fpg0_2 = buffer.data(fpg0 + 2);
    const auto *fpg0_3 = buffer.data(fpg0 + 3);
    const auto *fpg0_5 = buffer.data(fpg0 + 5);
    const auto *fpg0_10 = buffer.data(fpg0 + 10);
    const auto *fpg0_12 = buffer.data(fpg0 + 12);
    const auto *fpg0_13 = buffer.data(fpg0 + 13);
    const auto *fpg0_14 = buffer.data(fpg0 + 14);
    const auto *fpg0_35 = buffer.data(fpg0 + 35);
    const auto *fpg0_42 = buffer.data(fpg0 + 42);
    const auto *fpg0_43 = buffer.data(fpg0 + 43);
    const auto *fpg0_44 = buffer.data(fpg0 + 44);
    const auto *fpg0_48 = buffer.data(fpg0 + 48);
    const auto *fpg0_55 = buffer.data(fpg0 + 55);
    const auto *fpg0_56 = buffer.data(fpg0 + 56);
    const auto *fpg0_57 = buffer.data(fpg0 + 57);
    const auto *fpg0_60 = buffer.data(fpg0 + 60);
    const auto *fpg0_62 = buffer.data(fpg0 + 62);
    const auto *fpg0_63 = buffer.data(fpg0 + 63);
    const auto *fpg0_65 = buffer.data(fpg0 + 65);
    const auto *fpg0_66 = buffer.data(fpg0 + 66);
    const auto *fpg0_70 = buffer.data(fpg0 + 70);

    const auto *fpg1_0 = buffer.data(fpg1 + 0);
    const auto *fpg1_1 = buffer.data(fpg1 + 1);
    const auto *fpg1_2 = buffer.data(fpg1 + 2);
    const auto *fpg1_3 = buffer.data(fpg1 + 3);
    const auto *fpg1_5 = buffer.data(fpg1 + 5);
    const auto *fpg1_10 = buffer.data(fpg1 + 10);
    const auto *fpg1_12 = buffer.data(fpg1 + 12);
    const auto *fpg1_13 = buffer.data(fpg1 + 13);
    const auto *fpg1_14 = buffer.data(fpg1 + 14);
    const auto *fpg1_35 = buffer.data(fpg1 + 35);
    const auto *fpg1_42 = buffer.data(fpg1 + 42);
    const auto *fpg1_43 = buffer.data(fpg1 + 43);
    const auto *fpg1_44 = buffer.data(fpg1 + 44);
    const auto *fpg1_48 = buffer.data(fpg1 + 48);
    const auto *fpg1_55 = buffer.data(fpg1 + 55);
    const auto *fpg1_56 = buffer.data(fpg1 + 56);
    const auto *fpg1_57 = buffer.data(fpg1 + 57);
    const auto *fpg1_60 = buffer.data(fpg1 + 60);
    const auto *fpg1_62 = buffer.data(fpg1 + 62);
    const auto *fpg1_63 = buffer.data(fpg1 + 63);
    const auto *fpg1_65 = buffer.data(fpg1 + 65);
    const auto *fpg1_66 = buffer.data(fpg1 + 66);
    const auto *fpg1_70 = buffer.data(fpg1 + 70);

    const auto *fph_0 = buffer.data(fph + 0);
    const auto *fph_1 = buffer.data(fph + 1);
    const auto *fph_2 = buffer.data(fph + 2);
    const auto *fph_3 = buffer.data(fph + 3);
    const auto *fph_5 = buffer.data(fph + 5);
    const auto *fph_6 = buffer.data(fph + 6);
    const auto *fph_8 = buffer.data(fph + 8);
    const auto *fph_9 = buffer.data(fph + 9);
    const auto *fph_10 = buffer.data(fph + 10);
    const auto *fph_14 = buffer.data(fph + 14);
    const auto *fph_15 = buffer.data(fph + 15);
    const auto *fph_17 = buffer.data(fph + 17);
    const auto *fph_18 = buffer.data(fph + 18);
    const auto *fph_19 = buffer.data(fph + 19);
    const auto *fph_20 = buffer.data(fph + 20);
    const auto *fph_21 = buffer.data(fph + 21);
    const auto *fph_23 = buffer.data(fph + 23);
    const auto *fph_24 = buffer.data(fph + 24);
    const auto *fph_26 = buffer.data(fph + 26);
    const auto *fph_27 = buffer.data(fph + 27);
    const auto *fph_30 = buffer.data(fph + 30);
    const auto *fph_31 = buffer.data(fph + 31);
    const auto *fph_35 = buffer.data(fph + 35);
    const auto *fph_36 = buffer.data(fph + 36);
    const auto *fph_38 = buffer.data(fph + 38);
    const auto *fph_39 = buffer.data(fph + 39);
    const auto *fph_41 = buffer.data(fph + 41);
    const auto *fph_42 = buffer.data(fph + 42);
    const auto *fph_44 = buffer.data(fph + 44);
    const auto *fph_45 = buffer.data(fph + 45);
    const auto *fph_47 = buffer.data(fph + 47);
    const auto *fph_48 = buffer.data(fph + 48);
    const auto *fph_50 = buffer.data(fph + 50);
    const auto *fph_51 = buffer.data(fph + 51);
    const auto *fph_52 = buffer.data(fph + 52);
    const auto *fph_56 = buffer.data(fph + 56);
    const auto *fph_57 = buffer.data(fph + 57);
    const auto *fph_59 = buffer.data(fph + 59);
    const auto *fph_60 = buffer.data(fph + 60);
    const auto *fph_61 = buffer.data(fph + 61);
    const auto *fph_62 = buffer.data(fph + 62);
    const auto *fph_63 = buffer.data(fph + 63);
    const auto *fph_64 = buffer.data(fph + 64);
    const auto *fph_66 = buffer.data(fph + 66);
    const auto *fph_68 = buffer.data(fph + 68);
    const auto *fph_69 = buffer.data(fph + 69);
    const auto *fph_70 = buffer.data(fph + 70);
    const auto *fph_72 = buffer.data(fph + 72);
    const auto *fph_73 = buffer.data(fph + 73);
    const auto *fph_78 = buffer.data(fph + 78);
    const auto *fph_79 = buffer.data(fph + 79);
    const auto *fph_80 = buffer.data(fph + 80);
    const auto *fph_81 = buffer.data(fph + 81);
    const auto *fph_82 = buffer.data(fph + 82);
    const auto *fph_83 = buffer.data(fph + 83);
    const auto *fph_84 = buffer.data(fph + 84);
    const auto *fph_85 = buffer.data(fph + 85);
    const auto *fph_86 = buffer.data(fph + 86);
    const auto *fph_87 = buffer.data(fph + 87);
    const auto *fph_89 = buffer.data(fph + 89);
    const auto *fph_90 = buffer.data(fph + 90);
    const auto *fph_91 = buffer.data(fph + 91);
    const auto *fph_93 = buffer.data(fph + 93);
    const auto *fph_94 = buffer.data(fph + 94);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, dph_0, fsh_0, fpg0_0, \
                         fpg1_0, fph_0, fph_1, fph_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * dph_0[k]
                 + f_1 * fsh_0[k]
                 + f_2 * fpg0_0[k]
                 - f_3 * fpg1_0[k]
                 + f_4 * pc_x[k] * fph_0[k];

        t_1[k] = f_4 * pc_y[k] * fph_0[k];

        t_2[k] = f_4 * pc_z[k] * fph_0[k];

        t_3[k] = f_5 * fpg0_0[k]
                 - f_6 * fpg1_0[k]
                 + f_4 * pc_y[k] * fph_1[k];

        t_4[k] = f_4 * pc_y[k] * fph_2[k];

        t_5[k] = f_5 * fpg0_0[k]
                 - f_6 * fpg1_0[k]
                 + f_4 * pc_z[k] * fph_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, fpg0_1, fpg0_2, fpg0_3, fpg1_1, \
                         fpg1_2, fpg1_3, fph_3, fph_5, fph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * fpg0_1[k]
                 - f_8 * fpg1_1[k]
                 + f_4 * pc_y[k] * fph_3[k];

        t_7[k] = f_4 * pc_z[k] * fph_3[k];

        t_8[k] = f_4 * pc_y[k] * fph_5[k];

        t_9[k] = f_7 * fpg0_2[k]
                 - f_8 * fpg1_2[k]
                 + f_4 * pc_z[k] * fph_5[k];

        t_10[k] = f_9 * fpg0_3[k]
                  - f_10 * fpg1_3[k]
                  + f_4 * pc_y[k] * fph_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, dph_15, fsh_15, \
                         fpg0_5, fpg1_5, fph_6, fph_8, fph_9, fph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * fph_6[k];

        t_12[k] = f_5 * fpg0_5[k]
                  - f_6 * fpg1_5[k]
                  + f_4 * pc_y[k] * fph_8[k];

        t_13[k] = f_4 * pc_y[k] * fph_9[k];

        t_14[k] = f_9 * fpg0_5[k]
                  - f_10 * fpg1_5[k]
                  + f_4 * pc_z[k] * fph_9[k];

        t_15[k] = f_0 * dph_15[k]
                  + f_1 * fsh_15[k]
                  + f_4 * pc_x[k] * fph_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, dph_17, dph_18, fsh_17, \
                         fsh_18, fph_10, fph_14, fph_17, fph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * fph_10[k];

        t_17[k] = f_0 * dph_17[k]
                  + f_1 * fsh_17[k]
                  + f_4 * pc_x[k] * fph_17[k];

        t_18[k] = f_0 * dph_18[k]
                  + f_1 * fsh_18[k]
                  + f_4 * pc_x[k] * fph_18[k];

        t_19[k] = f_4 * pc_y[k] * fph_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, dph_20, fsh_20, fpg0_10, \
                         fpg0_12, fpg1_10, fpg1_12, fph_15, fph_17, \
                         fph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * dph_20[k]
                  + f_1 * fsh_20[k]
                  + f_4 * pc_x[k] * fph_20[k];

        t_21[k] = f_2 * fpg0_10[k]
                  - f_3 * fpg1_10[k]
                  + f_4 * pc_y[k] * fph_15[k];

        t_22[k] = f_4 * pc_z[k] * fph_15[k];

        t_23[k] = f_9 * fpg0_12[k]
                  - f_10 * fpg1_12[k]
                  + f_4 * pc_y[k] * fph_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, fpg0_13, fpg0_14, fpg1_13, \
                         fpg1_14, fph_18, fph_19, fph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * fpg0_13[k]
                  - f_8 * fpg1_13[k]
                  + f_4 * pc_y[k] * fph_18[k];

        t_25[k] = f_5 * fpg0_14[k]
                  - f_6 * fpg1_14[k]
                  + f_4 * pc_y[k] * fph_19[k];

        t_26[k] = f_4 * pc_y[k] * fph_20[k];

        t_27[k] = f_2 * fpg0_14[k]
                  - f_3 * fpg1_14[k]
                  + f_4 * pc_z[k] * fph_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, fsi0_0, fsi0_3, fsh_0, \
                         fsh_1, fsi1_0, fsi1_3, fph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * fsi0_0[k]
                  - f_11 * pc_y[k] * fsi1_0[k];

        t_29[k] = f_1 * fsh_0[k]
                  + f_4 * pc_y[k] * fph_21[k];

        t_30[k] = f_4 * pc_z[k] * fph_21[k];

        t_31[k] = pb_y[k] * fsi0_3[k]
                  + f_12 * fsh_1[k]
                  - f_11 * pc_y[k] * fsi1_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, fsi0_5, fsi0_6, fsh_2, \
                         fsh_3, fsi1_5, fsi1_6, fph_23, fph_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * fsh_2[k]
                  + f_4 * pc_y[k] * fph_23[k];

        t_33[k] = pb_y[k] * fsi0_5[k]
                  - f_11 * pc_y[k] * fsi1_5[k];

        t_34[k] = pb_y[k] * fsi0_6[k]
                  + f_0 * fsh_3[k]
                  - f_11 * pc_y[k] * fsi1_6[k];

        t_35[k] = f_4 * pc_z[k] * fph_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, fsi0_9, fsi0_10, fsh_5, \
                         fsh_6, fsi1_9, fsi1_10, fph_26, fph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * fsh_5[k]
                  + f_4 * pc_y[k] * fph_26[k];

        t_37[k] = pb_y[k] * fsi0_9[k]
                  - f_11 * pc_y[k] * fsi1_9[k];

        t_38[k] = pb_y[k] * fsi0_10[k]
                  + f_13 * fsh_6[k]
                  - f_11 * pc_y[k] * fsi1_10[k];

        t_39[k] = f_4 * pc_z[k] * fph_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_x, pc_y, dph_36, fsi0_12, fsi0_14, \
                         fsh_8, fsh_9, fsi1_12, fsi1_14, fph_30, \
                         fph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * fsi0_12[k]
                  + f_12 * fsh_8[k]
                  - f_11 * pc_y[k] * fsi1_12[k];

        t_41[k] = f_1 * fsh_9[k]
                  + f_4 * pc_y[k] * fph_30[k];

        t_42[k] = pb_y[k] * fsi0_14[k]
                  - f_11 * pc_y[k] * fsi1_14[k];

        t_43[k] = f_0 * dph_36[k]
                  + f_4 * pc_x[k] * fph_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, dph_38, dph_39, fsh_14, \
                         fph_31, fph_35, fph_38, fph_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * pc_z[k] * fph_31[k];

        t_45[k] = f_0 * dph_38[k]
                  + f_4 * pc_x[k] * fph_38[k];

        t_46[k] = f_0 * dph_39[k]
                  + f_4 * pc_x[k] * fph_39[k];

        t_47[k] = f_1 * fsh_14[k]
                  + f_4 * pc_y[k] * fph_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pc_x, pc_y, pc_z, dph_41, fsi0_21, fsh_15, \
                         fsi1_21, fph_36, fph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * dph_41[k]
                  + f_4 * pc_x[k] * fph_41[k];

        t_49[k] = pb_y[k] * fsi0_21[k]
                  + f_14 * fsh_15[k]
                  - f_11 * pc_y[k] * fsi1_21[k];

        t_50[k] = f_4 * pc_z[k] * fph_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pc_y, fsi0_23, fsi0_24, fsi0_25, fsh_17, \
                         fsh_18, fsh_19, fsi1_23, fsi1_24, fsi1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_y[k] * fsi0_23[k]
                  + f_13 * fsh_17[k]
                  - f_11 * pc_y[k] * fsi1_23[k];

        t_52[k] = pb_y[k] * fsi0_24[k]
                  + f_0 * fsh_18[k]
                  - f_11 * pc_y[k] * fsi1_24[k];

        t_53[k] = pb_y[k] * fsi0_25[k]
                  + f_12 * fsh_19[k]
                  - f_11 * pc_y[k] * fsi1_25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_y, pb_z, pc_y, pc_z, fsi0_0, fsi0_27, \
                         fsh_20, fsi1_0, fsi1_27, fph_41, fph_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * fsh_20[k]
                  + f_4 * pc_y[k] * fph_41[k];

        t_55[k] = pb_y[k] * fsi0_27[k]
                  - f_11 * pc_y[k] * fsi1_27[k];

        t_56[k] = pb_z[k] * fsi0_0[k]
                  - f_11 * pc_z[k] * fsi1_0[k];

        t_57[k] = f_4 * pc_y[k] * fph_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_z, pc_y, pc_z, fsi0_3, fsi0_5, fsh_0, \
                         fsh_2, fsi1_3, fsi1_5, fph_42, fph_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * fsh_0[k]
                  + f_4 * pc_z[k] * fph_42[k];

        t_59[k] = pb_z[k] * fsi0_3[k]
                  - f_11 * pc_z[k] * fsi1_3[k];

        t_60[k] = f_4 * pc_y[k] * fph_44[k];

        t_61[k] = pb_z[k] * fsi0_5[k]
                  + f_12 * fsh_2[k]
                  - f_11 * pc_z[k] * fsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_z, pc_y, pc_z, fsi0_6, fsi0_9, fsh_3, \
                         fsh_5, fsi1_6, fsi1_9, fph_45, fph_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * fsi0_6[k]
                  - f_11 * pc_z[k] * fsi1_6[k];

        t_63[k] = f_1 * fsh_3[k]
                  + f_4 * pc_z[k] * fph_45[k];

        t_64[k] = f_4 * pc_y[k] * fph_47[k];

        t_65[k] = pb_z[k] * fsi0_9[k]
                  + f_0 * fsh_5[k]
                  - f_11 * pc_z[k] * fsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_y, pc_z, fsi0_10, fsh_6, fsi1_10, \
                         fpg0_35, fpg1_35, fph_48, fph_50, fph_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * fsi0_10[k]
                  - f_11 * pc_z[k] * fsi1_10[k];

        t_67[k] = f_1 * fsh_6[k]
                  + f_4 * pc_z[k] * fph_48[k];

        t_68[k] = f_5 * fpg0_35[k]
                  - f_6 * fpg1_35[k]
                  + f_4 * pc_y[k] * fph_50[k];

        t_69[k] = f_4 * pc_y[k] * fph_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_x, pc_z, dph_57, dph_59, fsi0_14, \
                         fsh_9, fsh_10, fsi1_14, fph_52, fph_57, \
                         fph_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * fsi0_14[k]
                  + f_13 * fsh_9[k]
                  - f_11 * pc_z[k] * fsi1_14[k];

        t_71[k] = f_0 * dph_57[k]
                  + f_4 * pc_x[k] * fph_57[k];

        t_72[k] = f_1 * fsh_10[k]
                  + f_4 * pc_z[k] * fph_52[k];

        t_73[k] = f_0 * dph_59[k]
                  + f_4 * pc_x[k] * fph_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_x, pc_y, pc_z, dph_60, dph_62, \
                         fsi0_21, fsi1_21, fph_56, fph_60, fph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * dph_60[k]
                  + f_4 * pc_x[k] * fph_60[k];

        t_75[k] = f_4 * pc_y[k] * fph_56[k];

        t_76[k] = f_0 * dph_62[k]
                  + f_4 * pc_x[k] * fph_62[k];

        t_77[k] = pb_z[k] * fsi0_21[k]
                  - f_11 * pc_z[k] * fsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, fsh_15, fpg0_42, fpg0_43, fpg1_42, \
                         fpg1_43, fph_57, fph_59, fph_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * fsh_15[k]
                  + f_4 * pc_z[k] * fph_57[k];

        t_79[k] = f_9 * fpg0_42[k]
                  - f_10 * fpg1_42[k]
                  + f_4 * pc_y[k] * fph_59[k];

        t_80[k] = f_7 * fpg0_43[k]
                  - f_8 * fpg1_43[k]
                  + f_4 * pc_y[k] * fph_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, pc_y, pc_z, fsi0_27, fsh_20, fsi1_27, \
                         fpg0_44, fpg1_44, fph_61, fph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * fpg0_44[k]
                  - f_6 * fpg1_44[k]
                  + f_4 * pc_y[k] * fph_61[k];

        t_82[k] = f_4 * pc_y[k] * fph_62[k];

        t_83[k] = pb_z[k] * fsi0_27[k]
                  + f_14 * fsh_20[k]
                  - f_11 * pc_z[k] * fsi1_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pa_y, pc_y, pc_z, dpi0_0, dpi0_3, \
                         dph_0, dph_1, dpi1_0, dpi1_3, fph_63, fph_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * dpi0_0[k]
                  - f_11 * pc_y[k] * dpi1_0[k];

        t_85[k] = f_1 * dph_0[k]
                  + f_4 * pc_y[k] * fph_63[k];

        t_86[k] = f_4 * pc_z[k] * fph_63[k];

        t_87[k] = pa_y[k] * dpi0_3[k]
                  + f_12 * dph_1[k]
                  - f_11 * pc_y[k] * dpi1_3[k];

        t_88[k] = f_4 * pc_z[k] * fph_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pa_y, pc_y, pc_z, dpi0_5, dpi0_6, dph_3, \
                         dph_5, dpi1_5, dpi1_6, fph_66, fph_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * dpi0_5[k]
                  - f_11 * pc_y[k] * dpi1_5[k];

        t_90[k] = pa_y[k] * dpi0_6[k]
                  + f_0 * dph_3[k]
                  - f_11 * pc_y[k] * dpi1_6[k];

        t_91[k] = f_4 * pc_z[k] * fph_66[k];

        t_92[k] = f_1 * dph_5[k]
                  + f_4 * pc_y[k] * fph_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pa_y, pc_y, pc_z, dpi0_9, dpi0_10, dph_6, \
                         dpi1_9, dpi1_10, fpg0_48, fpg1_48, fph_69, \
                         fph_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * dpi0_9[k]
                  - f_11 * pc_y[k] * dpi1_9[k];

        t_94[k] = pa_y[k] * dpi0_10[k]
                  + f_13 * dph_6[k]
                  - f_11 * pc_y[k] * dpi1_10[k];

        t_95[k] = f_4 * pc_z[k] * fph_69[k];

        t_96[k] = f_5 * fpg0_48[k]
                  - f_6 * fpg1_48[k]
                  + f_4 * pc_z[k] * fph_70[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, pa_y, pc_x, pc_y, pc_z, dpi0_14, dph_9, \
                         dph_78, dpi1_14, fsh_36, fph_72, fph_73, \
                         fph_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_1 * dph_9[k]
                  + f_4 * pc_y[k] * fph_72[k];

        t_98[k] = pa_y[k] * dpi0_14[k]
                  - f_11 * pc_y[k] * dpi1_14[k];

        t_99[k] = f_12 * dph_78[k]
                  + f_1 * fsh_36[k]
                  + f_4 * pc_x[k] * fph_78[k];

        t_100[k] = f_4 * pc_z[k] * fph_73[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pc_x, dph_80, dph_81, dph_82, fsh_38, fsh_39, \
                         fsh_40, fph_80, fph_81, fph_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_12 * dph_80[k]
                   + f_1 * fsh_38[k]
                   + f_4 * pc_x[k] * fph_80[k];

        t_102[k] = f_12 * dph_81[k]
                   + f_1 * fsh_39[k]
                   + f_4 * pc_x[k] * fph_81[k];

        t_103[k] = f_12 * dph_82[k]
                   + f_1 * fsh_40[k]
                   + f_4 * pc_x[k] * fph_82[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, t_107, pa_y, pc_y, pc_z, dpi0_20, dph_15, \
                         dpi1_20, fpg0_55, fpg1_55, fph_78, fph_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * dpi0_20[k]
                   - f_11 * pc_y[k] * dpi1_20[k];

        t_105[k] = f_1 * dph_15[k]
                   + f_2 * fpg0_55[k]
                   - f_3 * fpg1_55[k]
                   + f_4 * pc_y[k] * fph_78[k];

        t_106[k] = f_4 * pc_z[k] * fph_78[k];

        t_107[k] = f_5 * fpg0_55[k]
                   - f_6 * fpg1_55[k]
                   + f_4 * pc_z[k] * fph_79[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pc_y, pc_z, dph_20, fpg0_56, fpg0_57, fpg1_56, \
                         fpg1_57, fph_80, fph_81, fph_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_7 * fpg0_56[k]
                   - f_8 * fpg1_56[k]
                   + f_4 * pc_z[k] * fph_80[k];

        t_109[k] = f_9 * fpg0_57[k]
                   - f_10 * fpg1_57[k]
                   + f_4 * pc_z[k] * fph_81[k];

        t_110[k] = f_1 * dph_20[k]
                   + f_4 * pc_y[k] * fph_83[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_y, pc_x, pc_y, pc_z, dpi0_27, dph_21, \
                         dph_84, dpi1_27, fsh_21, fpg0_60, fpg1_60, \
                         fph_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_y[k] * dpi0_27[k]
                   - f_11 * pc_y[k] * dpi1_27[k];

        t_112[k] = f_12 * dph_84[k]
                   + f_2 * fpg0_60[k]
                   - f_3 * fpg1_60[k]
                   + f_4 * pc_x[k] * fph_84[k];

        t_113[k] = f_1 * dph_21[k]
                   + f_1 * fsh_21[k]
                   + f_4 * pc_y[k] * fph_84[k];

        t_114[k] = f_4 * pc_z[k] * fph_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pc_x, pc_z, dph_87, fpg0_60, fpg0_63, fpg1_60, \
                         fpg1_63, fph_85, fph_86, fph_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_12 * dph_87[k]
                   + f_9 * fpg0_63[k]
                   - f_10 * fpg1_63[k]
                   + f_4 * pc_x[k] * fph_87[k];

        t_116[k] = f_4 * pc_z[k] * fph_85[k];

        t_117[k] = f_5 * fpg0_60[k]
                   - f_6 * fpg1_60[k]
                   + f_4 * pc_z[k] * fph_86[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pc_x, pc_y, pc_z, dph_26, dph_90, fsh_26, \
                         fpg0_66, fpg1_66, fph_87, fph_89, fph_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_12 * dph_90[k]
                   + f_7 * fpg0_66[k]
                   - f_8 * fpg1_66[k]
                   + f_4 * pc_x[k] * fph_90[k];

        t_119[k] = f_4 * pc_z[k] * fph_87[k];

        t_120[k] = f_1 * dph_26[k]
                   + f_1 * fsh_26[k]
                   + f_4 * pc_y[k] * fph_89[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pc_x, pc_z, dph_94, fpg0_62, fpg0_70, fpg1_62, \
                         fpg1_70, fph_89, fph_90, fph_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_7 * fpg0_62[k]
                   - f_8 * fpg1_62[k]
                   + f_4 * pc_z[k] * fph_89[k];

        t_122[k] = f_12 * dph_94[k]
                   + f_5 * fpg0_70[k]
                   - f_6 * fpg1_70[k]
                   + f_4 * pc_x[k] * fph_94[k];

        t_123[k] = f_4 * pc_z[k] * fph_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pc_y, pc_z, dph_30, fsh_30, fpg0_63, fpg0_65, \
                         fpg1_63, fpg1_65, fph_91, fph_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * fpg0_63[k]
                   - f_6 * fpg1_63[k]
                   + f_4 * pc_z[k] * fph_91[k];

        t_125[k] = f_1 * dph_30[k]
                   + f_1 * fsh_30[k]
                   + f_4 * pc_y[k] * fph_93[k];

        t_126[k] = f_9 * fpg0_65[k]
                   - f_10 * fpg1_65[k]
                   + f_4 * pc_z[k] * fph_93[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t ppi1,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 * gamma / (p * q);
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

    const auto *ppi0_133 = buffer.data(ppi0 + 133);

    const auto *ppi1_133 = buffer.data(ppi1 + 133);

    const auto *dpi0_0 = buffer.data(dpi0 + 0);
    const auto *dpi0_3 = buffer.data(dpi0 + 3);
    const auto *dpi0_5 = buffer.data(dpi0 + 5);
    const auto *dpi0_6 = buffer.data(dpi0 + 6);
    const auto *dpi0_9 = buffer.data(dpi0 + 9);
    const auto *dpi0_10 = buffer.data(dpi0 + 10);
    const auto *dpi0_14 = buffer.data(dpi0 + 14);
    const auto *dpi0_15 = buffer.data(dpi0 + 15);
    const auto *dpi0_21 = buffer.data(dpi0 + 21);
    const auto *dpi0_28 = buffer.data(dpi0 + 28);
    const auto *dpi0_31 = buffer.data(dpi0 + 31);
    const auto *dpi0_34 = buffer.data(dpi0 + 34);
    const auto *dpi0_35 = buffer.data(dpi0 + 35);
    const auto *dpi0_38 = buffer.data(dpi0 + 38);
    const auto *dpi0_39 = buffer.data(dpi0 + 39);
    const auto *dpi0_40 = buffer.data(dpi0 + 40);
    const auto *dpi0_49 = buffer.data(dpi0 + 49);
    const auto *dpi0_56 = buffer.data(dpi0 + 56);
    const auto *dpi0_61 = buffer.data(dpi0 + 61);
    const auto *dpi0_65 = buffer.data(dpi0 + 65);
    const auto *dpi0_68 = buffer.data(dpi0 + 68);
    const auto *dpi0_70 = buffer.data(dpi0 + 70);
    const auto *dpi0_83 = buffer.data(dpi0 + 83);
    const auto *dpi0_133 = buffer.data(dpi0 + 133);

    const auto *dph_0 = buffer.data(dph + 0);
    const auto *dph_2 = buffer.data(dph + 2);
    const auto *dph_5 = buffer.data(dph + 5);
    const auto *dph_9 = buffer.data(dph + 9);
    const auto *dph_20 = buffer.data(dph + 20);
    const auto *dph_21 = buffer.data(dph + 21);
    const auto *dph_24 = buffer.data(dph + 24);
    const auto *dph_27 = buffer.data(dph + 27);
    const auto *dph_28 = buffer.data(dph + 28);
    const auto *dph_41 = buffer.data(dph + 41);
    const auto *dph_42 = buffer.data(dph + 42);
    const auto *dph_47 = buffer.data(dph + 47);
    const auto *dph_50 = buffer.data(dph + 50);
    const auto *dph_51 = buffer.data(dph + 51);
    const auto *dph_62 = buffer.data(dph + 62);
    const auto *dph_99 = buffer.data(dph + 99);
    const auto *dph_101 = buffer.data(dph + 101);
    const auto *dph_102 = buffer.data(dph + 102);
    const auto *dph_103 = buffer.data(dph + 103);
    const auto *dph_104 = buffer.data(dph + 104);
    const auto *dph_120 = buffer.data(dph + 120);
    const auto *dph_122 = buffer.data(dph + 122);
    const auto *dph_123 = buffer.data(dph + 123);
    const auto *dph_124 = buffer.data(dph + 124);
    const auto *dph_125 = buffer.data(dph + 125);
    const auto *dph_142 = buffer.data(dph + 142);
    const auto *dph_143 = buffer.data(dph + 143);
    const auto *dph_144 = buffer.data(dph + 144);
    const auto *dph_146 = buffer.data(dph + 146);
    const auto *dph_162 = buffer.data(dph + 162);
    const auto *dph_163 = buffer.data(dph + 163);
    const auto *dph_164 = buffer.data(dph + 164);
    const auto *dph_165 = buffer.data(dph + 165);
    const auto *dph_167 = buffer.data(dph + 167);
    const auto *dph_168 = buffer.data(dph + 168);
    const auto *dph_173 = buffer.data(dph + 173);
    const auto *dph_177 = buffer.data(dph + 177);
    const auto *dph_182 = buffer.data(dph + 182);
    const auto *dph_183 = buffer.data(dph + 183);
    const auto *dph_184 = buffer.data(dph + 184);
    const auto *dph_185 = buffer.data(dph + 185);
    const auto *dph_186 = buffer.data(dph + 186);
    const auto *dph_188 = buffer.data(dph + 188);

    const auto *dpi1_0 = buffer.data(dpi1 + 0);
    const auto *dpi1_3 = buffer.data(dpi1 + 3);
    const auto *dpi1_5 = buffer.data(dpi1 + 5);
    const auto *dpi1_6 = buffer.data(dpi1 + 6);
    const auto *dpi1_9 = buffer.data(dpi1 + 9);
    const auto *dpi1_10 = buffer.data(dpi1 + 10);
    const auto *dpi1_14 = buffer.data(dpi1 + 14);
    const auto *dpi1_15 = buffer.data(dpi1 + 15);
    const auto *dpi1_21 = buffer.data(dpi1 + 21);
    const auto *dpi1_28 = buffer.data(dpi1 + 28);
    const auto *dpi1_31 = buffer.data(dpi1 + 31);
    const auto *dpi1_34 = buffer.data(dpi1 + 34);
    const auto *dpi1_35 = buffer.data(dpi1 + 35);
    const auto *dpi1_38 = buffer.data(dpi1 + 38);
    const auto *dpi1_39 = buffer.data(dpi1 + 39);
    const auto *dpi1_40 = buffer.data(dpi1 + 40);
    const auto *dpi1_49 = buffer.data(dpi1 + 49);
    const auto *dpi1_56 = buffer.data(dpi1 + 56);
    const auto *dpi1_61 = buffer.data(dpi1 + 61);
    const auto *dpi1_65 = buffer.data(dpi1 + 65);
    const auto *dpi1_68 = buffer.data(dpi1 + 68);
    const auto *dpi1_70 = buffer.data(dpi1 + 70);
    const auto *dpi1_83 = buffer.data(dpi1 + 83);
    const auto *dpi1_133 = buffer.data(dpi1 + 133);

    const auto *fsi0_31 = buffer.data(fsi0 + 31);
    const auto *fsi0_34 = buffer.data(fsi0 + 34);
    const auto *fsi0_38 = buffer.data(fsi0 + 38);
    const auto *fsi0_49 = buffer.data(fsi0 + 49);
    const auto *fsi0_51 = buffer.data(fsi0 + 51);
    const auto *fsi0_52 = buffer.data(fsi0 + 52);
    const auto *fsi0_53 = buffer.data(fsi0 + 53);
    const auto *fsi0_61 = buffer.data(fsi0 + 61);
    const auto *fsi0_65 = buffer.data(fsi0 + 65);
    const auto *fsi0_70 = buffer.data(fsi0 + 70);
    const auto *fsi0_78 = buffer.data(fsi0 + 78);
    const auto *fsi0_79 = buffer.data(fsi0 + 79);
    const auto *fsi0_80 = buffer.data(fsi0 + 80);
    const auto *fsi0_81 = buffer.data(fsi0 + 81);
    const auto *fsi0_83 = buffer.data(fsi0 + 83);

    const auto *fsh_21 = buffer.data(fsh + 21);
    const auto *fsh_22 = buffer.data(fsh + 22);
    const auto *fsh_24 = buffer.data(fsh + 24);
    const auto *fsh_27 = buffer.data(fsh + 27);
    const auto *fsh_31 = buffer.data(fsh + 31);
    const auto *fsh_36 = buffer.data(fsh + 36);
    const auto *fsh_37 = buffer.data(fsh + 37);
    const auto *fsh_38 = buffer.data(fsh + 38);
    const auto *fsh_39 = buffer.data(fsh + 39);
    const auto *fsh_41 = buffer.data(fsh + 41);
    const auto *fsh_42 = buffer.data(fsh + 42);
    const auto *fsh_44 = buffer.data(fsh + 44);
    const auto *fsh_47 = buffer.data(fsh + 47);
    const auto *fsh_51 = buffer.data(fsh + 51);
    const auto *fsh_56 = buffer.data(fsh + 56);
    const auto *fsh_58 = buffer.data(fsh + 58);
    const auto *fsh_59 = buffer.data(fsh + 59);
    const auto *fsh_60 = buffer.data(fsh + 60);
    const auto *fsh_61 = buffer.data(fsh + 61);
    const auto *fsh_62 = buffer.data(fsh + 62);

    const auto *fsi1_31 = buffer.data(fsi1 + 31);
    const auto *fsi1_34 = buffer.data(fsi1 + 34);
    const auto *fsi1_38 = buffer.data(fsi1 + 38);
    const auto *fsi1_49 = buffer.data(fsi1 + 49);
    const auto *fsi1_51 = buffer.data(fsi1 + 51);
    const auto *fsi1_52 = buffer.data(fsi1 + 52);
    const auto *fsi1_53 = buffer.data(fsi1 + 53);
    const auto *fsi1_61 = buffer.data(fsi1 + 61);
    const auto *fsi1_65 = buffer.data(fsi1 + 65);
    const auto *fsi1_70 = buffer.data(fsi1 + 70);
    const auto *fsi1_78 = buffer.data(fsi1 + 78);
    const auto *fsi1_79 = buffer.data(fsi1 + 79);
    const auto *fsi1_80 = buffer.data(fsi1 + 80);
    const auto *fsi1_81 = buffer.data(fsi1 + 81);
    const auto *fsi1_83 = buffer.data(fsi1 + 83);

    const auto *fpg0_70 = buffer.data(fpg0 + 70);
    const auto *fpg0_71 = buffer.data(fpg0 + 71);
    const auto *fpg0_72 = buffer.data(fpg0 + 72);
    const auto *fpg0_74 = buffer.data(fpg0 + 74);
    const auto *fpg0_92 = buffer.data(fpg0 + 92);
    const auto *fpg0_94 = buffer.data(fpg0 + 94);
    const auto *fpg0_95 = buffer.data(fpg0 + 95);
    const auto *fpg0_101 = buffer.data(fpg0 + 101);
    const auto *fpg0_102 = buffer.data(fpg0 + 102);
    const auto *fpg0_103 = buffer.data(fpg0 + 103);
    const auto *fpg0_104 = buffer.data(fpg0 + 104);
    const auto *fpg0_120 = buffer.data(fpg0 + 120);
    const auto *fpg0_121 = buffer.data(fpg0 + 121);
    const auto *fpg0_122 = buffer.data(fpg0 + 122);
    const auto *fpg0_123 = buffer.data(fpg0 + 123);
    const auto *fpg0_124 = buffer.data(fpg0 + 124);
    const auto *fpg0_125 = buffer.data(fpg0 + 125);
    const auto *fpg0_129 = buffer.data(fpg0 + 129);
    const auto *fpg0_130 = buffer.data(fpg0 + 130);
    const auto *fpg0_131 = buffer.data(fpg0 + 131);
    const auto *fpg0_132 = buffer.data(fpg0 + 132);
    const auto *fpg0_134 = buffer.data(fpg0 + 134);

    const auto *fpg1_70 = buffer.data(fpg1 + 70);
    const auto *fpg1_71 = buffer.data(fpg1 + 71);
    const auto *fpg1_72 = buffer.data(fpg1 + 72);
    const auto *fpg1_74 = buffer.data(fpg1 + 74);
    const auto *fpg1_92 = buffer.data(fpg1 + 92);
    const auto *fpg1_94 = buffer.data(fpg1 + 94);
    const auto *fpg1_95 = buffer.data(fpg1 + 95);
    const auto *fpg1_101 = buffer.data(fpg1 + 101);
    const auto *fpg1_102 = buffer.data(fpg1 + 102);
    const auto *fpg1_103 = buffer.data(fpg1 + 103);
    const auto *fpg1_104 = buffer.data(fpg1 + 104);
    const auto *fpg1_120 = buffer.data(fpg1 + 120);
    const auto *fpg1_121 = buffer.data(fpg1 + 121);
    const auto *fpg1_122 = buffer.data(fpg1 + 122);
    const auto *fpg1_123 = buffer.data(fpg1 + 123);
    const auto *fpg1_124 = buffer.data(fpg1 + 124);
    const auto *fpg1_125 = buffer.data(fpg1 + 125);
    const auto *fpg1_129 = buffer.data(fpg1 + 129);
    const auto *fpg1_130 = buffer.data(fpg1 + 130);
    const auto *fpg1_131 = buffer.data(fpg1 + 131);
    const auto *fpg1_132 = buffer.data(fpg1 + 132);
    const auto *fpg1_134 = buffer.data(fpg1 + 134);

    const auto *fph_94 = buffer.data(fph + 94);
    const auto *fph_99 = buffer.data(fph + 99);
    const auto *fph_100 = buffer.data(fph + 100);
    const auto *fph_101 = buffer.data(fph + 101);
    const auto *fph_102 = buffer.data(fph + 102);
    const auto *fph_103 = buffer.data(fph + 103);
    const auto *fph_104 = buffer.data(fph + 104);
    const auto *fph_105 = buffer.data(fph + 105);
    const auto *fph_106 = buffer.data(fph + 106);
    const auto *fph_108 = buffer.data(fph + 108);
    const auto *fph_110 = buffer.data(fph + 110);
    const auto *fph_111 = buffer.data(fph + 111);
    const auto *fph_114 = buffer.data(fph + 114);
    const auto *fph_115 = buffer.data(fph + 115);
    const auto *fph_120 = buffer.data(fph + 120);
    const auto *fph_122 = buffer.data(fph + 122);
    const auto *fph_123 = buffer.data(fph + 123);
    const auto *fph_124 = buffer.data(fph + 124);
    const auto *fph_125 = buffer.data(fph + 125);
    const auto *fph_126 = buffer.data(fph + 126);
    const auto *fph_128 = buffer.data(fph + 128);
    const auto *fph_130 = buffer.data(fph + 130);
    const auto *fph_131 = buffer.data(fph + 131);
    const auto *fph_133 = buffer.data(fph + 133);
    const auto *fph_134 = buffer.data(fph + 134);
    const auto *fph_135 = buffer.data(fph + 135);
    const auto *fph_140 = buffer.data(fph + 140);
    const auto *fph_142 = buffer.data(fph + 142);
    const auto *fph_143 = buffer.data(fph + 143);
    const auto *fph_144 = buffer.data(fph + 144);
    const auto *fph_145 = buffer.data(fph + 145);
    const auto *fph_146 = buffer.data(fph + 146);
    const auto *fph_147 = buffer.data(fph + 147);
    const auto *fph_149 = buffer.data(fph + 149);
    const auto *fph_152 = buffer.data(fph + 152);
    const auto *fph_156 = buffer.data(fph + 156);
    const auto *fph_161 = buffer.data(fph + 161);
    const auto *fph_162 = buffer.data(fph + 162);
    const auto *fph_163 = buffer.data(fph + 163);
    const auto *fph_164 = buffer.data(fph + 164);
    const auto *fph_165 = buffer.data(fph + 165);
    const auto *fph_167 = buffer.data(fph + 167);
    const auto *fph_168 = buffer.data(fph + 168);
    const auto *fph_169 = buffer.data(fph + 169);
    const auto *fph_170 = buffer.data(fph + 170);
    const auto *fph_171 = buffer.data(fph + 171);
    const auto *fph_172 = buffer.data(fph + 172);
    const auto *fph_173 = buffer.data(fph + 173);
    const auto *fph_174 = buffer.data(fph + 174);
    const auto *fph_175 = buffer.data(fph + 175);
    const auto *fph_176 = buffer.data(fph + 176);
    const auto *fph_177 = buffer.data(fph + 177);
    const auto *fph_182 = buffer.data(fph + 182);
    const auto *fph_183 = buffer.data(fph + 183);
    const auto *fph_184 = buffer.data(fph + 184);
    const auto *fph_185 = buffer.data(fph + 185);
    const auto *fph_186 = buffer.data(fph + 186);
    const auto *fph_188 = buffer.data(fph + 188);

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, pc_z, dph_99, dph_101, \
                         dph_102, dph_103, fph_94, fph_99, fph_101, fph_102, \
                         fph_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_12 * dph_99[k]
                   + f_4 * pc_x[k] * fph_99[k];

        t_128[k] = f_4 * pc_z[k] * fph_94[k];

        t_129[k] = f_12 * dph_101[k]
                   + f_4 * pc_x[k] * fph_101[k];

        t_130[k] = f_12 * dph_102[k]
                   + f_4 * pc_x[k] * fph_102[k];

        t_131[k] = f_12 * dph_103[k]
                   + f_4 * pc_x[k] * fph_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pa_x, pc_x, pc_z, ppi0_133, ppi1_133, dpi0_133, \
                         dph_104, dpi1_133, fph_99, fph_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_12 * dph_104[k]
                   + f_4 * pc_x[k] * fph_104[k];

        t_133[k] = f_15 * ppi0_133[k]
                   - f_16 * ppi1_133[k]
                   + pa_x[k] * dpi0_133[k]
                   - f_11 * pc_x[k] * dpi1_133[k];

        t_134[k] = f_4 * pc_z[k] * fph_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pc_z, fpg0_70, fpg0_71, fpg0_72, fpg1_70, \
                         fpg1_71, fpg1_72, fph_100, fph_101, fph_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_5 * fpg0_70[k]
                   - f_6 * fpg1_70[k]
                   + f_4 * pc_z[k] * fph_100[k];

        t_136[k] = f_7 * fpg0_71[k]
                   - f_8 * fpg1_71[k]
                   + f_4 * pc_z[k] * fph_101[k];

        t_137[k] = f_9 * fpg0_72[k]
                   - f_10 * fpg1_72[k]
                   + f_4 * pc_z[k] * fph_102[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, t_141, pa_y, pc_y, pc_z, dpi0_56, dph_41, \
                         dph_42, dpi1_56, fsh_41, fpg0_74, fpg1_74, fph_104, \
                         fph_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_1 * dph_41[k]
                   + f_1 * fsh_41[k]
                   + f_4 * pc_y[k] * fph_104[k];

        t_139[k] = f_2 * fpg0_74[k]
                   - f_3 * fpg1_74[k]
                   + f_4 * pc_z[k] * fph_104[k];

        t_140[k] = pa_y[k] * dpi0_56[k]
                   - f_11 * pc_y[k] * dpi1_56[k];

        t_141[k] = f_1 * dph_42[k]
                   + f_4 * pc_y[k] * fph_105[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, t_145, pa_y, pb_z, pc_y, pc_z, dpi0_61, dpi1_61, \
                         fsi0_31, fsh_21, fsh_22, fsi1_31, fph_105, \
                         fph_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_1 * fsh_21[k]
                   + f_4 * pc_z[k] * fph_105[k];

        t_143[k] = pb_z[k] * fsi0_31[k]
                   - f_11 * pc_z[k] * fsi1_31[k];

        t_144[k] = f_1 * fsh_22[k]
                   + f_4 * pc_z[k] * fph_106[k];

        t_145[k] = pa_y[k] * dpi0_61[k]
                   - f_11 * pc_y[k] * dpi1_61[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_y, pb_z, pc_y, pc_z, dpi0_65, dph_47, \
                         dpi1_65, fsi0_34, fsh_24, fsi1_34, fph_108, \
                         fph_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pb_z[k] * fsi0_34[k]
                   - f_11 * pc_z[k] * fsi1_34[k];

        t_147[k] = f_1 * fsh_24[k]
                   + f_4 * pc_z[k] * fph_108[k];

        t_148[k] = f_1 * dph_47[k]
                   + f_4 * pc_y[k] * fph_110[k];

        t_149[k] = pa_y[k] * dpi0_65[k]
                   - f_11 * pc_y[k] * dpi1_65[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_y, pb_z, pc_y, pc_z, dpi0_68, dph_50, \
                         dpi1_68, fsi0_38, fsh_27, fsi1_38, fph_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = pb_z[k] * fsi0_38[k]
                   - f_11 * pc_z[k] * fsi1_38[k];

        t_151[k] = f_1 * fsh_27[k]
                   + f_4 * pc_z[k] * fph_111[k];

        t_152[k] = pa_y[k] * dpi0_68[k]
                   + f_12 * dph_50[k]
                   - f_11 * pc_y[k] * dpi1_68[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_y, pc_x, pc_y, pc_z, dpi0_70, dph_51, \
                         dph_120, dpi1_70, fsh_31, fph_114, fph_115, \
                         fph_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_1 * dph_51[k]
                   + f_4 * pc_y[k] * fph_114[k];

        t_154[k] = pa_y[k] * dpi0_70[k]
                   - f_11 * pc_y[k] * dpi1_70[k];

        t_155[k] = f_12 * dph_120[k]
                   + f_4 * pc_x[k] * fph_120[k];

        t_156[k] = f_1 * fsh_31[k]
                   + f_4 * pc_z[k] * fph_115[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, pc_x, dph_122, dph_123, dph_124, dph_125, \
                         fph_122, fph_123, fph_124, fph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_12 * dph_122[k]
                   + f_4 * pc_x[k] * fph_122[k];

        t_158[k] = f_12 * dph_123[k]
                   + f_4 * pc_x[k] * fph_123[k];

        t_159[k] = f_12 * dph_124[k]
                   + f_4 * pc_x[k] * fph_124[k];

        t_160[k] = f_12 * dph_125[k]
                   + f_4 * pc_x[k] * fph_125[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, t_164, pb_z, pc_z, fsi0_49, fsi0_51, fsi0_52, \
                         fsh_36, fsh_37, fsh_38, fsi1_49, fsi1_51, fsi1_52, \
                         fph_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = pb_z[k] * fsi0_49[k]
                   - f_11 * pc_z[k] * fsi1_49[k];

        t_162[k] = f_1 * fsh_36[k]
                   + f_4 * pc_z[k] * fph_120[k];

        t_163[k] = pb_z[k] * fsi0_51[k]
                   + f_12 * fsh_37[k]
                   - f_11 * pc_z[k] * fsi1_51[k];

        t_164[k] = pb_z[k] * fsi0_52[k]
                   + f_0 * fsh_38[k]
                   - f_11 * pc_z[k] * fsi1_52[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pa_y, pb_z, pc_y, pc_z, dpi0_83, dph_62, \
                         dpi1_83, fsi0_53, fsh_39, fsi1_53, fph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = pb_z[k] * fsi0_53[k]
                   + f_13 * fsh_39[k]
                   - f_11 * pc_z[k] * fsi1_53[k];

        t_166[k] = f_1 * dph_62[k]
                   + f_4 * pc_y[k] * fph_125[k];

        t_167[k] = pa_y[k] * dpi0_83[k]
                   - f_11 * pc_y[k] * dpi1_83[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, t_171, t_172, pa_z, pc_y, pc_z, dpi0_0, dpi0_3, \
                         dph_0, dpi1_0, dpi1_3, fph_126, fph_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_z[k] * dpi0_0[k]
                   - f_11 * pc_z[k] * dpi1_0[k];

        t_169[k] = f_4 * pc_y[k] * fph_126[k];

        t_170[k] = f_1 * dph_0[k]
                   + f_4 * pc_z[k] * fph_126[k];

        t_171[k] = pa_z[k] * dpi0_3[k]
                   - f_11 * pc_z[k] * dpi1_3[k];

        t_172[k] = f_4 * pc_y[k] * fph_128[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pc_y, pc_z, dpi0_5, dpi0_6, dph_2, \
                         dpi1_5, dpi1_6, fpg0_92, fpg1_92, fph_130, \
                         fph_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = pa_z[k] * dpi0_5[k]
                   + f_12 * dph_2[k]
                   - f_11 * pc_z[k] * dpi1_5[k];

        t_174[k] = pa_z[k] * dpi0_6[k]
                   - f_11 * pc_z[k] * dpi1_6[k];

        t_175[k] = f_5 * fpg0_92[k]
                   - f_6 * fpg1_92[k]
                   + f_4 * pc_y[k] * fph_130[k];

        t_176[k] = f_4 * pc_y[k] * fph_131[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_z, pc_y, pc_z, dpi0_9, dpi0_10, dph_5, \
                         dpi1_9, dpi1_10, fpg0_94, fpg1_94, fph_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_z[k] * dpi0_9[k]
                   + f_0 * dph_5[k]
                   - f_11 * pc_z[k] * dpi1_9[k];

        t_178[k] = pa_z[k] * dpi0_10[k]
                   - f_11 * pc_z[k] * dpi1_10[k];

        t_179[k] = f_7 * fpg0_94[k]
                   - f_8 * fpg1_94[k]
                   + f_4 * pc_y[k] * fph_133[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, t_183, pa_z, pc_y, pc_z, dpi0_14, dpi0_15, \
                         dph_9, dpi1_14, dpi1_15, fpg0_95, fpg1_95, fph_134, \
                         fph_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_5 * fpg0_95[k]
                   - f_6 * fpg1_95[k]
                   + f_4 * pc_y[k] * fph_134[k];

        t_181[k] = f_4 * pc_y[k] * fph_135[k];

        t_182[k] = pa_z[k] * dpi0_14[k]
                   + f_13 * dph_9[k]
                   - f_11 * pc_z[k] * dpi1_14[k];

        t_183[k] = pa_z[k] * dpi0_15[k]
                   - f_11 * pc_z[k] * dpi1_15[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, t_187, pc_x, pc_y, dph_142, dph_143, dph_144, \
                         fsh_58, fsh_59, fsh_60, fph_140, fph_142, fph_143, \
                         fph_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_12 * dph_142[k]
                   + f_1 * fsh_58[k]
                   + f_4 * pc_x[k] * fph_142[k];

        t_185[k] = f_12 * dph_143[k]
                   + f_1 * fsh_59[k]
                   + f_4 * pc_x[k] * fph_143[k];

        t_186[k] = f_12 * dph_144[k]
                   + f_1 * fsh_60[k]
                   + f_4 * pc_x[k] * fph_144[k];

        t_187[k] = f_4 * pc_y[k] * fph_140[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_z, pc_x, pc_y, pc_z, dpi0_21, dph_146, \
                         dpi1_21, fsh_62, fpg0_101, fpg1_101, fph_142, \
                         fph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_12 * dph_146[k]
                   + f_1 * fsh_62[k]
                   + f_4 * pc_x[k] * fph_146[k];

        t_189[k] = pa_z[k] * dpi0_21[k]
                   - f_11 * pc_z[k] * dpi1_21[k];

        t_190[k] = f_17 * fpg0_101[k]
                   - f_18 * fpg1_101[k]
                   + f_4 * pc_y[k] * fph_142[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, t_194, pc_y, fpg0_102, fpg0_103, fpg0_104, \
                         fpg1_102, fpg1_103, fpg1_104, fph_143, fph_144, fph_145, \
                         fph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_9 * fpg0_102[k]
                   - f_10 * fpg1_102[k]
                   + f_4 * pc_y[k] * fph_143[k];

        t_192[k] = f_7 * fpg0_103[k]
                   - f_8 * fpg1_103[k]
                   + f_4 * pc_y[k] * fph_144[k];

        t_193[k] = f_5 * fpg0_104[k]
                   - f_6 * fpg1_104[k]
                   + f_4 * pc_y[k] * fph_145[k];

        t_194[k] = f_4 * pc_y[k] * fph_146[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, t_198, pa_z, pc_y, pc_z, dpi0_28, dph_20, \
                         dph_21, dpi1_28, fsh_42, fpg0_104, fpg1_104, fph_146, \
                         fph_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_1 * dph_20[k]
                   + f_2 * fpg0_104[k]
                   - f_3 * fpg1_104[k]
                   + f_4 * pc_z[k] * fph_146[k];

        t_196[k] = pa_z[k] * dpi0_28[k]
                   - f_11 * pc_z[k] * dpi1_28[k];

        t_197[k] = f_1 * fsh_42[k]
                   + f_4 * pc_y[k] * fph_147[k];

        t_198[k] = f_1 * dph_21[k]
                   + f_4 * pc_z[k] * fph_147[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, t_202, pa_z, pb_y, pc_y, pc_z, dpi0_31, dpi0_34, \
                         dpi1_31, dpi1_34, fsi0_61, fsh_44, fsi1_61, \
                         fph_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = pa_z[k] * dpi0_31[k]
                   - f_11 * pc_z[k] * dpi1_31[k];

        t_200[k] = f_1 * fsh_44[k]
                   + f_4 * pc_y[k] * fph_149[k];

        t_201[k] = pb_y[k] * fsi0_61[k]
                   - f_11 * pc_y[k] * fsi1_61[k];

        t_202[k] = pa_z[k] * dpi0_34[k]
                   - f_11 * pc_z[k] * dpi1_34[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pa_z, pb_y, pc_y, pc_z, dpi0_35, dph_24, \
                         dpi1_35, fsi0_65, fsh_47, fsi1_65, fph_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pa_z[k] * dpi0_35[k]
                   + f_1 * dph_24[k]
                   - f_11 * pc_z[k] * dpi1_35[k];

        t_204[k] = f_1 * fsh_47[k]
                   + f_4 * pc_y[k] * fph_152[k];

        t_205[k] = pb_y[k] * fsi0_65[k]
                   - f_11 * pc_y[k] * fsi1_65[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_z, pc_z, dpi0_38, dpi0_39, dpi0_40, dph_27, \
                         dph_28, dpi1_38, dpi1_39, dpi1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_z[k] * dpi0_38[k]
                   - f_11 * pc_z[k] * dpi1_38[k];

        t_207[k] = pa_z[k] * dpi0_39[k]
                   + f_1 * dph_27[k]
                   - f_11 * pc_z[k] * dpi1_39[k];

        t_208[k] = pa_z[k] * dpi0_40[k]
                   + f_12 * dph_28[k]
                   - f_11 * pc_z[k] * dpi1_40[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_y, pc_x, pc_y, dph_162, dph_163, \
                         fsi0_70, fsh_51, fsi1_70, fph_156, fph_162, \
                         fph_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_1 * fsh_51[k]
                   + f_4 * pc_y[k] * fph_156[k];

        t_210[k] = pb_y[k] * fsi0_70[k]
                   - f_11 * pc_y[k] * fsi1_70[k];

        t_211[k] = f_12 * dph_162[k]
                   + f_4 * pc_x[k] * fph_162[k];

        t_212[k] = f_12 * dph_163[k]
                   + f_4 * pc_x[k] * fph_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, dph_164, dph_165, dph_167, \
                         fsh_56, fph_161, fph_164, fph_165, fph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_12 * dph_164[k]
                   + f_4 * pc_x[k] * fph_164[k];

        t_214[k] = f_12 * dph_165[k]
                   + f_4 * pc_x[k] * fph_165[k];

        t_215[k] = f_1 * fsh_56[k]
                   + f_4 * pc_y[k] * fph_161[k];

        t_216[k] = f_12 * dph_167[k]
                   + f_4 * pc_x[k] * fph_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, pa_z, pb_y, pc_y, pc_z, dpi0_49, dpi1_49, \
                         fsi0_78, fsi0_79, fsh_58, fsh_59, fsi1_78, \
                         fsi1_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_z[k] * dpi0_49[k]
                   - f_11 * pc_z[k] * dpi1_49[k];

        t_218[k] = pb_y[k] * fsi0_78[k]
                   + f_19 * fsh_58[k]
                   - f_11 * pc_y[k] * fsi1_78[k];

        t_219[k] = pb_y[k] * fsi0_79[k]
                   + f_13 * fsh_59[k]
                   - f_11 * pc_y[k] * fsi1_79[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pb_y, pc_y, fsi0_80, fsi0_81, fsi0_83, \
                         fsh_60, fsh_61, fsh_62, fsi1_80, fsi1_81, fsi1_83, \
                         fph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pb_y[k] * fsi0_80[k]
                   + f_0 * fsh_60[k]
                   - f_11 * pc_y[k] * fsi1_80[k];

        t_221[k] = pb_y[k] * fsi0_81[k]
                   + f_12 * fsh_61[k]
                   - f_11 * pc_y[k] * fsi1_81[k];

        t_222[k] = f_1 * fsh_62[k]
                   + f_4 * pc_y[k] * fph_167[k];

        t_223[k] = pb_y[k] * fsi0_83[k]
                   - f_11 * pc_y[k] * fsi1_83[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pc_x, pc_y, pc_z, dph_42, dph_168, \
                         fsh_42, fpg0_120, fpg1_120, fph_168, fph_169, \
                         fph_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_12 * dph_168[k]
                   + f_2 * fpg0_120[k]
                   - f_3 * fpg1_120[k]
                   + f_4 * pc_x[k] * fph_168[k];

        t_225[k] = f_4 * pc_y[k] * fph_168[k];

        t_226[k] = f_1 * dph_42[k]
                   + f_1 * fsh_42[k]
                   + f_4 * pc_z[k] * fph_168[k];

        t_227[k] = f_5 * fpg0_120[k]
                   - f_6 * fpg1_120[k]
                   + f_4 * pc_y[k] * fph_169[k];

        t_228[k] = f_4 * pc_y[k] * fph_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pc_x, pc_y, dph_173, fpg0_121, fpg0_122, \
                         fpg0_125, fpg1_121, fpg1_122, fpg1_125, fph_171, fph_172, \
                         fph_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_12 * dph_173[k]
                   + f_9 * fpg0_125[k]
                   - f_10 * fpg1_125[k]
                   + f_4 * pc_x[k] * fph_173[k];

        t_230[k] = f_7 * fpg0_121[k]
                   - f_8 * fpg1_121[k]
                   + f_4 * pc_y[k] * fph_171[k];

        t_231[k] = f_5 * fpg0_122[k]
                   - f_6 * fpg1_122[k]
                   + f_4 * pc_y[k] * fph_172[k];

        t_232[k] = f_4 * pc_y[k] * fph_173[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, pc_x, pc_y, dph_177, fpg0_123, fpg0_124, \
                         fpg0_129, fpg1_123, fpg1_124, fpg1_129, fph_174, fph_175, \
                         fph_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_12 * dph_177[k]
                   + f_7 * fpg0_129[k]
                   - f_8 * fpg1_129[k]
                   + f_4 * pc_x[k] * fph_177[k];

        t_234[k] = f_9 * fpg0_123[k]
                   - f_10 * fpg1_123[k]
                   + f_4 * pc_y[k] * fph_174[k];

        t_235[k] = f_7 * fpg0_124[k]
                   - f_8 * fpg1_124[k]
                   + f_4 * pc_y[k] * fph_175[k];
    }

#pragma omp simd aligned(t_236, t_237, t_238, t_239, pc_x, pc_y, dph_182, dph_183, fpg0_125, \
                         fpg0_134, fpg1_125, fpg1_134, fph_176, fph_177, fph_182, \
                         fph_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_236[k] = f_5 * fpg0_125[k]
                   - f_6 * fpg1_125[k]
                   + f_4 * pc_y[k] * fph_176[k];

        t_237[k] = f_4 * pc_y[k] * fph_177[k];

        t_238[k] = f_12 * dph_182[k]
                   + f_5 * fpg0_134[k]
                   - f_6 * fpg1_134[k]
                   + f_4 * pc_x[k] * fph_182[k];

        t_239[k] = f_12 * dph_183[k]
                   + f_4 * pc_x[k] * fph_183[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_y, dph_184, dph_185, \
                         dph_186, dph_188, fph_182, fph_184, fph_185, fph_186, \
                         fph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_12 * dph_184[k]
                   + f_4 * pc_x[k] * fph_184[k];

        t_241[k] = f_12 * dph_185[k]
                   + f_4 * pc_x[k] * fph_185[k];

        t_242[k] = f_12 * dph_186[k]
                   + f_4 * pc_x[k] * fph_186[k];

        t_243[k] = f_4 * pc_y[k] * fph_182[k];

        t_244[k] = f_12 * dph_188[k]
                   + f_4 * pc_x[k] * fph_188[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pc_y, fpg0_130, fpg0_131, fpg0_132, fpg1_130, \
                         fpg1_131, fpg1_132, fph_183, fph_184, \
                         fph_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_2 * fpg0_130[k]
                   - f_3 * fpg1_130[k]
                   + f_4 * pc_y[k] * fph_183[k];

        t_246[k] = f_17 * fpg0_131[k]
                   - f_18 * fpg1_131[k]
                   + f_4 * pc_y[k] * fph_184[k];

        t_247[k] = f_9 * fpg0_132[k]
                   - f_10 * fpg1_132[k]
                   + f_4 * pc_y[k] * fph_185[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t ppi1,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *dpi0_87 = buffer.data(dpi0 + 87);
    const auto *dpi0_90 = buffer.data(dpi0 + 90);
    const auto *dpi0_94 = buffer.data(dpi0 + 94);
    const auto *dpi0_99 = buffer.data(dpi0 + 99);
    const auto *dpi0_113 = buffer.data(dpi0 + 113);
    const auto *dpi0_168 = buffer.data(dpi0 + 168);
    const auto *dpi0_173 = buffer.data(dpi0 + 173);
    const auto *dpi0_177 = buffer.data(dpi0 + 177);
    const auto *dpi0_182 = buffer.data(dpi0 + 182);
    const auto *dpi0_188 = buffer.data(dpi0 + 188);
    const auto *dpi0_251 = buffer.data(dpi0 + 251);
    const auto *dpi0_280 = buffer.data(dpi0 + 280);
    const auto *dpi0_283 = buffer.data(dpi0 + 283);
    const auto *dpi0_286 = buffer.data(dpi0 + 286);
    const auto *dpi0_290 = buffer.data(dpi0 + 290);
    const auto *dpi0_301 = buffer.data(dpi0 + 301);
    const auto *dpi0_303 = buffer.data(dpi0 + 303);
    const auto *dpi0_304 = buffer.data(dpi0 + 304);
    const auto *dpi0_305 = buffer.data(dpi0 + 305);
    const auto *dpi0_306 = buffer.data(dpi0 + 306);
    const auto *dpi0_307 = buffer.data(dpi0 + 307);
    const auto *dpi0_313 = buffer.data(dpi0 + 313);
    const auto *dpi0_317 = buffer.data(dpi0 + 317);
    const auto *dpi0_320 = buffer.data(dpi0 + 320);
    const auto *dpi0_322 = buffer.data(dpi0 + 322);
    const auto *dpi0_329 = buffer.data(dpi0 + 329);
    const auto *dpi0_331 = buffer.data(dpi0 + 331);
    const auto *dpi0_332 = buffer.data(dpi0 + 332);
    const auto *dpi0_333 = buffer.data(dpi0 + 333);
    const auto *dpi0_335 = buffer.data(dpi0 + 335);

    const auto *dph_63 = buffer.data(dph + 63);
    const auto *dph_66 = buffer.data(dph + 66);
    const auto *dph_68 = buffer.data(dph + 68);
    const auto *dph_69 = buffer.data(dph + 69);
    const auto *dph_72 = buffer.data(dph + 72);
    const auto *dph_78 = buffer.data(dph + 78);
    const auto *dph_83 = buffer.data(dph + 83);
    const auto *dph_84 = buffer.data(dph + 84);
    const auto *dph_89 = buffer.data(dph + 89);
    const auto *dph_93 = buffer.data(dph + 93);
    const auto *dph_105 = buffer.data(dph + 105);
    const auto *dph_110 = buffer.data(dph + 110);
    const auto *dph_114 = buffer.data(dph + 114);
    const auto *dph_125 = buffer.data(dph + 125);
    const auto *dph_126 = buffer.data(dph + 126);
    const auto *dph_128 = buffer.data(dph + 128);
    const auto *dph_131 = buffer.data(dph + 131);
    const auto *dph_135 = buffer.data(dph + 135);
    const auto *dph_141 = buffer.data(dph + 141);
    const auto *dph_143 = buffer.data(dph + 143);
    const auto *dph_144 = buffer.data(dph + 144);
    const auto *dph_145 = buffer.data(dph + 145);
    const auto *dph_146 = buffer.data(dph + 146);
    const auto *dph_189 = buffer.data(dph + 189);
    const auto *dph_192 = buffer.data(dph + 192);
    const auto *dph_195 = buffer.data(dph + 195);
    const auto *dph_199 = buffer.data(dph + 199);
    const auto *dph_204 = buffer.data(dph + 204);
    const auto *dph_206 = buffer.data(dph + 206);
    const auto *dph_207 = buffer.data(dph + 207);
    const auto *dph_208 = buffer.data(dph + 208);
    const auto *dph_209 = buffer.data(dph + 209);
    const auto *dph_210 = buffer.data(dph + 210);
    const auto *dph_213 = buffer.data(dph + 213);
    const auto *dph_216 = buffer.data(dph + 216);
    const auto *dph_220 = buffer.data(dph + 220);
    const auto *dph_225 = buffer.data(dph + 225);
    const auto *dph_227 = buffer.data(dph + 227);
    const auto *dph_228 = buffer.data(dph + 228);
    const auto *dph_229 = buffer.data(dph + 229);
    const auto *dph_230 = buffer.data(dph + 230);
    const auto *dph_236 = buffer.data(dph + 236);
    const auto *dph_240 = buffer.data(dph + 240);
    const auto *dph_243 = buffer.data(dph + 243);
    const auto *dph_245 = buffer.data(dph + 245);
    const auto *dph_246 = buffer.data(dph + 246);
    const auto *dph_248 = buffer.data(dph + 248);
    const auto *dph_249 = buffer.data(dph + 249);
    const auto *dph_250 = buffer.data(dph + 250);
    const auto *dph_251 = buffer.data(dph + 251);
    const auto *dph_264 = buffer.data(dph + 264);
    const auto *dph_268 = buffer.data(dph + 268);
    const auto *dph_269 = buffer.data(dph + 269);
    const auto *dph_270 = buffer.data(dph + 270);
    const auto *dph_271 = buffer.data(dph + 271);
    const auto *dph_273 = buffer.data(dph + 273);

    const auto *dpi1_87 = buffer.data(dpi1 + 87);
    const auto *dpi1_90 = buffer.data(dpi1 + 90);
    const auto *dpi1_94 = buffer.data(dpi1 + 94);
    const auto *dpi1_99 = buffer.data(dpi1 + 99);
    const auto *dpi1_113 = buffer.data(dpi1 + 113);
    const auto *dpi1_168 = buffer.data(dpi1 + 168);
    const auto *dpi1_173 = buffer.data(dpi1 + 173);
    const auto *dpi1_177 = buffer.data(dpi1 + 177);
    const auto *dpi1_182 = buffer.data(dpi1 + 182);
    const auto *dpi1_188 = buffer.data(dpi1 + 188);
    const auto *dpi1_251 = buffer.data(dpi1 + 251);
    const auto *dpi1_280 = buffer.data(dpi1 + 280);
    const auto *dpi1_283 = buffer.data(dpi1 + 283);
    const auto *dpi1_286 = buffer.data(dpi1 + 286);
    const auto *dpi1_290 = buffer.data(dpi1 + 290);
    const auto *dpi1_301 = buffer.data(dpi1 + 301);
    const auto *dpi1_303 = buffer.data(dpi1 + 303);
    const auto *dpi1_304 = buffer.data(dpi1 + 304);
    const auto *dpi1_305 = buffer.data(dpi1 + 305);
    const auto *dpi1_306 = buffer.data(dpi1 + 306);
    const auto *dpi1_307 = buffer.data(dpi1 + 307);
    const auto *dpi1_313 = buffer.data(dpi1 + 313);
    const auto *dpi1_317 = buffer.data(dpi1 + 317);
    const auto *dpi1_320 = buffer.data(dpi1 + 320);
    const auto *dpi1_322 = buffer.data(dpi1 + 322);
    const auto *dpi1_329 = buffer.data(dpi1 + 329);
    const auto *dpi1_331 = buffer.data(dpi1 + 331);
    const auto *dpi1_332 = buffer.data(dpi1 + 332);
    const auto *dpi1_333 = buffer.data(dpi1 + 333);
    const auto *dpi1_335 = buffer.data(dpi1 + 335);

    const auto *fsi0_84 = buffer.data(fsi0 + 84);
    const auto *fsi0_87 = buffer.data(fsi0 + 87);
    const auto *fsi0_90 = buffer.data(fsi0 + 90);
    const auto *fsi0_94 = buffer.data(fsi0 + 94);

    const auto *fsh_63 = buffer.data(fsh + 63);
    const auto *fsh_64 = buffer.data(fsh + 64);
    const auto *fsh_66 = buffer.data(fsh + 66);
    const auto *fsh_68 = buffer.data(fsh + 68);
    const auto *fsh_69 = buffer.data(fsh + 69);
    const auto *fsh_72 = buffer.data(fsh + 72);
    const auto *fsh_73 = buffer.data(fsh + 73);
    const auto *fsh_78 = buffer.data(fsh + 78);
    const auto *fsh_80 = buffer.data(fsh + 80);
    const auto *fsh_81 = buffer.data(fsh + 81);
    const auto *fsh_82 = buffer.data(fsh + 82);
    const auto *fsh_83 = buffer.data(fsh + 83);
    const auto *fsh_96 = buffer.data(fsh + 96);
    const auto *fsh_100 = buffer.data(fsh + 100);
    const auto *fsh_101 = buffer.data(fsh + 101);
    const auto *fsh_102 = buffer.data(fsh + 102);
    const auto *fsh_103 = buffer.data(fsh + 103);

    const auto *fsi1_84 = buffer.data(fsi1 + 84);
    const auto *fsi1_87 = buffer.data(fsi1 + 87);
    const auto *fsi1_90 = buffer.data(fsi1 + 90);
    const auto *fsi1_94 = buffer.data(fsi1 + 94);

    const auto *fpg0_133 = buffer.data(fpg0 + 133);
    const auto *fpg0_134 = buffer.data(fpg0 + 134);
    const auto *fpg0_135 = buffer.data(fpg0 + 135);
    const auto *fpg0_137 = buffer.data(fpg0 + 137);
    const auto *fpg0_138 = buffer.data(fpg0 + 138);
    const auto *fpg0_140 = buffer.data(fpg0 + 140);
    const auto *fpg0_141 = buffer.data(fpg0 + 141);
    const auto *fpg0_145 = buffer.data(fpg0 + 145);
    const auto *fpg0_146 = buffer.data(fpg0 + 146);
    const auto *fpg0_147 = buffer.data(fpg0 + 147);
    const auto *fpg0_149 = buffer.data(fpg0 + 149);
    const auto *fpg0_150 = buffer.data(fpg0 + 150);
    const auto *fpg0_152 = buffer.data(fpg0 + 152);
    const auto *fpg0_153 = buffer.data(fpg0 + 153);
    const auto *fpg0_155 = buffer.data(fpg0 + 155);
    const auto *fpg0_190 = buffer.data(fpg0 + 190);
    const auto *fpg0_192 = buffer.data(fpg0 + 192);
    const auto *fpg0_193 = buffer.data(fpg0 + 193);
    const auto *fpg0_194 = buffer.data(fpg0 + 194);
    const auto *fpg0_195 = buffer.data(fpg0 + 195);

    const auto *fpg1_133 = buffer.data(fpg1 + 133);
    const auto *fpg1_134 = buffer.data(fpg1 + 134);
    const auto *fpg1_135 = buffer.data(fpg1 + 135);
    const auto *fpg1_137 = buffer.data(fpg1 + 137);
    const auto *fpg1_138 = buffer.data(fpg1 + 138);
    const auto *fpg1_140 = buffer.data(fpg1 + 140);
    const auto *fpg1_141 = buffer.data(fpg1 + 141);
    const auto *fpg1_145 = buffer.data(fpg1 + 145);
    const auto *fpg1_146 = buffer.data(fpg1 + 146);
    const auto *fpg1_147 = buffer.data(fpg1 + 147);
    const auto *fpg1_149 = buffer.data(fpg1 + 149);
    const auto *fpg1_150 = buffer.data(fpg1 + 150);
    const auto *fpg1_152 = buffer.data(fpg1 + 152);
    const auto *fpg1_153 = buffer.data(fpg1 + 153);
    const auto *fpg1_155 = buffer.data(fpg1 + 155);
    const auto *fpg1_190 = buffer.data(fpg1 + 190);
    const auto *fpg1_192 = buffer.data(fpg1 + 192);
    const auto *fpg1_193 = buffer.data(fpg1 + 193);
    const auto *fpg1_194 = buffer.data(fpg1 + 194);
    const auto *fpg1_195 = buffer.data(fpg1 + 195);

    const auto *fph_186 = buffer.data(fph + 186);
    const auto *fph_187 = buffer.data(fph + 187);
    const auto *fph_188 = buffer.data(fph + 188);
    const auto *fph_189 = buffer.data(fph + 189);
    const auto *fph_190 = buffer.data(fph + 190);
    const auto *fph_191 = buffer.data(fph + 191);
    const auto *fph_192 = buffer.data(fph + 192);
    const auto *fph_194 = buffer.data(fph + 194);
    const auto *fph_195 = buffer.data(fph + 195);
    const auto *fph_196 = buffer.data(fph + 196);
    const auto *fph_198 = buffer.data(fph + 198);
    const auto *fph_199 = buffer.data(fph + 199);
    const auto *fph_204 = buffer.data(fph + 204);
    const auto *fph_205 = buffer.data(fph + 205);
    const auto *fph_206 = buffer.data(fph + 206);
    const auto *fph_207 = buffer.data(fph + 207);
    const auto *fph_208 = buffer.data(fph + 208);
    const auto *fph_209 = buffer.data(fph + 209);
    const auto *fph_210 = buffer.data(fph + 210);
    const auto *fph_211 = buffer.data(fph + 211);
    const auto *fph_212 = buffer.data(fph + 212);
    const auto *fph_213 = buffer.data(fph + 213);
    const auto *fph_215 = buffer.data(fph + 215);
    const auto *fph_216 = buffer.data(fph + 216);
    const auto *fph_217 = buffer.data(fph + 217);
    const auto *fph_219 = buffer.data(fph + 219);
    const auto *fph_220 = buffer.data(fph + 220);
    const auto *fph_225 = buffer.data(fph + 225);
    const auto *fph_227 = buffer.data(fph + 227);
    const auto *fph_228 = buffer.data(fph + 228);
    const auto *fph_229 = buffer.data(fph + 229);
    const auto *fph_230 = buffer.data(fph + 230);
    const auto *fph_231 = buffer.data(fph + 231);
    const auto *fph_232 = buffer.data(fph + 232);
    const auto *fph_234 = buffer.data(fph + 234);
    const auto *fph_236 = buffer.data(fph + 236);
    const auto *fph_237 = buffer.data(fph + 237);
    const auto *fph_240 = buffer.data(fph + 240);
    const auto *fph_241 = buffer.data(fph + 241);
    const auto *fph_246 = buffer.data(fph + 246);
    const auto *fph_248 = buffer.data(fph + 248);
    const auto *fph_249 = buffer.data(fph + 249);
    const auto *fph_250 = buffer.data(fph + 250);
    const auto *fph_251 = buffer.data(fph + 251);
    const auto *fph_252 = buffer.data(fph + 252);
    const auto *fph_254 = buffer.data(fph + 254);
    const auto *fph_255 = buffer.data(fph + 255);
    const auto *fph_257 = buffer.data(fph + 257);
    const auto *fph_258 = buffer.data(fph + 258);
    const auto *fph_261 = buffer.data(fph + 261);
    const auto *fph_264 = buffer.data(fph + 264);
    const auto *fph_267 = buffer.data(fph + 267);
    const auto *fph_268 = buffer.data(fph + 268);
    const auto *fph_269 = buffer.data(fph + 269);
    const auto *fph_270 = buffer.data(fph + 270);
    const auto *fph_271 = buffer.data(fph + 271);
    const auto *fph_272 = buffer.data(fph + 272);
    const auto *fph_273 = buffer.data(fph + 273);

#pragma omp simd aligned(t_248, t_249, t_250, pc_y, fpg0_133, fpg0_134, fpg1_133, fpg1_134, \
                         fph_186, fph_187, fph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_7 * fpg0_133[k]
                   - f_8 * fpg1_133[k]
                   + f_4 * pc_y[k] * fph_186[k];

        t_249[k] = f_5 * fpg0_134[k]
                   - f_6 * fpg1_134[k]
                   + f_4 * pc_y[k] * fph_187[k];

        t_250[k] = f_4 * pc_y[k] * fph_188[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pa_x, pc_x, pc_y, ppi0_251, ppi1_251, dpi0_251, \
                         dph_63, dph_189, dpi1_251, fsh_63, fpg0_135, fpg1_135, \
                         fph_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_15 * ppi0_251[k]
                   - f_16 * ppi1_251[k]
                   + pa_x[k] * dpi0_251[k]
                   - f_11 * pc_x[k] * dpi1_251[k];

        t_252[k] = f_1 * dph_189[k]
                   + f_1 * fsh_63[k]
                   + f_2 * fpg0_135[k]
                   - f_3 * fpg1_135[k]
                   + f_4 * pc_x[k] * fph_189[k];

        t_253[k] = f_12 * dph_63[k]
                   + f_4 * pc_y[k] * fph_189[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pc_x, pc_z, dph_192, fsh_66, fpg0_135, \
                         fpg0_138, fpg1_135, fpg1_138, fph_189, fph_190, fph_191, \
                         fph_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_4 * pc_z[k] * fph_189[k];

        t_255[k] = f_1 * dph_192[k]
                   + f_1 * fsh_66[k]
                   + f_9 * fpg0_138[k]
                   - f_10 * fpg1_138[k]
                   + f_4 * pc_x[k] * fph_192[k];

        t_256[k] = f_4 * pc_z[k] * fph_190[k];

        t_257[k] = f_5 * fpg0_135[k]
                   - f_6 * fpg1_135[k]
                   + f_4 * pc_z[k] * fph_191[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pc_x, pc_y, pc_z, dph_68, dph_195, fsh_69, \
                         fpg0_141, fpg1_141, fph_192, fph_194, \
                         fph_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_1 * dph_195[k]
                   + f_1 * fsh_69[k]
                   + f_7 * fpg0_141[k]
                   - f_8 * fpg1_141[k]
                   + f_4 * pc_x[k] * fph_195[k];

        t_259[k] = f_4 * pc_z[k] * fph_192[k];

        t_260[k] = f_12 * dph_68[k]
                   + f_4 * pc_y[k] * fph_194[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pc_x, pc_z, dph_199, fsh_73, fpg0_137, fpg0_145, \
                         fpg1_137, fpg1_145, fph_194, fph_195, \
                         fph_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_7 * fpg0_137[k]
                   - f_8 * fpg1_137[k]
                   + f_4 * pc_z[k] * fph_194[k];

        t_262[k] = f_1 * dph_199[k]
                   + f_1 * fsh_73[k]
                   + f_5 * fpg0_145[k]
                   - f_6 * fpg1_145[k]
                   + f_4 * pc_x[k] * fph_199[k];

        t_263[k] = f_4 * pc_z[k] * fph_195[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pc_y, pc_z, dph_72, fpg0_138, fpg0_140, \
                         fpg1_138, fpg1_140, fph_196, fph_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_5 * fpg0_138[k]
                   - f_6 * fpg1_138[k]
                   + f_4 * pc_z[k] * fph_196[k];

        t_265[k] = f_12 * dph_72[k]
                   + f_4 * pc_y[k] * fph_198[k];

        t_266[k] = f_9 * fpg0_140[k]
                   - f_10 * fpg1_140[k]
                   + f_4 * pc_z[k] * fph_198[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, pc_x, pc_z, dph_204, dph_206, dph_207, \
                         fsh_78, fsh_80, fsh_81, fph_199, fph_204, fph_206, \
                         fph_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_1 * dph_204[k]
                   + f_1 * fsh_78[k]
                   + f_4 * pc_x[k] * fph_204[k];

        t_268[k] = f_4 * pc_z[k] * fph_199[k];

        t_269[k] = f_1 * dph_206[k]
                   + f_1 * fsh_80[k]
                   + f_4 * pc_x[k] * fph_206[k];

        t_270[k] = f_1 * dph_207[k]
                   + f_1 * fsh_81[k]
                   + f_4 * pc_x[k] * fph_207[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pc_x, pc_y, dph_78, dph_208, dph_209, fsh_82, \
                         fsh_83, fpg0_145, fpg1_145, fph_204, fph_208, \
                         fph_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = f_1 * dph_208[k]
                   + f_1 * fsh_82[k]
                   + f_4 * pc_x[k] * fph_208[k];

        t_272[k] = f_1 * dph_209[k]
                   + f_1 * fsh_83[k]
                   + f_4 * pc_x[k] * fph_209[k];

        t_273[k] = f_12 * dph_78[k]
                   + f_2 * fpg0_145[k]
                   - f_3 * fpg1_145[k]
                   + f_4 * pc_y[k] * fph_204[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pc_z, fpg0_145, fpg0_146, fpg0_147, \
                         fpg1_145, fpg1_146, fpg1_147, fph_204, fph_205, fph_206, \
                         fph_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * pc_z[k] * fph_204[k];

        t_275[k] = f_5 * fpg0_145[k]
                   - f_6 * fpg1_145[k]
                   + f_4 * pc_z[k] * fph_205[k];

        t_276[k] = f_7 * fpg0_146[k]
                   - f_8 * fpg1_146[k]
                   + f_4 * pc_z[k] * fph_206[k];

        t_277[k] = f_9 * fpg0_147[k]
                   - f_10 * fpg1_147[k]
                   + f_4 * pc_z[k] * fph_207[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pa_x, pc_x, pc_y, pc_z, dpi0_280, dph_83, \
                         dph_210, dpi1_280, fpg0_149, fpg1_149, \
                         fph_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_12 * dph_83[k]
                   + f_4 * pc_y[k] * fph_209[k];

        t_279[k] = f_2 * fpg0_149[k]
                   - f_3 * fpg1_149[k]
                   + f_4 * pc_z[k] * fph_209[k];

        t_280[k] = pa_x[k] * dpi0_280[k]
                   + f_14 * dph_210[k]
                   - f_11 * pc_x[k] * dpi1_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pa_x, pc_x, pc_y, pc_z, dpi0_283, dph_84, \
                         dph_213, dpi1_283, fsh_63, fph_210, fph_211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_12 * dph_84[k]
                   + f_1 * fsh_63[k]
                   + f_4 * pc_y[k] * fph_210[k];

        t_282[k] = f_4 * pc_z[k] * fph_210[k];

        t_283[k] = pa_x[k] * dpi0_283[k]
                   + f_13 * dph_213[k]
                   - f_11 * pc_x[k] * dpi1_283[k];

        t_284[k] = f_4 * pc_z[k] * fph_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_x, pc_x, pc_z, dpi0_286, dph_216, dpi1_286, \
                         fpg0_150, fpg1_150, fph_212, fph_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_5 * fpg0_150[k]
                   - f_6 * fpg1_150[k]
                   + f_4 * pc_z[k] * fph_212[k];

        t_286[k] = pa_x[k] * dpi0_286[k]
                   + f_0 * dph_216[k]
                   - f_11 * pc_x[k] * dpi1_286[k];

        t_287[k] = f_4 * pc_z[k] * fph_213[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_x, pc_x, pc_y, pc_z, dpi0_290, dph_89, \
                         dph_220, dpi1_290, fsh_68, fpg0_152, fpg1_152, \
                         fph_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_12 * dph_89[k]
                   + f_1 * fsh_68[k]
                   + f_4 * pc_y[k] * fph_215[k];

        t_289[k] = f_7 * fpg0_152[k]
                   - f_8 * fpg1_152[k]
                   + f_4 * pc_z[k] * fph_215[k];

        t_290[k] = pa_x[k] * dpi0_290[k]
                   + f_12 * dph_220[k]
                   - f_11 * pc_x[k] * dpi1_290[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, t_294, pc_y, pc_z, dph_93, fsh_72, fpg0_153, \
                         fpg0_155, fpg1_153, fpg1_155, fph_216, fph_217, \
                         fph_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_4 * pc_z[k] * fph_216[k];

        t_292[k] = f_5 * fpg0_153[k]
                   - f_6 * fpg1_153[k]
                   + f_4 * pc_z[k] * fph_217[k];

        t_293[k] = f_12 * dph_93[k]
                   + f_1 * fsh_72[k]
                   + f_4 * pc_y[k] * fph_219[k];

        t_294[k] = f_9 * fpg0_155[k]
                   - f_10 * fpg1_155[k]
                   + f_4 * pc_z[k] * fph_219[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, t_298, t_299, pc_x, pc_z, dph_225, dph_227, \
                         dph_228, dph_229, fph_220, fph_225, fph_227, fph_228, \
                         fph_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = f_1 * dph_225[k]
                   + f_4 * pc_x[k] * fph_225[k];

        t_296[k] = f_4 * pc_z[k] * fph_220[k];

        t_297[k] = f_1 * dph_227[k]
                   + f_4 * pc_x[k] * fph_227[k];

        t_298[k] = f_1 * dph_228[k]
                   + f_4 * pc_x[k] * fph_228[k];

        t_299[k] = f_1 * dph_229[k]
                   + f_4 * pc_x[k] * fph_229[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, t_303, pa_x, pc_x, pc_z, dpi0_301, dpi0_303, \
                         dph_230, dpi1_301, dpi1_303, fph_225, \
                         fph_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_1 * dph_230[k]
                   + f_4 * pc_x[k] * fph_230[k];

        t_301[k] = pa_x[k] * dpi0_301[k]
                   - f_11 * pc_x[k] * dpi1_301[k];

        t_302[k] = f_4 * pc_z[k] * fph_225[k];

        t_303[k] = pa_x[k] * dpi0_303[k]
                   - f_11 * pc_x[k] * dpi1_303[k];
    }

#pragma omp simd aligned(t_304, t_305, t_306, t_307, pa_x, pc_x, dpi0_304, dpi0_305, dpi0_306, \
                         dpi0_307, dpi1_304, dpi1_305, dpi1_306, \
                         dpi1_307 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_304[k] = pa_x[k] * dpi0_304[k]
                   - f_11 * pc_x[k] * dpi1_304[k];

        t_305[k] = pa_x[k] * dpi0_305[k]
                   - f_11 * pc_x[k] * dpi1_305[k];

        t_306[k] = pa_x[k] * dpi0_306[k]
                   - f_11 * pc_x[k] * dpi1_306[k];

        t_307[k] = pa_x[k] * dpi0_307[k]
                   - f_11 * pc_x[k] * dpi1_307[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pb_z, pc_y, pc_z, dph_105, fsi0_84, \
                         fsi0_87, fsh_63, fsi1_84, fsi1_87, fph_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = pb_z[k] * fsi0_84[k]
                   - f_11 * pc_z[k] * fsi1_84[k];

        t_309[k] = f_12 * dph_105[k]
                   + f_4 * pc_y[k] * fph_231[k];

        t_310[k] = f_1 * fsh_63[k]
                   + f_4 * pc_z[k] * fph_231[k];

        t_311[k] = pb_z[k] * fsi0_87[k]
                   - f_11 * pc_z[k] * fsi1_87[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pb_z, pc_x, pc_z, dpi0_313, dph_236, \
                         dpi1_313, fsi0_90, fsh_64, fsi1_90, fph_232 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = f_1 * fsh_64[k]
                   + f_4 * pc_z[k] * fph_232[k];

        t_313[k] = pa_x[k] * dpi0_313[k]
                   + f_13 * dph_236[k]
                   - f_11 * pc_x[k] * dpi1_313[k];

        t_314[k] = pb_z[k] * fsi0_90[k]
                   - f_11 * pc_z[k] * fsi1_90[k];
    }

#pragma omp simd aligned(t_315, t_316, t_317, pa_x, pc_x, pc_y, pc_z, dpi0_317, dph_110, \
                         dph_240, dpi1_317, fsh_66, fph_234, fph_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = f_1 * fsh_66[k]
                   + f_4 * pc_z[k] * fph_234[k];

        t_316[k] = f_12 * dph_110[k]
                   + f_4 * pc_y[k] * fph_236[k];

        t_317[k] = pa_x[k] * dpi0_317[k]
                   + f_0 * dph_240[k]
                   - f_11 * pc_x[k] * dpi1_317[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pa_x, pb_z, pc_x, pc_z, dpi0_320, dph_243, \
                         dpi1_320, fsi0_94, fsh_69, fsi1_94, fph_237 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pb_z[k] * fsi0_94[k]
                   - f_11 * pc_z[k] * fsi1_94[k];

        t_319[k] = f_1 * fsh_69[k]
                   + f_4 * pc_z[k] * fph_237[k];

        t_320[k] = pa_x[k] * dpi0_320[k]
                   + f_12 * dph_243[k]
                   - f_11 * pc_x[k] * dpi1_320[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, pa_x, pc_x, pc_y, dpi0_322, dph_114, dph_245, \
                         dph_246, dpi1_322, fph_240, fph_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_12 * dph_114[k]
                   + f_4 * pc_y[k] * fph_240[k];

        t_322[k] = pa_x[k] * dpi0_322[k]
                   + f_12 * dph_245[k]
                   - f_11 * pc_x[k] * dpi1_322[k];

        t_323[k] = f_1 * dph_246[k]
                   + f_4 * pc_x[k] * fph_246[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, pc_x, pc_z, dph_248, dph_249, dph_250, \
                         fsh_73, fph_241, fph_248, fph_249, fph_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * fsh_73[k]
                   + f_4 * pc_z[k] * fph_241[k];

        t_325[k] = f_1 * dph_248[k]
                   + f_4 * pc_x[k] * fph_248[k];

        t_326[k] = f_1 * dph_249[k]
                   + f_4 * pc_x[k] * fph_249[k];

        t_327[k] = f_1 * dph_250[k]
                   + f_4 * pc_x[k] * fph_250[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pa_x, pc_x, pc_z, dpi0_329, dpi0_331, \
                         dph_251, dpi1_329, dpi1_331, fsh_78, fph_246, \
                         fph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_1 * dph_251[k]
                   + f_4 * pc_x[k] * fph_251[k];

        t_329[k] = pa_x[k] * dpi0_329[k]
                   - f_11 * pc_x[k] * dpi1_329[k];

        t_330[k] = f_1 * fsh_78[k]
                   + f_4 * pc_z[k] * fph_246[k];

        t_331[k] = pa_x[k] * dpi0_331[k]
                   - f_11 * pc_x[k] * dpi1_331[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, t_335, pa_x, pc_x, pc_y, dpi0_332, dpi0_333, \
                         dpi0_335, dph_125, dpi1_332, dpi1_333, dpi1_335, \
                         fph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = pa_x[k] * dpi0_332[k]
                   - f_11 * pc_x[k] * dpi1_332[k];

        t_333[k] = pa_x[k] * dpi0_333[k]
                   - f_11 * pc_x[k] * dpi1_333[k];

        t_334[k] = f_12 * dph_125[k]
                   + f_4 * pc_y[k] * fph_251[k];

        t_335[k] = pa_x[k] * dpi0_335[k]
                   - f_11 * pc_x[k] * dpi1_335[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_y, pa_z, pc_y, pc_z, dpi0_87, \
                         dpi0_168, dph_63, dph_126, dpi1_87, dpi1_168, \
                         fph_252 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_y[k] * dpi0_168[k]
                   - f_11 * pc_y[k] * dpi1_168[k];

        t_337[k] = f_1 * dph_126[k]
                   + f_4 * pc_y[k] * fph_252[k];

        t_338[k] = f_1 * dph_63[k]
                   + f_4 * pc_z[k] * fph_252[k];

        t_339[k] = pa_z[k] * dpi0_87[k]
                   - f_11 * pc_z[k] * dpi1_87[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, t_343, pa_y, pa_z, pc_y, pc_z, dpi0_90, \
                         dpi0_173, dph_66, dph_128, dpi1_90, dpi1_173, fph_254, \
                         fph_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_1 * dph_128[k]
                   + f_4 * pc_y[k] * fph_254[k];

        t_341[k] = pa_y[k] * dpi0_173[k]
                   - f_11 * pc_y[k] * dpi1_173[k];

        t_342[k] = pa_z[k] * dpi0_90[k]
                   - f_11 * pc_z[k] * dpi1_90[k];

        t_343[k] = f_1 * dph_66[k]
                   + f_4 * pc_z[k] * fph_255[k];
    }

#pragma omp simd aligned(t_344, t_345, t_346, t_347, pa_y, pa_z, pc_y, pc_z, dpi0_94, \
                         dpi0_177, dph_69, dph_131, dpi1_94, dpi1_177, fph_257, \
                         fph_258 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = f_1 * dph_131[k]
                   + f_4 * pc_y[k] * fph_257[k];

        t_345[k] = pa_y[k] * dpi0_177[k]
                   - f_11 * pc_y[k] * dpi1_177[k];

        t_346[k] = pa_z[k] * dpi0_94[k]
                   - f_11 * pc_z[k] * dpi1_94[k];

        t_347[k] = f_1 * dph_69[k]
                   + f_4 * pc_z[k] * fph_258[k];
    }

#pragma omp simd aligned(t_348, t_349, t_350, pa_y, pc_x, pc_y, dpi0_182, dph_135, dph_264, \
                         dpi1_182, fsh_96, fpg0_192, fpg1_192, fph_261, \
                         fph_264 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_348[k] = f_1 * dph_264[k]
                   + f_1 * fsh_96[k]
                   + f_5 * fpg0_192[k]
                   - f_6 * fpg1_192[k]
                   + f_4 * pc_x[k] * fph_264[k];

        t_349[k] = f_1 * dph_135[k]
                   + f_4 * pc_y[k] * fph_261[k];

        t_350[k] = pa_y[k] * dpi0_182[k]
                   - f_11 * pc_y[k] * dpi1_182[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pa_z, pc_x, pc_z, dpi0_99, dph_268, dph_269, \
                         dpi1_99, fsh_100, fsh_101, fph_268, fph_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = pa_z[k] * dpi0_99[k]
                   - f_11 * pc_z[k] * dpi1_99[k];

        t_352[k] = f_1 * dph_268[k]
                   + f_1 * fsh_100[k]
                   + f_4 * pc_x[k] * fph_268[k];

        t_353[k] = f_1 * dph_269[k]
                   + f_1 * fsh_101[k]
                   + f_4 * pc_x[k] * fph_269[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_y, pc_x, pc_y, dpi0_188, dph_270, dph_271, \
                         dpi1_188, fsh_102, fsh_103, fph_270, fph_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_1 * dph_270[k]
                   + f_1 * fsh_102[k]
                   + f_4 * pc_x[k] * fph_270[k];

        t_355[k] = f_1 * dph_271[k]
                   + f_1 * fsh_103[k]
                   + f_4 * pc_x[k] * fph_271[k];

        t_356[k] = pa_y[k] * dpi0_188[k]
                   - f_11 * pc_y[k] * dpi1_188[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pc_y, pc_z, dph_78, dph_141, dph_143, fpg0_190, \
                         fpg0_192, fpg1_190, fpg1_192, fph_267, \
                         fph_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = f_1 * dph_141[k]
                   + f_2 * fpg0_190[k]
                   - f_3 * fpg1_190[k]
                   + f_4 * pc_y[k] * fph_267[k];

        t_358[k] = f_1 * dph_78[k]
                   + f_4 * pc_z[k] * fph_267[k];

        t_359[k] = f_1 * dph_143[k]
                   + f_9 * fpg0_192[k]
                   - f_10 * fpg1_192[k]
                   + f_4 * pc_y[k] * fph_269[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pc_y, dph_144, dph_145, dph_146, fpg0_193, \
                         fpg0_194, fpg1_193, fpg1_194, fph_270, fph_271, \
                         fph_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = f_1 * dph_144[k]
                   + f_7 * fpg0_193[k]
                   - f_8 * fpg1_193[k]
                   + f_4 * pc_y[k] * fph_270[k];

        t_361[k] = f_1 * dph_145[k]
                   + f_5 * fpg0_194[k]
                   - f_6 * fpg1_194[k]
                   + f_4 * pc_y[k] * fph_271[k];

        t_362[k] = f_1 * dph_146[k]
                   + f_4 * pc_y[k] * fph_272[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pa_z, pc_x, pc_z, dpi0_113, dph_83, dph_273, \
                         dpi1_113, fpg0_194, fpg0_195, fpg1_194, fpg1_195, fph_272, \
                         fph_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_1 * dph_83[k]
                   + f_2 * fpg0_194[k]
                   - f_3 * fpg1_194[k]
                   + f_4 * pc_z[k] * fph_272[k];

        t_364[k] = f_1 * dph_273[k]
                   + f_2 * fpg0_195[k]
                   - f_3 * fpg1_195[k]
                   + f_4 * pc_x[k] * fph_273[k];

        t_365[k] = pa_z[k] * dpi0_113[k]
                   - f_11 * pc_z[k] * dpi1_113[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_115 = buffer.data(dpi0 + 115);
    const auto *dpi0_118 = buffer.data(dpi0 + 118);
    const auto *dpi0_122 = buffer.data(dpi0 + 122);
    const auto *dpi0_224 = buffer.data(dpi0 + 224);
    const auto *dpi0_226 = buffer.data(dpi0 + 226);
    const auto *dpi0_229 = buffer.data(dpi0 + 229);
    const auto *dpi0_233 = buffer.data(dpi0 + 233);
    const auto *dpi0_238 = buffer.data(dpi0 + 238);
    const auto *dpi0_376 = buffer.data(dpi0 + 376);
    const auto *dpi0_385 = buffer.data(dpi0 + 385);
    const auto *dpi0_387 = buffer.data(dpi0 + 387);
    const auto *dpi0_388 = buffer.data(dpi0 + 388);
    const auto *dpi0_389 = buffer.data(dpi0 + 389);
    const auto *dpi0_390 = buffer.data(dpi0 + 390);
    const auto *dpi0_391 = buffer.data(dpi0 + 391);
    const auto *dpi0_404 = buffer.data(dpi0 + 404);
    const auto *dpi0_413 = buffer.data(dpi0 + 413);
    const auto *dpi0_414 = buffer.data(dpi0 + 414);
    const auto *dpi0_415 = buffer.data(dpi0 + 415);
    const auto *dpi0_416 = buffer.data(dpi0 + 416);
    const auto *dpi0_417 = buffer.data(dpi0 + 417);
    const auto *dpi0_419 = buffer.data(dpi0 + 419);
    const auto *dpi0_451 = buffer.data(dpi0 + 451);
    const auto *dpi0_454 = buffer.data(dpi0 + 454);
    const auto *dpi0_455 = buffer.data(dpi0 + 455);
    const auto *dpi0_458 = buffer.data(dpi0 + 458);
    const auto *dpi0_459 = buffer.data(dpi0 + 459);
    const auto *dpi0_460 = buffer.data(dpi0 + 460);
    const auto *dpi0_469 = buffer.data(dpi0 + 469);
    const auto *dpi0_470 = buffer.data(dpi0 + 470);
    const auto *dpi0_471 = buffer.data(dpi0 + 471);
    const auto *dpi0_472 = buffer.data(dpi0 + 472);
    const auto *dpi0_473 = buffer.data(dpi0 + 473);
    const auto *dpi0_475 = buffer.data(dpi0 + 475);
    const auto *dpi0_476 = buffer.data(dpi0 + 476);
    const auto *dpi0_481 = buffer.data(dpi0 + 481);

    const auto *dph_84 = buffer.data(dph + 84);
    const auto *dph_87 = buffer.data(dph + 87);
    const auto *dph_90 = buffer.data(dph + 90);
    const auto *dph_99 = buffer.data(dph + 99);
    const auto *dph_108 = buffer.data(dph + 108);
    const auto *dph_111 = buffer.data(dph + 111);
    const auto *dph_126 = buffer.data(dph + 126);
    const auto *dph_146 = buffer.data(dph + 146);
    const auto *dph_147 = buffer.data(dph + 147);
    const auto *dph_149 = buffer.data(dph + 149);
    const auto *dph_152 = buffer.data(dph + 152);
    const auto *dph_156 = buffer.data(dph + 156);
    const auto *dph_168 = buffer.data(dph + 168);
    const auto *dph_170 = buffer.data(dph + 170);
    const auto *dph_173 = buffer.data(dph + 173);
    const auto *dph_177 = buffer.data(dph + 177);
    const auto *dph_188 = buffer.data(dph + 188);
    const auto *dph_278 = buffer.data(dph + 278);
    const auto *dph_282 = buffer.data(dph + 282);
    const auto *dph_285 = buffer.data(dph + 285);
    const auto *dph_287 = buffer.data(dph + 287);
    const auto *dph_288 = buffer.data(dph + 288);
    const auto *dph_289 = buffer.data(dph + 289);
    const auto *dph_290 = buffer.data(dph + 290);
    const auto *dph_291 = buffer.data(dph + 291);
    const auto *dph_292 = buffer.data(dph + 292);
    const auto *dph_293 = buffer.data(dph + 293);
    const auto *dph_297 = buffer.data(dph + 297);
    const auto *dph_300 = buffer.data(dph + 300);
    const auto *dph_304 = buffer.data(dph + 304);
    const auto *dph_306 = buffer.data(dph + 306);
    const auto *dph_309 = buffer.data(dph + 309);
    const auto *dph_310 = buffer.data(dph + 310);
    const auto *dph_311 = buffer.data(dph + 311);
    const auto *dph_312 = buffer.data(dph + 312);
    const auto *dph_313 = buffer.data(dph + 313);
    const auto *dph_314 = buffer.data(dph + 314);
    const auto *dph_315 = buffer.data(dph + 315);
    const auto *dph_320 = buffer.data(dph + 320);
    const auto *dph_324 = buffer.data(dph + 324);
    const auto *dph_329 = buffer.data(dph + 329);
    const auto *dph_330 = buffer.data(dph + 330);
    const auto *dph_331 = buffer.data(dph + 331);
    const auto *dph_332 = buffer.data(dph + 332);
    const auto *dph_333 = buffer.data(dph + 333);
    const auto *dph_335 = buffer.data(dph + 335);
    const auto *dph_339 = buffer.data(dph + 339);
    const auto *dph_342 = buffer.data(dph + 342);
    const auto *dph_343 = buffer.data(dph + 343);
    const auto *dph_346 = buffer.data(dph + 346);
    const auto *dph_347 = buffer.data(dph + 347);
    const auto *dph_348 = buffer.data(dph + 348);
    const auto *dph_351 = buffer.data(dph + 351);
    const auto *dph_352 = buffer.data(dph + 352);
    const auto *dph_353 = buffer.data(dph + 353);
    const auto *dph_354 = buffer.data(dph + 354);
    const auto *dph_356 = buffer.data(dph + 356);
    const auto *dph_357 = buffer.data(dph + 357);
    const auto *dph_362 = buffer.data(dph + 362);

    const auto *dpi1_115 = buffer.data(dpi1 + 115);
    const auto *dpi1_118 = buffer.data(dpi1 + 118);
    const auto *dpi1_122 = buffer.data(dpi1 + 122);
    const auto *dpi1_224 = buffer.data(dpi1 + 224);
    const auto *dpi1_226 = buffer.data(dpi1 + 226);
    const auto *dpi1_229 = buffer.data(dpi1 + 229);
    const auto *dpi1_233 = buffer.data(dpi1 + 233);
    const auto *dpi1_238 = buffer.data(dpi1 + 238);
    const auto *dpi1_376 = buffer.data(dpi1 + 376);
    const auto *dpi1_385 = buffer.data(dpi1 + 385);
    const auto *dpi1_387 = buffer.data(dpi1 + 387);
    const auto *dpi1_388 = buffer.data(dpi1 + 388);
    const auto *dpi1_389 = buffer.data(dpi1 + 389);
    const auto *dpi1_390 = buffer.data(dpi1 + 390);
    const auto *dpi1_391 = buffer.data(dpi1 + 391);
    const auto *dpi1_404 = buffer.data(dpi1 + 404);
    const auto *dpi1_413 = buffer.data(dpi1 + 413);
    const auto *dpi1_414 = buffer.data(dpi1 + 414);
    const auto *dpi1_415 = buffer.data(dpi1 + 415);
    const auto *dpi1_416 = buffer.data(dpi1 + 416);
    const auto *dpi1_417 = buffer.data(dpi1 + 417);
    const auto *dpi1_419 = buffer.data(dpi1 + 419);
    const auto *dpi1_451 = buffer.data(dpi1 + 451);
    const auto *dpi1_454 = buffer.data(dpi1 + 454);
    const auto *dpi1_455 = buffer.data(dpi1 + 455);
    const auto *dpi1_458 = buffer.data(dpi1 + 458);
    const auto *dpi1_459 = buffer.data(dpi1 + 459);
    const auto *dpi1_460 = buffer.data(dpi1 + 460);
    const auto *dpi1_469 = buffer.data(dpi1 + 469);
    const auto *dpi1_470 = buffer.data(dpi1 + 470);
    const auto *dpi1_471 = buffer.data(dpi1 + 471);
    const auto *dpi1_472 = buffer.data(dpi1 + 472);
    const auto *dpi1_473 = buffer.data(dpi1 + 473);
    const auto *dpi1_475 = buffer.data(dpi1 + 475);
    const auto *dpi1_476 = buffer.data(dpi1 + 476);
    const auto *dpi1_481 = buffer.data(dpi1 + 481);

    const auto *fsi0_140 = buffer.data(fsi0 + 140);
    const auto *fsi0_145 = buffer.data(fsi0 + 145);
    const auto *fsi0_149 = buffer.data(fsi0 + 149);
    const auto *fsi0_154 = buffer.data(fsi0 + 154);

    const auto *fsh_86 = buffer.data(fsh + 86);
    const auto *fsh_87 = buffer.data(fsh + 87);
    const auto *fsh_89 = buffer.data(fsh + 89);
    const auto *fsh_90 = buffer.data(fsh + 90);
    const auto *fsh_93 = buffer.data(fsh + 93);
    const auto *fsh_105 = buffer.data(fsh + 105);
    const auto *fsh_107 = buffer.data(fsh + 107);
    const auto *fsh_110 = buffer.data(fsh + 110);
    const auto *fsh_114 = buffer.data(fsh + 114);
    const auto *fsh_119 = buffer.data(fsh + 119);
    const auto *fsh_120 = buffer.data(fsh + 120);
    const auto *fsh_121 = buffer.data(fsh + 121);
    const auto *fsh_122 = buffer.data(fsh + 122);
    const auto *fsh_123 = buffer.data(fsh + 123);
    const auto *fsh_125 = buffer.data(fsh + 125);

    const auto *fsi1_140 = buffer.data(fsi1 + 140);
    const auto *fsi1_145 = buffer.data(fsi1 + 145);
    const auto *fsi1_149 = buffer.data(fsi1 + 149);
    const auto *fsi1_154 = buffer.data(fsi1 + 154);

    const auto *fpg0_200 = buffer.data(fpg0 + 200);
    const auto *fpg0_204 = buffer.data(fpg0 + 204);
    const auto *fpg0_209 = buffer.data(fpg0 + 209);
    const auto *fpg0_213 = buffer.data(fpg0 + 213);
    const auto *fpg0_216 = buffer.data(fpg0 + 216);
    const auto *fpg0_220 = buffer.data(fpg0 + 220);
    const auto *fpg0_225 = buffer.data(fpg0 + 225);
    const auto *fpg0_226 = buffer.data(fpg0 + 226);
    const auto *fpg0_227 = buffer.data(fpg0 + 227);
    const auto *fpg0_228 = buffer.data(fpg0 + 228);
    const auto *fpg0_229 = buffer.data(fpg0 + 229);
    const auto *fpg0_230 = buffer.data(fpg0 + 230);
    const auto *fpg0_234 = buffer.data(fpg0 + 234);
    const auto *fpg0_235 = buffer.data(fpg0 + 235);
    const auto *fpg0_236 = buffer.data(fpg0 + 236);
    const auto *fpg0_237 = buffer.data(fpg0 + 237);
    const auto *fpg0_238 = buffer.data(fpg0 + 238);
    const auto *fpg0_239 = buffer.data(fpg0 + 239);
    const auto *fpg0_255 = buffer.data(fpg0 + 255);

    const auto *fpg1_200 = buffer.data(fpg1 + 200);
    const auto *fpg1_204 = buffer.data(fpg1 + 204);
    const auto *fpg1_209 = buffer.data(fpg1 + 209);
    const auto *fpg1_213 = buffer.data(fpg1 + 213);
    const auto *fpg1_216 = buffer.data(fpg1 + 216);
    const auto *fpg1_220 = buffer.data(fpg1 + 220);
    const auto *fpg1_225 = buffer.data(fpg1 + 225);
    const auto *fpg1_226 = buffer.data(fpg1 + 226);
    const auto *fpg1_227 = buffer.data(fpg1 + 227);
    const auto *fpg1_228 = buffer.data(fpg1 + 228);
    const auto *fpg1_229 = buffer.data(fpg1 + 229);
    const auto *fpg1_230 = buffer.data(fpg1 + 230);
    const auto *fpg1_234 = buffer.data(fpg1 + 234);
    const auto *fpg1_235 = buffer.data(fpg1 + 235);
    const auto *fpg1_236 = buffer.data(fpg1 + 236);
    const auto *fpg1_237 = buffer.data(fpg1 + 237);
    const auto *fpg1_238 = buffer.data(fpg1 + 238);
    const auto *fpg1_239 = buffer.data(fpg1 + 239);
    const auto *fpg1_255 = buffer.data(fpg1 + 255);

    const auto *fph_273 = buffer.data(fph + 273);
    const auto *fph_275 = buffer.data(fph + 275);
    const auto *fph_276 = buffer.data(fph + 276);
    const auto *fph_278 = buffer.data(fph + 278);
    const auto *fph_279 = buffer.data(fph + 279);
    const auto *fph_282 = buffer.data(fph + 282);
    const auto *fph_287 = buffer.data(fph + 287);
    const auto *fph_288 = buffer.data(fph + 288);
    const auto *fph_289 = buffer.data(fph + 289);
    const auto *fph_290 = buffer.data(fph + 290);
    const auto *fph_291 = buffer.data(fph + 291);
    const auto *fph_292 = buffer.data(fph + 292);
    const auto *fph_293 = buffer.data(fph + 293);
    const auto *fph_294 = buffer.data(fph + 294);
    const auto *fph_296 = buffer.data(fph + 296);
    const auto *fph_297 = buffer.data(fph + 297);
    const auto *fph_299 = buffer.data(fph + 299);
    const auto *fph_300 = buffer.data(fph + 300);
    const auto *fph_303 = buffer.data(fph + 303);
    const auto *fph_304 = buffer.data(fph + 304);
    const auto *fph_309 = buffer.data(fph + 309);
    const auto *fph_310 = buffer.data(fph + 310);
    const auto *fph_311 = buffer.data(fph + 311);
    const auto *fph_312 = buffer.data(fph + 312);
    const auto *fph_313 = buffer.data(fph + 313);
    const auto *fph_314 = buffer.data(fph + 314);
    const auto *fph_315 = buffer.data(fph + 315);
    const auto *fph_316 = buffer.data(fph + 316);
    const auto *fph_317 = buffer.data(fph + 317);
    const auto *fph_318 = buffer.data(fph + 318);
    const auto *fph_319 = buffer.data(fph + 319);
    const auto *fph_320 = buffer.data(fph + 320);
    const auto *fph_321 = buffer.data(fph + 321);
    const auto *fph_322 = buffer.data(fph + 322);
    const auto *fph_323 = buffer.data(fph + 323);
    const auto *fph_324 = buffer.data(fph + 324);
    const auto *fph_329 = buffer.data(fph + 329);
    const auto *fph_330 = buffer.data(fph + 330);
    const auto *fph_331 = buffer.data(fph + 331);
    const auto *fph_332 = buffer.data(fph + 332);
    const auto *fph_333 = buffer.data(fph + 333);
    const auto *fph_334 = buffer.data(fph + 334);
    const auto *fph_335 = buffer.data(fph + 335);
    const auto *fph_336 = buffer.data(fph + 336);
    const auto *fph_338 = buffer.data(fph + 338);
    const auto *fph_341 = buffer.data(fph + 341);
    const auto *fph_345 = buffer.data(fph + 345);
    const auto *fph_350 = buffer.data(fph + 350);
    const auto *fph_351 = buffer.data(fph + 351);
    const auto *fph_352 = buffer.data(fph + 352);
    const auto *fph_353 = buffer.data(fph + 353);
    const auto *fph_354 = buffer.data(fph + 354);
    const auto *fph_356 = buffer.data(fph + 356);
    const auto *fph_357 = buffer.data(fph + 357);
    const auto *fph_358 = buffer.data(fph + 358);
    const auto *fph_359 = buffer.data(fph + 359);

#pragma omp simd aligned(t_366, t_367, t_368, pa_z, pc_y, pc_z, dpi0_115, dph_84, dph_149, \
                         dpi1_115, fsh_86, fph_273, fph_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_1 * dph_84[k]
                   + f_4 * pc_z[k] * fph_273[k];

        t_367[k] = pa_z[k] * dpi0_115[k]
                   - f_11 * pc_z[k] * dpi1_115[k];

        t_368[k] = f_1 * dph_149[k]
                   + f_1 * fsh_86[k]
                   + f_4 * pc_y[k] * fph_275[k];
    }

#pragma omp simd aligned(t_369, t_370, t_371, pa_z, pc_x, pc_z, dpi0_118, dph_87, dph_278, \
                         dpi1_118, fpg0_200, fpg1_200, fph_276, \
                         fph_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_369[k] = f_1 * dph_278[k]
                   + f_9 * fpg0_200[k]
                   - f_10 * fpg1_200[k]
                   + f_4 * pc_x[k] * fph_278[k];

        t_370[k] = pa_z[k] * dpi0_118[k]
                   - f_11 * pc_z[k] * dpi1_118[k];

        t_371[k] = f_1 * dph_87[k]
                   + f_4 * pc_z[k] * fph_276[k];
    }

#pragma omp simd aligned(t_372, t_373, t_374, pa_z, pc_x, pc_y, pc_z, dpi0_122, dph_152, \
                         dph_282, dpi1_122, fsh_89, fpg0_204, fpg1_204, fph_278, \
                         fph_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_372[k] = f_1 * dph_152[k]
                   + f_1 * fsh_89[k]
                   + f_4 * pc_y[k] * fph_278[k];

        t_373[k] = f_1 * dph_282[k]
                   + f_7 * fpg0_204[k]
                   - f_8 * fpg1_204[k]
                   + f_4 * pc_x[k] * fph_282[k];

        t_374[k] = pa_z[k] * dpi0_122[k]
                   - f_11 * pc_z[k] * dpi1_122[k];
    }

#pragma omp simd aligned(t_375, t_376, t_377, pa_x, pc_x, pc_y, pc_z, dpi0_376, dph_90, \
                         dph_156, dph_285, dpi1_376, fsh_93, fph_279, \
                         fph_282 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_1 * dph_90[k]
                   + f_4 * pc_z[k] * fph_279[k];

        t_376[k] = pa_x[k] * dpi0_376[k]
                   + f_12 * dph_285[k]
                   - f_11 * pc_x[k] * dpi1_376[k];

        t_377[k] = f_1 * dph_156[k]
                   + f_1 * fsh_93[k]
                   + f_4 * pc_y[k] * fph_282[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pc_x, dph_287, dph_288, dph_289, dph_290, \
                         fpg0_209, fpg1_209, fph_287, fph_288, fph_289, \
                         fph_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = f_1 * dph_287[k]
                   + f_5 * fpg0_209[k]
                   - f_6 * fpg1_209[k]
                   + f_4 * pc_x[k] * fph_287[k];

        t_379[k] = f_1 * dph_288[k]
                   + f_4 * pc_x[k] * fph_288[k];

        t_380[k] = f_1 * dph_289[k]
                   + f_4 * pc_x[k] * fph_289[k];

        t_381[k] = f_1 * dph_290[k]
                   + f_4 * pc_x[k] * fph_290[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_x, pc_x, dpi0_385, dph_291, dph_292, \
                         dph_293, dpi1_385, fph_291, fph_292, fph_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_1 * dph_291[k]
                   + f_4 * pc_x[k] * fph_291[k];

        t_383[k] = f_1 * dph_292[k]
                   + f_4 * pc_x[k] * fph_292[k];

        t_384[k] = f_1 * dph_293[k]
                   + f_4 * pc_x[k] * fph_293[k];

        t_385[k] = pa_x[k] * dpi0_385[k]
                   - f_11 * pc_x[k] * dpi1_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pa_x, pc_x, pc_z, dpi0_387, dpi0_388, \
                         dpi0_389, dph_99, dpi1_387, dpi1_388, dpi1_389, \
                         fph_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_1 * dph_99[k]
                   + f_4 * pc_z[k] * fph_288[k];

        t_387[k] = pa_x[k] * dpi0_387[k]
                   - f_11 * pc_x[k] * dpi1_387[k];

        t_388[k] = pa_x[k] * dpi0_388[k]
                   - f_11 * pc_x[k] * dpi1_388[k];

        t_389[k] = pa_x[k] * dpi0_389[k]
                   - f_11 * pc_x[k] * dpi1_389[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, t_393, pa_x, pa_y, pc_x, pc_y, dpi0_224, \
                         dpi0_390, dpi0_391, dph_168, dpi1_224, dpi1_390, dpi1_391, \
                         fph_294 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = pa_x[k] * dpi0_390[k]
                   - f_11 * pc_x[k] * dpi1_390[k];

        t_391[k] = pa_x[k] * dpi0_391[k]
                   - f_11 * pc_x[k] * dpi1_391[k];

        t_392[k] = pa_y[k] * dpi0_224[k]
                   - f_11 * pc_y[k] * dpi1_224[k];

        t_393[k] = f_1 * dph_168[k]
                   + f_4 * pc_y[k] * fph_294[k];
    }

#pragma omp simd aligned(t_394, t_395, t_396, pa_y, pc_x, pc_y, dpi0_226, dph_170, dph_297, \
                         dpi1_226, fpg0_213, fpg1_213, fph_296, \
                         fph_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_394[k] = pa_y[k] * dpi0_226[k]
                   - f_11 * pc_y[k] * dpi1_226[k];

        t_395[k] = f_1 * dph_297[k]
                   + f_9 * fpg0_213[k]
                   - f_10 * fpg1_213[k]
                   + f_4 * pc_x[k] * fph_297[k];

        t_396[k] = f_1 * dph_170[k]
                   + f_4 * pc_y[k] * fph_296[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pa_y, pc_x, pc_y, pc_z, dpi0_229, dph_108, \
                         dph_300, dpi1_229, fsh_87, fpg0_216, fpg1_216, fph_297, \
                         fph_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = pa_y[k] * dpi0_229[k]
                   - f_11 * pc_y[k] * dpi1_229[k];

        t_398[k] = f_1 * dph_300[k]
                   + f_7 * fpg0_216[k]
                   - f_8 * fpg1_216[k]
                   + f_4 * pc_x[k] * fph_300[k];

        t_399[k] = f_1 * dph_108[k]
                   + f_1 * fsh_87[k]
                   + f_4 * pc_z[k] * fph_297[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, pa_y, pc_x, pc_y, dpi0_233, dph_173, dph_304, \
                         dpi1_233, fpg0_220, fpg1_220, fph_299, \
                         fph_304 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_1 * dph_173[k]
                   + f_4 * pc_y[k] * fph_299[k];

        t_401[k] = pa_y[k] * dpi0_233[k]
                   - f_11 * pc_y[k] * dpi1_233[k];

        t_402[k] = f_1 * dph_304[k]
                   + f_5 * fpg0_220[k]
                   - f_6 * fpg1_220[k]
                   + f_4 * pc_x[k] * fph_304[k];
    }

#pragma omp simd aligned(t_403, t_404, t_405, pa_x, pc_x, pc_y, pc_z, dpi0_404, dph_111, \
                         dph_177, dph_306, dpi1_404, fsh_90, fph_300, \
                         fph_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_403[k] = f_1 * dph_111[k]
                   + f_1 * fsh_90[k]
                   + f_4 * pc_z[k] * fph_300[k];

        t_404[k] = pa_x[k] * dpi0_404[k]
                   + f_12 * dph_306[k]
                   - f_11 * pc_x[k] * dpi1_404[k];

        t_405[k] = f_1 * dph_177[k]
                   + f_4 * pc_y[k] * fph_303[k];
    }

#pragma omp simd aligned(t_406, t_407, t_408, t_409, pa_y, pc_x, pc_y, dpi0_238, dph_309, \
                         dph_310, dph_311, dpi1_238, fph_309, fph_310, \
                         fph_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_406[k] = pa_y[k] * dpi0_238[k]
                   - f_11 * pc_y[k] * dpi1_238[k];

        t_407[k] = f_1 * dph_309[k]
                   + f_4 * pc_x[k] * fph_309[k];

        t_408[k] = f_1 * dph_310[k]
                   + f_4 * pc_x[k] * fph_310[k];

        t_409[k] = f_1 * dph_311[k]
                   + f_4 * pc_x[k] * fph_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pa_x, pc_x, dpi0_413, dph_312, dph_313, \
                         dph_314, dpi1_413, fph_312, fph_313, fph_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_1 * dph_312[k]
                   + f_4 * pc_x[k] * fph_312[k];

        t_411[k] = f_1 * dph_313[k]
                   + f_4 * pc_x[k] * fph_313[k];

        t_412[k] = f_1 * dph_314[k]
                   + f_4 * pc_x[k] * fph_314[k];

        t_413[k] = pa_x[k] * dpi0_413[k]
                   - f_11 * pc_x[k] * dpi1_413[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, t_417, pa_x, pc_x, dpi0_414, dpi0_415, dpi0_416, \
                         dpi0_417, dpi1_414, dpi1_415, dpi1_416, \
                         dpi1_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = pa_x[k] * dpi0_414[k]
                   - f_11 * pc_x[k] * dpi1_414[k];

        t_415[k] = pa_x[k] * dpi0_415[k]
                   - f_11 * pc_x[k] * dpi1_415[k];

        t_416[k] = pa_x[k] * dpi0_416[k]
                   - f_11 * pc_x[k] * dpi1_416[k];

        t_417[k] = pa_x[k] * dpi0_417[k]
                   - f_11 * pc_x[k] * dpi1_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, t_421, pa_x, pc_x, pc_y, dpi0_419, dph_188, \
                         dph_315, dpi1_419, fsh_105, fpg0_225, fpg1_225, fph_314, \
                         fph_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_1 * dph_188[k]
                   + f_4 * pc_y[k] * fph_314[k];

        t_419[k] = pa_x[k] * dpi0_419[k]
                   - f_11 * pc_x[k] * dpi1_419[k];

        t_420[k] = f_1 * dph_315[k]
                   + f_1 * fsh_105[k]
                   + f_2 * fpg0_225[k]
                   - f_3 * fpg1_225[k]
                   + f_4 * pc_x[k] * fph_315[k];

        t_421[k] = f_4 * pc_y[k] * fph_315[k];
    }

#pragma omp simd aligned(t_422, t_423, t_424, pc_y, pc_z, dph_126, fpg0_225, fpg1_225, \
                         fph_315, fph_316, fph_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_422[k] = f_12 * dph_126[k]
                   + f_4 * pc_z[k] * fph_315[k];

        t_423[k] = f_5 * fpg0_225[k]
                   - f_6 * fpg1_225[k]
                   + f_4 * pc_y[k] * fph_316[k];

        t_424[k] = f_4 * pc_y[k] * fph_317[k];
    }

#pragma omp simd aligned(t_425, t_426, t_427, pc_x, pc_y, dph_320, fsh_110, fpg0_226, \
                         fpg0_227, fpg0_230, fpg1_226, fpg1_227, fpg1_230, fph_318, fph_319, \
                         fph_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_425[k] = f_1 * dph_320[k]
                   + f_1 * fsh_110[k]
                   + f_9 * fpg0_230[k]
                   - f_10 * fpg1_230[k]
                   + f_4 * pc_x[k] * fph_320[k];

        t_426[k] = f_7 * fpg0_226[k]
                   - f_8 * fpg1_226[k]
                   + f_4 * pc_y[k] * fph_318[k];

        t_427[k] = f_5 * fpg0_227[k]
                   - f_6 * fpg1_227[k]
                   + f_4 * pc_y[k] * fph_319[k];
    }

#pragma omp simd aligned(t_428, t_429, t_430, pc_x, pc_y, dph_324, fsh_114, fpg0_228, \
                         fpg0_234, fpg1_228, fpg1_234, fph_320, fph_321, \
                         fph_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_428[k] = f_4 * pc_y[k] * fph_320[k];

        t_429[k] = f_1 * dph_324[k]
                   + f_1 * fsh_114[k]
                   + f_7 * fpg0_234[k]
                   - f_8 * fpg1_234[k]
                   + f_4 * pc_x[k] * fph_324[k];

        t_430[k] = f_9 * fpg0_228[k]
                   - f_10 * fpg1_228[k]
                   + f_4 * pc_y[k] * fph_321[k];
    }

#pragma omp simd aligned(t_431, t_432, t_433, pc_y, fpg0_229, fpg0_230, fpg1_229, fpg1_230, \
                         fph_322, fph_323, fph_324 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_431[k] = f_7 * fpg0_229[k]
                   - f_8 * fpg1_229[k]
                   + f_4 * pc_y[k] * fph_322[k];

        t_432[k] = f_5 * fpg0_230[k]
                   - f_6 * fpg1_230[k]
                   + f_4 * pc_y[k] * fph_323[k];

        t_433[k] = f_4 * pc_y[k] * fph_324[k];
    }

#pragma omp simd aligned(t_434, t_435, t_436, pc_x, dph_329, dph_330, dph_331, fsh_119, \
                         fsh_120, fsh_121, fpg0_239, fpg1_239, fph_329, fph_330, \
                         fph_331 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_434[k] = f_1 * dph_329[k]
                   + f_1 * fsh_119[k]
                   + f_5 * fpg0_239[k]
                   - f_6 * fpg1_239[k]
                   + f_4 * pc_x[k] * fph_329[k];

        t_435[k] = f_1 * dph_330[k]
                   + f_1 * fsh_120[k]
                   + f_4 * pc_x[k] * fph_330[k];

        t_436[k] = f_1 * dph_331[k]
                   + f_1 * fsh_121[k]
                   + f_4 * pc_x[k] * fph_331[k];
    }

#pragma omp simd aligned(t_437, t_438, t_439, t_440, pc_x, pc_y, dph_332, dph_333, dph_335, \
                         fsh_122, fsh_123, fsh_125, fph_329, fph_332, fph_333, \
                         fph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_437[k] = f_1 * dph_332[k]
                   + f_1 * fsh_122[k]
                   + f_4 * pc_x[k] * fph_332[k];

        t_438[k] = f_1 * dph_333[k]
                   + f_1 * fsh_123[k]
                   + f_4 * pc_x[k] * fph_333[k];

        t_439[k] = f_4 * pc_y[k] * fph_329[k];

        t_440[k] = f_1 * dph_335[k]
                   + f_1 * fsh_125[k]
                   + f_4 * pc_x[k] * fph_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, pc_y, fpg0_235, fpg0_236, fpg0_237, fpg1_235, \
                         fpg1_236, fpg1_237, fph_330, fph_331, \
                         fph_332 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = f_2 * fpg0_235[k]
                   - f_3 * fpg1_235[k]
                   + f_4 * pc_y[k] * fph_330[k];

        t_442[k] = f_17 * fpg0_236[k]
                   - f_18 * fpg1_236[k]
                   + f_4 * pc_y[k] * fph_331[k];

        t_443[k] = f_9 * fpg0_237[k]
                   - f_10 * fpg1_237[k]
                   + f_4 * pc_y[k] * fph_332[k];
    }

#pragma omp simd aligned(t_444, t_445, t_446, t_447, pc_y, pc_z, dph_146, fpg0_238, fpg0_239, \
                         fpg1_238, fpg1_239, fph_333, fph_334, \
                         fph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_444[k] = f_7 * fpg0_238[k]
                   - f_8 * fpg1_238[k]
                   + f_4 * pc_y[k] * fph_333[k];

        t_445[k] = f_5 * fpg0_239[k]
                   - f_6 * fpg1_239[k]
                   + f_4 * pc_y[k] * fph_334[k];

        t_446[k] = f_4 * pc_y[k] * fph_335[k];

        t_447[k] = f_12 * dph_146[k]
                   + f_2 * fpg0_239[k]
                   - f_3 * fpg1_239[k]
                   + f_4 * pc_z[k] * fph_335[k];
    }

#pragma omp simd aligned(t_448, t_449, t_450, pb_y, pc_y, pc_z, dph_147, fsi0_140, fsh_105, \
                         fsi1_140, fph_336 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_448[k] = pb_y[k] * fsi0_140[k]
                   - f_11 * pc_y[k] * fsi1_140[k];

        t_449[k] = f_1 * fsh_105[k]
                   + f_4 * pc_y[k] * fph_336[k];

        t_450[k] = f_12 * dph_147[k]
                   + f_4 * pc_z[k] * fph_336[k];
    }

#pragma omp simd aligned(t_451, t_452, t_453, pa_x, pb_y, pc_x, pc_y, dpi0_451, dph_339, \
                         dpi1_451, fsi0_145, fsh_107, fsi1_145, \
                         fph_338 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_451[k] = pa_x[k] * dpi0_451[k]
                   + f_13 * dph_339[k]
                   - f_11 * pc_x[k] * dpi1_451[k];

        t_452[k] = f_1 * fsh_107[k]
                   + f_4 * pc_y[k] * fph_338[k];

        t_453[k] = pb_y[k] * fsi0_145[k]
                   - f_11 * pc_y[k] * fsi1_145[k];
    }

#pragma omp simd aligned(t_454, t_455, t_456, pa_x, pc_x, pc_y, dpi0_454, dpi0_455, dph_342, \
                         dph_343, dpi1_454, dpi1_455, fsh_110, \
                         fph_341 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_454[k] = pa_x[k] * dpi0_454[k]
                   + f_0 * dph_342[k]
                   - f_11 * pc_x[k] * dpi1_454[k];

        t_455[k] = pa_x[k] * dpi0_455[k]
                   + f_0 * dph_343[k]
                   - f_11 * pc_x[k] * dpi1_455[k];

        t_456[k] = f_1 * fsh_110[k]
                   + f_4 * pc_y[k] * fph_341[k];
    }

#pragma omp simd aligned(t_457, t_458, t_459, pa_x, pb_y, pc_x, pc_y, dpi0_458, dpi0_459, \
                         dph_346, dph_347, dpi1_458, dpi1_459, fsi0_149, \
                         fsi1_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_457[k] = pb_y[k] * fsi0_149[k]
                   - f_11 * pc_y[k] * fsi1_149[k];

        t_458[k] = pa_x[k] * dpi0_458[k]
                   + f_12 * dph_346[k]
                   - f_11 * pc_x[k] * dpi1_458[k];

        t_459[k] = pa_x[k] * dpi0_459[k]
                   + f_12 * dph_347[k]
                   - f_11 * pc_x[k] * dpi1_459[k];
    }

#pragma omp simd aligned(t_460, t_461, t_462, pa_x, pb_y, pc_x, pc_y, dpi0_460, dph_348, \
                         dpi1_460, fsi0_154, fsh_114, fsi1_154, \
                         fph_345 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_460[k] = pa_x[k] * dpi0_460[k]
                   + f_12 * dph_348[k]
                   - f_11 * pc_x[k] * dpi1_460[k];

        t_461[k] = f_1 * fsh_114[k]
                   + f_4 * pc_y[k] * fph_345[k];

        t_462[k] = pb_y[k] * fsi0_154[k]
                   - f_11 * pc_y[k] * fsi1_154[k];
    }

#pragma omp simd aligned(t_463, t_464, t_465, t_466, pc_x, dph_351, dph_352, dph_353, dph_354, \
                         fph_351, fph_352, fph_353, fph_354 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_463[k] = f_1 * dph_351[k]
                   + f_4 * pc_x[k] * fph_351[k];

        t_464[k] = f_1 * dph_352[k]
                   + f_4 * pc_x[k] * fph_352[k];

        t_465[k] = f_1 * dph_353[k]
                   + f_4 * pc_x[k] * fph_353[k];

        t_466[k] = f_1 * dph_354[k]
                   + f_4 * pc_x[k] * fph_354[k];
    }

#pragma omp simd aligned(t_467, t_468, t_469, t_470, pa_x, pc_x, pc_y, dpi0_469, dpi0_470, \
                         dph_356, dpi1_469, dpi1_470, fsh_119, fph_350, \
                         fph_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_467[k] = f_1 * fsh_119[k]
                   + f_4 * pc_y[k] * fph_350[k];

        t_468[k] = f_1 * dph_356[k]
                   + f_4 * pc_x[k] * fph_356[k];

        t_469[k] = pa_x[k] * dpi0_469[k]
                   - f_11 * pc_x[k] * dpi1_469[k];

        t_470[k] = pa_x[k] * dpi0_470[k]
                   - f_11 * pc_x[k] * dpi1_470[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, t_474, pa_x, pc_x, pc_y, dpi0_471, dpi0_472, \
                         dpi0_473, dpi1_471, dpi1_472, dpi1_473, fsh_125, \
                         fph_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pa_x[k] * dpi0_471[k]
                   - f_11 * pc_x[k] * dpi1_471[k];

        t_472[k] = pa_x[k] * dpi0_472[k]
                   - f_11 * pc_x[k] * dpi1_472[k];

        t_473[k] = pa_x[k] * dpi0_473[k]
                   - f_11 * pc_x[k] * dpi1_473[k];

        t_474[k] = f_1 * fsh_125[k]
                   + f_4 * pc_y[k] * fph_356[k];
    }

#pragma omp simd aligned(t_475, t_476, t_477, t_478, pa_x, pc_x, pc_y, pc_z, dpi0_475, \
                         dpi0_476, dph_168, dph_357, dpi1_475, dpi1_476, fsh_105, \
                         fph_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_475[k] = pa_x[k] * dpi0_475[k]
                   - f_11 * pc_x[k] * dpi1_475[k];

        t_476[k] = pa_x[k] * dpi0_476[k]
                   + f_14 * dph_357[k]
                   - f_11 * pc_x[k] * dpi1_476[k];

        t_477[k] = f_4 * pc_y[k] * fph_357[k];

        t_478[k] = f_12 * dph_168[k]
                   + f_1 * fsh_105[k]
                   + f_4 * pc_z[k] * fph_357[k];
    }

#pragma omp simd aligned(t_479, t_480, t_481, pa_x, pc_x, pc_y, dpi0_481, dph_362, dpi1_481, \
                         fpg0_255, fpg1_255, fph_358, fph_359 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_479[k] = f_5 * fpg0_255[k]
                   - f_6 * fpg1_255[k]
                   + f_4 * pc_y[k] * fph_358[k];

        t_480[k] = f_4 * pc_y[k] * fph_359[k];

        t_481[k] = pa_x[k] * dpi0_481[k]
                   + f_13 * dph_362[k]
                   - f_11 * pc_x[k] * dpi1_481[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_252 = buffer.data(dpi0 + 252);
    const auto *dpi0_253 = buffer.data(dpi0 + 253);
    const auto *dpi0_255 = buffer.data(dpi0 + 255);
    const auto *dpi0_256 = buffer.data(dpi0 + 256);
    const auto *dpi0_258 = buffer.data(dpi0 + 258);
    const auto *dpi0_259 = buffer.data(dpi0 + 259);
    const auto *dpi0_260 = buffer.data(dpi0 + 260);
    const auto *dpi0_262 = buffer.data(dpi0 + 262);
    const auto *dpi0_263 = buffer.data(dpi0 + 263);
    const auto *dpi0_264 = buffer.data(dpi0 + 264);
    const auto *dpi0_265 = buffer.data(dpi0 + 265);
    const auto *dpi0_485 = buffer.data(dpi0 + 485);
    const auto *dpi0_490 = buffer.data(dpi0 + 490);
    const auto *dpi0_497 = buffer.data(dpi0 + 497);
    const auto *dpi0_498 = buffer.data(dpi0 + 498);
    const auto *dpi0_499 = buffer.data(dpi0 + 499);
    const auto *dpi0_500 = buffer.data(dpi0 + 500);
    const auto *dpi0_501 = buffer.data(dpi0 + 501);
    const auto *dpi0_503 = buffer.data(dpi0 + 503);

    const auto *dph_190 = buffer.data(dph + 190);
    const auto *dph_192 = buffer.data(dph + 192);
    const auto *dph_193 = buffer.data(dph + 193);
    const auto *dph_195 = buffer.data(dph + 195);
    const auto *dph_196 = buffer.data(dph + 196);
    const auto *dph_197 = buffer.data(dph + 197);
    const auto *dph_209 = buffer.data(dph + 209);
    const auto *dph_225 = buffer.data(dph + 225);
    const auto *dph_230 = buffer.data(dph + 230);
    const auto *dph_251 = buffer.data(dph + 251);
    const auto *dph_366 = buffer.data(dph + 366);
    const auto *dph_371 = buffer.data(dph + 371);
    const auto *dph_372 = buffer.data(dph + 372);
    const auto *dph_373 = buffer.data(dph + 373);
    const auto *dph_374 = buffer.data(dph + 374);
    const auto *dph_375 = buffer.data(dph + 375);
    const auto *dph_377 = buffer.data(dph + 377);

    const auto *dpi1_252 = buffer.data(dpi1 + 252);
    const auto *dpi1_253 = buffer.data(dpi1 + 253);
    const auto *dpi1_255 = buffer.data(dpi1 + 255);
    const auto *dpi1_256 = buffer.data(dpi1 + 256);
    const auto *dpi1_258 = buffer.data(dpi1 + 258);
    const auto *dpi1_259 = buffer.data(dpi1 + 259);
    const auto *dpi1_260 = buffer.data(dpi1 + 260);
    const auto *dpi1_262 = buffer.data(dpi1 + 262);
    const auto *dpi1_263 = buffer.data(dpi1 + 263);
    const auto *dpi1_264 = buffer.data(dpi1 + 264);
    const auto *dpi1_265 = buffer.data(dpi1 + 265);
    const auto *dpi1_485 = buffer.data(dpi1 + 485);
    const auto *dpi1_490 = buffer.data(dpi1 + 490);
    const auto *dpi1_497 = buffer.data(dpi1 + 497);
    const auto *dpi1_498 = buffer.data(dpi1 + 498);
    const auto *dpi1_499 = buffer.data(dpi1 + 499);
    const auto *dpi1_500 = buffer.data(dpi1 + 500);
    const auto *dpi1_501 = buffer.data(dpi1 + 501);
    const auto *dpi1_503 = buffer.data(dpi1 + 503);

    const auto *fsi0_168 = buffer.data(fsi0 + 168);
    const auto *fsi0_169 = buffer.data(fsi0 + 169);
    const auto *fsi0_171 = buffer.data(fsi0 + 171);
    const auto *fsi0_173 = buffer.data(fsi0 + 173);
    const auto *fsi0_174 = buffer.data(fsi0 + 174);
    const auto *fsi0_176 = buffer.data(fsi0 + 176);
    const auto *fsi0_177 = buffer.data(fsi0 + 177);
    const auto *fsi0_178 = buffer.data(fsi0 + 178);
    const auto *fsi0_180 = buffer.data(fsi0 + 180);
    const auto *fsi0_181 = buffer.data(fsi0 + 181);
    const auto *fsi0_182 = buffer.data(fsi0 + 182);
    const auto *fsi0_189 = buffer.data(fsi0 + 189);
    const auto *fsi0_191 = buffer.data(fsi0 + 191);
    const auto *fsi0_192 = buffer.data(fsi0 + 192);
    const auto *fsi0_193 = buffer.data(fsi0 + 193);
    const auto *fsi0_195 = buffer.data(fsi0 + 195);
    const auto *fsi0_198 = buffer.data(fsi0 + 198);
    const auto *fsi0_201 = buffer.data(fsi0 + 201);
    const auto *fsi0_205 = buffer.data(fsi0 + 205);
    const auto *fsi0_210 = buffer.data(fsi0 + 210);

    const auto *fsh_126 = buffer.data(fsh + 126);
    const auto *fsh_127 = buffer.data(fsh + 127);
    const auto *fsh_129 = buffer.data(fsh + 129);
    const auto *fsh_131 = buffer.data(fsh + 131);
    const auto *fsh_132 = buffer.data(fsh + 132);
    const auto *fsh_134 = buffer.data(fsh + 134);
    const auto *fsh_135 = buffer.data(fsh + 135);
    const auto *fsh_136 = buffer.data(fsh + 136);
    const auto *fsh_138 = buffer.data(fsh + 138);
    const auto *fsh_139 = buffer.data(fsh + 139);
    const auto *fsh_140 = buffer.data(fsh + 140);
    const auto *fsh_141 = buffer.data(fsh + 141);
    const auto *fsh_142 = buffer.data(fsh + 142);
    const auto *fsh_143 = buffer.data(fsh + 143);
    const auto *fsh_144 = buffer.data(fsh + 144);
    const auto *fsh_145 = buffer.data(fsh + 145);
    const auto *fsh_146 = buffer.data(fsh + 146);
    const auto *fsh_149 = buffer.data(fsh + 149);
    const auto *fsh_152 = buffer.data(fsh + 152);
    const auto *fsh_156 = buffer.data(fsh + 156);
    const auto *fsh_161 = buffer.data(fsh + 161);

    const auto *fsi1_168 = buffer.data(fsi1 + 168);
    const auto *fsi1_169 = buffer.data(fsi1 + 169);
    const auto *fsi1_171 = buffer.data(fsi1 + 171);
    const auto *fsi1_173 = buffer.data(fsi1 + 173);
    const auto *fsi1_174 = buffer.data(fsi1 + 174);
    const auto *fsi1_176 = buffer.data(fsi1 + 176);
    const auto *fsi1_177 = buffer.data(fsi1 + 177);
    const auto *fsi1_178 = buffer.data(fsi1 + 178);
    const auto *fsi1_180 = buffer.data(fsi1 + 180);
    const auto *fsi1_181 = buffer.data(fsi1 + 181);
    const auto *fsi1_182 = buffer.data(fsi1 + 182);
    const auto *fsi1_189 = buffer.data(fsi1 + 189);
    const auto *fsi1_191 = buffer.data(fsi1 + 191);
    const auto *fsi1_192 = buffer.data(fsi1 + 192);
    const auto *fsi1_193 = buffer.data(fsi1 + 193);
    const auto *fsi1_195 = buffer.data(fsi1 + 195);
    const auto *fsi1_198 = buffer.data(fsi1 + 198);
    const auto *fsi1_201 = buffer.data(fsi1 + 201);
    const auto *fsi1_205 = buffer.data(fsi1 + 205);
    const auto *fsi1_210 = buffer.data(fsi1 + 210);

    const auto *fpg0_256 = buffer.data(fpg0 + 256);
    const auto *fpg0_257 = buffer.data(fpg0 + 257);
    const auto *fpg0_258 = buffer.data(fpg0 + 258);
    const auto *fpg0_259 = buffer.data(fpg0 + 259);
    const auto *fpg0_260 = buffer.data(fpg0 + 260);
    const auto *fpg0_285 = buffer.data(fpg0 + 285);
    const auto *fpg0_286 = buffer.data(fpg0 + 286);
    const auto *fpg0_288 = buffer.data(fpg0 + 288);
    const auto *fpg0_290 = buffer.data(fpg0 + 290);
    const auto *fpg0_291 = buffer.data(fpg0 + 291);
    const auto *fpg0_293 = buffer.data(fpg0 + 293);
    const auto *fpg0_294 = buffer.data(fpg0 + 294);
    const auto *fpg0_295 = buffer.data(fpg0 + 295);
    const auto *fpg0_296 = buffer.data(fpg0 + 296);
    const auto *fpg0_297 = buffer.data(fpg0 + 297);
    const auto *fpg0_298 = buffer.data(fpg0 + 298);
    const auto *fpg0_299 = buffer.data(fpg0 + 299);
    const auto *fpg0_305 = buffer.data(fpg0 + 305);
    const auto *fpg0_308 = buffer.data(fpg0 + 308);
    const auto *fpg0_309 = buffer.data(fpg0 + 309);
    const auto *fpg0_312 = buffer.data(fpg0 + 312);
    const auto *fpg0_313 = buffer.data(fpg0 + 313);
    const auto *fpg0_314 = buffer.data(fpg0 + 314);

    const auto *fpg1_256 = buffer.data(fpg1 + 256);
    const auto *fpg1_257 = buffer.data(fpg1 + 257);
    const auto *fpg1_258 = buffer.data(fpg1 + 258);
    const auto *fpg1_259 = buffer.data(fpg1 + 259);
    const auto *fpg1_260 = buffer.data(fpg1 + 260);
    const auto *fpg1_285 = buffer.data(fpg1 + 285);
    const auto *fpg1_286 = buffer.data(fpg1 + 286);
    const auto *fpg1_288 = buffer.data(fpg1 + 288);
    const auto *fpg1_290 = buffer.data(fpg1 + 290);
    const auto *fpg1_291 = buffer.data(fpg1 + 291);
    const auto *fpg1_293 = buffer.data(fpg1 + 293);
    const auto *fpg1_294 = buffer.data(fpg1 + 294);
    const auto *fpg1_295 = buffer.data(fpg1 + 295);
    const auto *fpg1_296 = buffer.data(fpg1 + 296);
    const auto *fpg1_297 = buffer.data(fpg1 + 297);
    const auto *fpg1_298 = buffer.data(fpg1 + 298);
    const auto *fpg1_299 = buffer.data(fpg1 + 299);
    const auto *fpg1_305 = buffer.data(fpg1 + 305);
    const auto *fpg1_308 = buffer.data(fpg1 + 308);
    const auto *fpg1_309 = buffer.data(fpg1 + 309);
    const auto *fpg1_312 = buffer.data(fpg1 + 312);
    const auto *fpg1_313 = buffer.data(fpg1 + 313);
    const auto *fpg1_314 = buffer.data(fpg1 + 314);

    const auto *fph_360 = buffer.data(fph + 360);
    const auto *fph_361 = buffer.data(fph + 361);
    const auto *fph_362 = buffer.data(fph + 362);
    const auto *fph_363 = buffer.data(fph + 363);
    const auto *fph_364 = buffer.data(fph + 364);
    const auto *fph_365 = buffer.data(fph + 365);
    const auto *fph_366 = buffer.data(fph + 366);
    const auto *fph_371 = buffer.data(fph + 371);
    const auto *fph_372 = buffer.data(fph + 372);
    const auto *fph_373 = buffer.data(fph + 373);
    const auto *fph_374 = buffer.data(fph + 374);
    const auto *fph_375 = buffer.data(fph + 375);
    const auto *fph_377 = buffer.data(fph + 377);
    const auto *fph_378 = buffer.data(fph + 378);
    const auto *fph_379 = buffer.data(fph + 379);
    const auto *fph_381 = buffer.data(fph + 381);
    const auto *fph_384 = buffer.data(fph + 384);
    const auto *fph_393 = buffer.data(fph + 393);
    const auto *fph_394 = buffer.data(fph + 394);
    const auto *fph_395 = buffer.data(fph + 395);
    const auto *fph_396 = buffer.data(fph + 396);
    const auto *fph_397 = buffer.data(fph + 397);
    const auto *fph_398 = buffer.data(fph + 398);
    const auto *fph_399 = buffer.data(fph + 399);
    const auto *fph_400 = buffer.data(fph + 400);
    const auto *fph_402 = buffer.data(fph + 402);
    const auto *fph_404 = buffer.data(fph + 404);
    const auto *fph_405 = buffer.data(fph + 405);
    const auto *fph_407 = buffer.data(fph + 407);
    const auto *fph_408 = buffer.data(fph + 408);
    const auto *fph_409 = buffer.data(fph + 409);
    const auto *fph_411 = buffer.data(fph + 411);
    const auto *fph_412 = buffer.data(fph + 412);
    const auto *fph_413 = buffer.data(fph + 413);
    const auto *fph_414 = buffer.data(fph + 414);
    const auto *fph_415 = buffer.data(fph + 415);
    const auto *fph_416 = buffer.data(fph + 416);
    const auto *fph_417 = buffer.data(fph + 417);
    const auto *fph_418 = buffer.data(fph + 418);
    const auto *fph_419 = buffer.data(fph + 419);
    const auto *fph_420 = buffer.data(fph + 420);
    const auto *fph_421 = buffer.data(fph + 421);
    const auto *fph_423 = buffer.data(fph + 423);
    const auto *fph_425 = buffer.data(fph + 425);
    const auto *fph_426 = buffer.data(fph + 426);
    const auto *fph_428 = buffer.data(fph + 428);
    const auto *fph_429 = buffer.data(fph + 429);
    const auto *fph_432 = buffer.data(fph + 432);
    const auto *fph_433 = buffer.data(fph + 433);
    const auto *fph_434 = buffer.data(fph + 434);
    const auto *fph_435 = buffer.data(fph + 435);
    const auto *fph_436 = buffer.data(fph + 436);
    const auto *fph_437 = buffer.data(fph + 437);
    const auto *fph_438 = buffer.data(fph + 438);
    const auto *fph_439 = buffer.data(fph + 439);
    const auto *fph_440 = buffer.data(fph + 440);

#pragma omp simd aligned(t_482, t_483, t_484, pc_y, fpg0_256, fpg0_257, fpg1_256, fpg1_257, \
                         fph_360, fph_361, fph_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_7 * fpg0_256[k]
                   - f_8 * fpg1_256[k]
                   + f_4 * pc_y[k] * fph_360[k];

        t_483[k] = f_5 * fpg0_257[k]
                   - f_6 * fpg1_257[k]
                   + f_4 * pc_y[k] * fph_361[k];

        t_484[k] = f_4 * pc_y[k] * fph_362[k];
    }

#pragma omp simd aligned(t_485, t_486, t_487, pa_x, pc_x, pc_y, dpi0_485, dph_366, dpi1_485, \
                         fpg0_258, fpg0_259, fpg1_258, fpg1_259, fph_363, \
                         fph_364 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_485[k] = pa_x[k] * dpi0_485[k]
                   + f_0 * dph_366[k]
                   - f_11 * pc_x[k] * dpi1_485[k];

        t_486[k] = f_9 * fpg0_258[k]
                   - f_10 * fpg1_258[k]
                   + f_4 * pc_y[k] * fph_363[k];

        t_487[k] = f_7 * fpg0_259[k]
                   - f_8 * fpg1_259[k]
                   + f_4 * pc_y[k] * fph_364[k];
    }

#pragma omp simd aligned(t_488, t_489, t_490, t_491, pa_x, pc_x, pc_y, dpi0_490, dph_371, \
                         dph_372, dpi1_490, fpg0_260, fpg1_260, fph_365, fph_366, \
                         fph_372 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_488[k] = f_5 * fpg0_260[k]
                   - f_6 * fpg1_260[k]
                   + f_4 * pc_y[k] * fph_365[k];

        t_489[k] = f_4 * pc_y[k] * fph_366[k];

        t_490[k] = pa_x[k] * dpi0_490[k]
                   + f_12 * dph_371[k]
                   - f_11 * pc_x[k] * dpi1_490[k];

        t_491[k] = f_1 * dph_372[k]
                   + f_4 * pc_x[k] * fph_372[k];
    }

#pragma omp simd aligned(t_492, t_493, t_494, t_495, t_496, pc_x, pc_y, dph_373, dph_374, \
                         dph_375, dph_377, fph_371, fph_373, fph_374, fph_375, \
                         fph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_492[k] = f_1 * dph_373[k]
                   + f_4 * pc_x[k] * fph_373[k];

        t_493[k] = f_1 * dph_374[k]
                   + f_4 * pc_x[k] * fph_374[k];

        t_494[k] = f_1 * dph_375[k]
                   + f_4 * pc_x[k] * fph_375[k];

        t_495[k] = f_4 * pc_y[k] * fph_371[k];

        t_496[k] = f_1 * dph_377[k]
                   + f_4 * pc_x[k] * fph_377[k];
    }

#pragma omp simd aligned(t_497, t_498, t_499, t_500, pa_x, pc_x, dpi0_497, dpi0_498, dpi0_499, \
                         dpi0_500, dpi1_497, dpi1_498, dpi1_499, \
                         dpi1_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_497[k] = pa_x[k] * dpi0_497[k]
                   - f_11 * pc_x[k] * dpi1_497[k];

        t_498[k] = pa_x[k] * dpi0_498[k]
                   - f_11 * pc_x[k] * dpi1_498[k];

        t_499[k] = pa_x[k] * dpi0_499[k]
                   - f_11 * pc_x[k] * dpi1_499[k];

        t_500[k] = pa_x[k] * dpi0_500[k]
                   - f_11 * pc_x[k] * dpi1_500[k];
    }

#pragma omp simd aligned(t_501, t_502, t_503, t_504, pa_x, pb_x, pc_x, pc_y, dpi0_501, \
                         dpi0_503, dpi1_501, dpi1_503, fsi0_168, fsh_126, fsi1_168, \
                         fph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_501[k] = pa_x[k] * dpi0_501[k]
                   - f_11 * pc_x[k] * dpi1_501[k];

        t_502[k] = f_4 * pc_y[k] * fph_377[k];

        t_503[k] = pa_x[k] * dpi0_503[k]
                   - f_11 * pc_x[k] * dpi1_503[k];

        t_504[k] = pb_x[k] * fsi0_168[k]
                   + f_14 * fsh_126[k]
                   - f_11 * pc_x[k] * fsi1_168[k];
    }

#pragma omp simd aligned(t_505, t_506, t_507, t_508, pb_x, pc_x, pc_z, fsi0_169, fsi0_171, \
                         fsh_127, fsh_129, fsi1_169, fsi1_171, fph_378, \
                         fph_379 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_505[k] = pb_x[k] * fsi0_169[k]
                   + f_19 * fsh_127[k]
                   - f_11 * pc_x[k] * fsi1_169[k];

        t_506[k] = f_4 * pc_z[k] * fph_378[k];

        t_507[k] = pb_x[k] * fsi0_171[k]
                   + f_13 * fsh_129[k]
                   - f_11 * pc_x[k] * fsi1_171[k];

        t_508[k] = f_4 * pc_z[k] * fph_379[k];
    }

#pragma omp simd aligned(t_509, t_510, t_511, pb_x, pc_x, pc_z, fsi0_173, fsi0_174, fsh_131, \
                         fsh_132, fsi1_173, fsi1_174, fph_381 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_509[k] = pb_x[k] * fsi0_173[k]
                   + f_13 * fsh_131[k]
                   - f_11 * pc_x[k] * fsi1_173[k];

        t_510[k] = pb_x[k] * fsi0_174[k]
                   + f_0 * fsh_132[k]
                   - f_11 * pc_x[k] * fsi1_174[k];

        t_511[k] = f_4 * pc_z[k] * fph_381[k];
    }

#pragma omp simd aligned(t_512, t_513, t_514, pb_x, pc_x, fsi0_176, fsi0_177, fsi0_178, \
                         fsh_134, fsh_135, fsh_136, fsi1_176, fsi1_177, \
                         fsi1_178 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_512[k] = pb_x[k] * fsi0_176[k]
                   + f_0 * fsh_134[k]
                   - f_11 * pc_x[k] * fsi1_176[k];

        t_513[k] = pb_x[k] * fsi0_177[k]
                   + f_0 * fsh_135[k]
                   - f_11 * pc_x[k] * fsi1_177[k];

        t_514[k] = pb_x[k] * fsi0_178[k]
                   + f_12 * fsh_136[k]
                   - f_11 * pc_x[k] * fsi1_178[k];
    }

#pragma omp simd aligned(t_515, t_516, t_517, pb_x, pc_x, pc_z, fsi0_180, fsi0_181, fsh_138, \
                         fsh_139, fsi1_180, fsi1_181, fph_384 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_515[k] = f_4 * pc_z[k] * fph_384[k];

        t_516[k] = pb_x[k] * fsi0_180[k]
                   + f_12 * fsh_138[k]
                   - f_11 * pc_x[k] * fsi1_180[k];

        t_517[k] = pb_x[k] * fsi0_181[k]
                   + f_12 * fsh_139[k]
                   - f_11 * pc_x[k] * fsi1_181[k];
    }

#pragma omp simd aligned(t_518, t_519, t_520, t_521, pb_x, pc_x, fsi0_182, fsh_140, fsh_141, \
                         fsh_142, fsh_143, fsi1_182, fph_393, fph_394, \
                         fph_395 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_518[k] = pb_x[k] * fsi0_182[k]
                   + f_12 * fsh_140[k]
                   - f_11 * pc_x[k] * fsi1_182[k];

        t_519[k] = f_1 * fsh_141[k]
                   + f_4 * pc_x[k] * fph_393[k];

        t_520[k] = f_1 * fsh_142[k]
                   + f_4 * pc_x[k] * fph_394[k];

        t_521[k] = f_1 * fsh_143[k]
                   + f_4 * pc_x[k] * fph_395[k];
    }

#pragma omp simd aligned(t_522, t_523, t_524, t_525, pb_x, pc_x, fsi0_189, fsh_144, fsh_145, \
                         fsh_146, fsi1_189, fph_396, fph_397, fph_398 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_522[k] = f_1 * fsh_144[k]
                   + f_4 * pc_x[k] * fph_396[k];

        t_523[k] = f_1 * fsh_145[k]
                   + f_4 * pc_x[k] * fph_397[k];

        t_524[k] = f_1 * fsh_146[k]
                   + f_4 * pc_x[k] * fph_398[k];

        t_525[k] = pb_x[k] * fsi0_189[k]
                   - f_11 * pc_x[k] * fsi1_189[k];
    }

#pragma omp simd aligned(t_526, t_527, t_528, t_529, pb_x, pc_x, pc_z, fsi0_191, fsi0_192, \
                         fsi0_193, fsi1_191, fsi1_192, fsi1_193, \
                         fph_393 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_526[k] = f_4 * pc_z[k] * fph_393[k];

        t_527[k] = pb_x[k] * fsi0_191[k]
                   - f_11 * pc_x[k] * fsi1_191[k];

        t_528[k] = pb_x[k] * fsi0_192[k]
                   - f_11 * pc_x[k] * fsi1_192[k];

        t_529[k] = pb_x[k] * fsi0_193[k]
                   - f_11 * pc_x[k] * fsi1_193[k];
    }

#pragma omp simd aligned(t_530, t_531, t_532, pb_x, pc_x, pc_y, dph_209, fsi0_195, fsi1_195, \
                         fpg0_285, fpg1_285, fph_398, fph_399 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_530[k] = f_0 * dph_209[k]
                   + f_4 * pc_y[k] * fph_398[k];

        t_531[k] = pb_x[k] * fsi0_195[k]
                   - f_11 * pc_x[k] * fsi1_195[k];

        t_532[k] = f_2 * fpg0_285[k]
                   - f_3 * fpg1_285[k]
                   + f_4 * pc_x[k] * fph_399[k];
    }

#pragma omp simd aligned(t_533, t_534, t_535, t_536, pc_x, pc_z, fpg0_286, fpg0_288, fpg1_286, \
                         fpg1_288, fph_399, fph_400, fph_402 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_533[k] = f_17 * fpg0_286[k]
                   - f_18 * fpg1_286[k]
                   + f_4 * pc_x[k] * fph_400[k];

        t_534[k] = f_4 * pc_z[k] * fph_399[k];

        t_535[k] = f_9 * fpg0_288[k]
                   - f_10 * fpg1_288[k]
                   + f_4 * pc_x[k] * fph_402[k];

        t_536[k] = f_4 * pc_z[k] * fph_400[k];
    }

#pragma omp simd aligned(t_537, t_538, t_539, t_540, pc_x, pc_z, fpg0_290, fpg0_291, fpg0_293, \
                         fpg1_290, fpg1_291, fpg1_293, fph_402, fph_404, fph_405, \
                         fph_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_537[k] = f_9 * fpg0_290[k]
                   - f_10 * fpg1_290[k]
                   + f_4 * pc_x[k] * fph_404[k];

        t_538[k] = f_7 * fpg0_291[k]
                   - f_8 * fpg1_291[k]
                   + f_4 * pc_x[k] * fph_405[k];

        t_539[k] = f_4 * pc_z[k] * fph_402[k];

        t_540[k] = f_7 * fpg0_293[k]
                   - f_8 * fpg1_293[k]
                   + f_4 * pc_x[k] * fph_407[k];
    }

#pragma omp simd aligned(t_541, t_542, t_543, t_544, pc_x, pc_z, fpg0_294, fpg0_295, fpg0_297, \
                         fpg1_294, fpg1_295, fpg1_297, fph_405, fph_408, fph_409, \
                         fph_411 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_541[k] = f_7 * fpg0_294[k]
                   - f_8 * fpg1_294[k]
                   + f_4 * pc_x[k] * fph_408[k];

        t_542[k] = f_5 * fpg0_295[k]
                   - f_6 * fpg1_295[k]
                   + f_4 * pc_x[k] * fph_409[k];

        t_543[k] = f_4 * pc_z[k] * fph_405[k];

        t_544[k] = f_5 * fpg0_297[k]
                   - f_6 * fpg1_297[k]
                   + f_4 * pc_x[k] * fph_411[k];
    }

#pragma omp simd aligned(t_545, t_546, t_547, t_548, t_549, pc_x, fpg0_298, fpg0_299, \
                         fpg1_298, fpg1_299, fph_412, fph_413, fph_414, fph_415, \
                         fph_416 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_545[k] = f_5 * fpg0_298[k]
                   - f_6 * fpg1_298[k]
                   + f_4 * pc_x[k] * fph_412[k];

        t_546[k] = f_5 * fpg0_299[k]
                   - f_6 * fpg1_299[k]
                   + f_4 * pc_x[k] * fph_413[k];

        t_547[k] = f_4 * pc_x[k] * fph_414[k];

        t_548[k] = f_4 * pc_x[k] * fph_415[k];

        t_549[k] = f_4 * pc_x[k] * fph_416[k];
    }

#pragma omp simd aligned(t_550, t_551, t_552, t_553, t_554, pc_x, pc_y, pc_z, dph_225, \
                         fsh_141, fpg0_295, fpg1_295, fph_414, fph_417, fph_418, \
                         fph_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_550[k] = f_4 * pc_x[k] * fph_417[k];

        t_551[k] = f_4 * pc_x[k] * fph_418[k];

        t_552[k] = f_4 * pc_x[k] * fph_419[k];

        t_553[k] = f_0 * dph_225[k]
                   + f_1 * fsh_141[k]
                   + f_2 * fpg0_295[k]
                   - f_3 * fpg1_295[k]
                   + f_4 * pc_y[k] * fph_414[k];

        t_554[k] = f_4 * pc_z[k] * fph_414[k];
    }

#pragma omp simd aligned(t_555, t_556, t_557, pc_z, fpg0_295, fpg0_296, fpg0_297, fpg1_295, \
                         fpg1_296, fpg1_297, fph_415, fph_416, \
                         fph_417 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_555[k] = f_5 * fpg0_295[k]
                   - f_6 * fpg1_295[k]
                   + f_4 * pc_z[k] * fph_415[k];

        t_556[k] = f_7 * fpg0_296[k]
                   - f_8 * fpg1_296[k]
                   + f_4 * pc_z[k] * fph_416[k];

        t_557[k] = f_9 * fpg0_297[k]
                   - f_10 * fpg1_297[k]
                   + f_4 * pc_z[k] * fph_417[k];
    }

#pragma omp simd aligned(t_558, t_559, t_560, t_561, pb_z, pc_y, pc_z, dph_230, fsi0_168, \
                         fsi0_169, fsh_146, fsi1_168, fsi1_169, fpg0_299, fpg1_299, \
                         fph_419 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_558[k] = f_0 * dph_230[k]
                   + f_1 * fsh_146[k]
                   + f_4 * pc_y[k] * fph_419[k];

        t_559[k] = f_2 * fpg0_299[k]
                   - f_3 * fpg1_299[k]
                   + f_4 * pc_z[k] * fph_419[k];

        t_560[k] = pb_z[k] * fsi0_168[k]
                   - f_11 * pc_z[k] * fsi1_168[k];

        t_561[k] = pb_z[k] * fsi0_169[k]
                   - f_11 * pc_z[k] * fsi1_169[k];
    }

#pragma omp simd aligned(t_562, t_563, t_564, t_565, pb_z, pc_x, pc_z, fsi0_171, fsh_126, \
                         fsh_127, fsi1_171, fpg0_305, fpg1_305, fph_420, fph_421, \
                         fph_425 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_562[k] = f_1 * fsh_126[k]
                   + f_4 * pc_z[k] * fph_420[k];

        t_563[k] = pb_z[k] * fsi0_171[k]
                   - f_11 * pc_z[k] * fsi1_171[k];

        t_564[k] = f_1 * fsh_127[k]
                   + f_4 * pc_z[k] * fph_421[k];

        t_565[k] = f_9 * fpg0_305[k]
                   - f_10 * fpg1_305[k]
                   + f_4 * pc_x[k] * fph_425[k];
    }

#pragma omp simd aligned(t_566, t_567, t_568, pb_z, pc_x, pc_z, fsi0_174, fsh_129, fsi1_174, \
                         fpg0_308, fpg1_308, fph_423, fph_428 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_566[k] = pb_z[k] * fsi0_174[k]
                   - f_11 * pc_z[k] * fsi1_174[k];

        t_567[k] = f_1 * fsh_129[k]
                   + f_4 * pc_z[k] * fph_423[k];

        t_568[k] = f_7 * fpg0_308[k]
                   - f_8 * fpg1_308[k]
                   + f_4 * pc_x[k] * fph_428[k];
    }

#pragma omp simd aligned(t_569, t_570, t_571, pb_z, pc_x, pc_z, fsi0_178, fsh_132, fsi1_178, \
                         fpg0_309, fpg1_309, fph_426, fph_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_569[k] = f_7 * fpg0_309[k]
                   - f_8 * fpg1_309[k]
                   + f_4 * pc_x[k] * fph_429[k];

        t_570[k] = pb_z[k] * fsi0_178[k]
                   - f_11 * pc_z[k] * fsi1_178[k];

        t_571[k] = f_1 * fsh_132[k]
                   + f_4 * pc_z[k] * fph_426[k];
    }

#pragma omp simd aligned(t_572, t_573, t_574, t_575, pc_x, fpg0_312, fpg0_313, fpg0_314, \
                         fpg1_312, fpg1_313, fpg1_314, fph_432, fph_433, fph_434, \
                         fph_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_572[k] = f_5 * fpg0_312[k]
                   - f_6 * fpg1_312[k]
                   + f_4 * pc_x[k] * fph_432[k];

        t_573[k] = f_5 * fpg0_313[k]
                   - f_6 * fpg1_313[k]
                   + f_4 * pc_x[k] * fph_433[k];

        t_574[k] = f_5 * fpg0_314[k]
                   - f_6 * fpg1_314[k]
                   + f_4 * pc_x[k] * fph_434[k];

        t_575[k] = f_4 * pc_x[k] * fph_435[k];
    }

#pragma omp simd aligned(t_576, t_577, t_578, t_579, t_580, t_581, pb_z, pc_x, pc_z, fsi0_189, \
                         fsi1_189, fph_436, fph_437, fph_438, fph_439, \
                         fph_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_576[k] = f_4 * pc_x[k] * fph_436[k];

        t_577[k] = f_4 * pc_x[k] * fph_437[k];

        t_578[k] = f_4 * pc_x[k] * fph_438[k];

        t_579[k] = f_4 * pc_x[k] * fph_439[k];

        t_580[k] = f_4 * pc_x[k] * fph_440[k];

        t_581[k] = pb_z[k] * fsi0_189[k]
                   - f_11 * pc_z[k] * fsi1_189[k];
    }

#pragma omp simd aligned(t_582, t_583, t_584, pb_z, pc_z, fsi0_191, fsi0_192, fsh_141, \
                         fsh_142, fsh_143, fsi1_191, fsi1_192, \
                         fph_435 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_582[k] = f_1 * fsh_141[k]
                   + f_4 * pc_z[k] * fph_435[k];

        t_583[k] = pb_z[k] * fsi0_191[k]
                   + f_12 * fsh_142[k]
                   - f_11 * pc_z[k] * fsi1_191[k];

        t_584[k] = pb_z[k] * fsi0_192[k]
                   + f_0 * fsh_143[k]
                   - f_11 * pc_z[k] * fsi1_192[k];
    }

#pragma omp simd aligned(t_585, t_586, t_587, pb_z, pc_y, pc_z, dph_251, fsi0_193, fsi0_195, \
                         fsh_144, fsh_146, fsi1_193, fsi1_195, \
                         fph_440 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_585[k] = pb_z[k] * fsi0_193[k]
                   + f_13 * fsh_144[k]
                   - f_11 * pc_z[k] * fsi1_193[k];

        t_586[k] = f_0 * dph_251[k]
                   + f_4 * pc_y[k] * fph_440[k];

        t_587[k] = pb_z[k] * fsi0_195[k]
                   + f_14 * fsh_146[k]
                   - f_11 * pc_z[k] * fsi1_195[k];
    }

#pragma omp simd aligned(t_588, t_589, t_590, pa_z, pb_x, pc_x, pc_z, dpi0_252, dpi0_253, \
                         dpi1_252, dpi1_253, fsi0_198, fsh_149, \
                         fsi1_198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_588[k] = pa_z[k] * dpi0_252[k]
                   - f_11 * pc_z[k] * dpi1_252[k];

        t_589[k] = pa_z[k] * dpi0_253[k]
                   - f_11 * pc_z[k] * dpi1_253[k];

        t_590[k] = pb_x[k] * fsi0_198[k]
                   + f_19 * fsh_149[k]
                   - f_11 * pc_x[k] * fsi1_198[k];
    }

#pragma omp simd aligned(t_591, t_592, t_593, pa_z, pb_x, pc_x, pc_z, dpi0_255, dpi0_256, \
                         dph_190, dpi1_255, dpi1_256, fsi0_201, fsh_152, \
                         fsi1_201 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_591[k] = pa_z[k] * dpi0_255[k]
                   - f_11 * pc_z[k] * dpi1_255[k];

        t_592[k] = pa_z[k] * dpi0_256[k]
                   + f_1 * dph_190[k]
                   - f_11 * pc_z[k] * dpi1_256[k];

        t_593[k] = pb_x[k] * fsi0_201[k]
                   + f_13 * fsh_152[k]
                   - f_11 * pc_x[k] * fsi1_201[k];
    }

#pragma omp simd aligned(t_594, t_595, t_596, pa_z, pc_z, dpi0_258, dpi0_259, dpi0_260, \
                         dph_192, dph_193, dpi1_258, dpi1_259, \
                         dpi1_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_594[k] = pa_z[k] * dpi0_258[k]
                   - f_11 * pc_z[k] * dpi1_258[k];

        t_595[k] = pa_z[k] * dpi0_259[k]
                   + f_1 * dph_192[k]
                   - f_11 * pc_z[k] * dpi1_259[k];

        t_596[k] = pa_z[k] * dpi0_260[k]
                   + f_12 * dph_193[k]
                   - f_11 * pc_z[k] * dpi1_260[k];
    }

#pragma omp simd aligned(t_597, t_598, t_599, pa_z, pb_x, pc_x, pc_z, dpi0_262, dpi0_263, \
                         dph_195, dpi1_262, dpi1_263, fsi0_205, fsh_156, \
                         fsi1_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_597[k] = pb_x[k] * fsi0_205[k]
                   + f_0 * fsh_156[k]
                   - f_11 * pc_x[k] * fsi1_205[k];

        t_598[k] = pa_z[k] * dpi0_262[k]
                   - f_11 * pc_z[k] * dpi1_262[k];

        t_599[k] = pa_z[k] * dpi0_263[k]
                   + f_1 * dph_195[k]
                   - f_11 * pc_z[k] * dpi1_263[k];
    }

#pragma omp simd aligned(t_600, t_601, t_602, pa_z, pb_x, pc_x, pc_z, dpi0_264, dpi0_265, \
                         dph_196, dph_197, dpi1_264, dpi1_265, fsi0_210, fsh_161, \
                         fsi1_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_600[k] = pa_z[k] * dpi0_264[k]
                   + f_12 * dph_196[k]
                   - f_11 * pc_z[k] * dpi1_264[k];

        t_601[k] = pa_z[k] * dpi0_265[k]
                   + f_0 * dph_197[k]
                   - f_11 * pc_z[k] * dpi1_265[k];

        t_602[k] = pb_x[k] * fsi0_210[k]
                   + f_12 * fsh_161[k]
                   - f_11 * pc_x[k] * fsi1_210[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece5(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t ppi1,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_15 = 0.5 / p;
    const auto f_16 = 0.5 * gamma / (p * q);
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);

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
    auto *t_720 = buffer.data(target + 720);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *dpi0_273 = buffer.data(dpi0 + 273);
    const auto *dpi0_280 = buffer.data(dpi0 + 280);
    const auto *dpi0_281 = buffer.data(dpi0 + 281);
    const auto *dpi0_283 = buffer.data(dpi0 + 283);
    const auto *dpi0_286 = buffer.data(dpi0 + 286);
    const auto *dpi0_290 = buffer.data(dpi0 + 290);
    const auto *dpi0_301 = buffer.data(dpi0 + 301);
    const auto *dpi0_303 = buffer.data(dpi0 + 303);
    const auto *dpi0_304 = buffer.data(dpi0 + 304);
    const auto *dpi0_305 = buffer.data(dpi0 + 305);
    const auto *dpi0_419 = buffer.data(dpi0 + 419);
    const auto *dpi0_420 = buffer.data(dpi0 + 420);
    const auto *dpi0_421 = buffer.data(dpi0 + 421);
    const auto *dpi0_422 = buffer.data(dpi0 + 422);
    const auto *dpi0_423 = buffer.data(dpi0 + 423);
    const auto *dpi0_424 = buffer.data(dpi0 + 424);
    const auto *dpi0_425 = buffer.data(dpi0 + 425);
    const auto *dpi0_426 = buffer.data(dpi0 + 426);
    const auto *dpi0_427 = buffer.data(dpi0 + 427);
    const auto *dpi0_428 = buffer.data(dpi0 + 428);
    const auto *dpi0_429 = buffer.data(dpi0 + 429);
    const auto *dpi0_430 = buffer.data(dpi0 + 430);
    const auto *dpi0_431 = buffer.data(dpi0 + 431);
    const auto *dpi0_432 = buffer.data(dpi0 + 432);
    const auto *dpi0_433 = buffer.data(dpi0 + 433);
    const auto *dpi0_434 = buffer.data(dpi0 + 434);
    const auto *dpi0_447 = buffer.data(dpi0 + 447);

    const auto *dph_204 = buffer.data(dph + 204);
    const auto *dph_225 = buffer.data(dph + 225);
    const auto *dph_226 = buffer.data(dph + 226);
    const auto *dph_227 = buffer.data(dph + 227);
    const auto *dph_228 = buffer.data(dph + 228);
    const auto *dph_230 = buffer.data(dph + 230);
    const auto *dph_246 = buffer.data(dph + 246);
    const auto *dph_267 = buffer.data(dph + 267);
    const auto *dph_272 = buffer.data(dph + 272);
    const auto *dph_293 = buffer.data(dph + 293);
    const auto *dph_309 = buffer.data(dph + 309);
    const auto *dph_311 = buffer.data(dph + 311);
    const auto *dph_312 = buffer.data(dph + 312);
    const auto *dph_313 = buffer.data(dph + 313);
    const auto *dph_314 = buffer.data(dph + 314);
    const auto *dph_315 = buffer.data(dph + 315);
    const auto *dph_316 = buffer.data(dph + 316);
    const auto *dph_317 = buffer.data(dph + 317);
    const auto *dph_318 = buffer.data(dph + 318);
    const auto *dph_319 = buffer.data(dph + 319);
    const auto *dph_320 = buffer.data(dph + 320);
    const auto *dph_321 = buffer.data(dph + 321);
    const auto *dph_322 = buffer.data(dph + 322);
    const auto *dph_323 = buffer.data(dph + 323);
    const auto *dph_324 = buffer.data(dph + 324);
    const auto *dph_335 = buffer.data(dph + 335);

    const auto *dpi1_273 = buffer.data(dpi1 + 273);
    const auto *dpi1_280 = buffer.data(dpi1 + 280);
    const auto *dpi1_281 = buffer.data(dpi1 + 281);
    const auto *dpi1_283 = buffer.data(dpi1 + 283);
    const auto *dpi1_286 = buffer.data(dpi1 + 286);
    const auto *dpi1_290 = buffer.data(dpi1 + 290);
    const auto *dpi1_301 = buffer.data(dpi1 + 301);
    const auto *dpi1_303 = buffer.data(dpi1 + 303);
    const auto *dpi1_304 = buffer.data(dpi1 + 304);
    const auto *dpi1_305 = buffer.data(dpi1 + 305);
    const auto *dpi1_419 = buffer.data(dpi1 + 419);
    const auto *dpi1_420 = buffer.data(dpi1 + 420);
    const auto *dpi1_421 = buffer.data(dpi1 + 421);
    const auto *dpi1_422 = buffer.data(dpi1 + 422);
    const auto *dpi1_423 = buffer.data(dpi1 + 423);
    const auto *dpi1_424 = buffer.data(dpi1 + 424);
    const auto *dpi1_425 = buffer.data(dpi1 + 425);
    const auto *dpi1_426 = buffer.data(dpi1 + 426);
    const auto *dpi1_427 = buffer.data(dpi1 + 427);
    const auto *dpi1_428 = buffer.data(dpi1 + 428);
    const auto *dpi1_429 = buffer.data(dpi1 + 429);
    const auto *dpi1_430 = buffer.data(dpi1 + 430);
    const auto *dpi1_431 = buffer.data(dpi1 + 431);
    const auto *dpi1_432 = buffer.data(dpi1 + 432);
    const auto *dpi1_433 = buffer.data(dpi1 + 433);
    const auto *dpi1_434 = buffer.data(dpi1 + 434);
    const auto *dpi1_447 = buffer.data(dpi1 + 447);

    const auto *fsi0_219 = buffer.data(fsi0 + 219);
    const auto *fsi0_220 = buffer.data(fsi0 + 220);
    const auto *fsi0_221 = buffer.data(fsi0 + 221);
    const auto *fsi0_223 = buffer.data(fsi0 + 223);
    const auto *fsi0_245 = buffer.data(fsi0 + 245);
    const auto *fsi0_247 = buffer.data(fsi0 + 247);
    const auto *fsi0_248 = buffer.data(fsi0 + 248);
    const auto *fsi0_249 = buffer.data(fsi0 + 249);

    const auto *fsh_162 = buffer.data(fsh + 162);
    const auto *fsh_163 = buffer.data(fsh + 163);
    const auto *fsh_164 = buffer.data(fsh + 164);
    const auto *fsh_165 = buffer.data(fsh + 165);
    const auto *fsh_166 = buffer.data(fsh + 166);
    const auto *fsh_167 = buffer.data(fsh + 167);
    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_184 = buffer.data(fsh + 184);
    const auto *fsh_185 = buffer.data(fsh + 185);
    const auto *fsh_186 = buffer.data(fsh + 186);
    const auto *fsh_187 = buffer.data(fsh + 187);
    const auto *fsh_188 = buffer.data(fsh + 188);

    const auto *fsi1_219 = buffer.data(fsi1 + 219);
    const auto *fsi1_220 = buffer.data(fsi1 + 220);
    const auto *fsi1_221 = buffer.data(fsi1 + 221);
    const auto *fsi1_223 = buffer.data(fsi1 + 223);
    const auto *fsi1_245 = buffer.data(fsi1 + 245);
    const auto *fsi1_247 = buffer.data(fsi1 + 247);
    const auto *fsi1_248 = buffer.data(fsi1 + 248);
    const auto *fsi1_249 = buffer.data(fsi1 + 249);

    const auto *fpg0_332 = buffer.data(fpg0 + 332);
    const auto *fpg0_334 = buffer.data(fpg0 + 334);
    const auto *fpg0_335 = buffer.data(fpg0 + 335);
    const auto *fpg0_337 = buffer.data(fpg0 + 337);
    const auto *fpg0_338 = buffer.data(fpg0 + 338);
    const auto *fpg0_339 = buffer.data(fpg0 + 339);
    const auto *fpg0_341 = buffer.data(fpg0 + 341);
    const auto *fpg0_342 = buffer.data(fpg0 + 342);
    const auto *fpg0_343 = buffer.data(fpg0 + 343);
    const auto *fpg0_344 = buffer.data(fpg0 + 344);
    const auto *fpg0_345 = buffer.data(fpg0 + 345);
    const auto *fpg0_346 = buffer.data(fpg0 + 346);
    const auto *fpg0_347 = buffer.data(fpg0 + 347);
    const auto *fpg0_348 = buffer.data(fpg0 + 348);
    const auto *fpg0_349 = buffer.data(fpg0 + 349);
    const auto *fpg0_350 = buffer.data(fpg0 + 350);
    const auto *fpg0_351 = buffer.data(fpg0 + 351);
    const auto *fpg0_352 = buffer.data(fpg0 + 352);
    const auto *fpg0_353 = buffer.data(fpg0 + 353);
    const auto *fpg0_354 = buffer.data(fpg0 + 354);
    const auto *fpg0_355 = buffer.data(fpg0 + 355);
    const auto *fpg0_356 = buffer.data(fpg0 + 356);
    const auto *fpg0_357 = buffer.data(fpg0 + 357);
    const auto *fpg0_358 = buffer.data(fpg0 + 358);
    const auto *fpg0_359 = buffer.data(fpg0 + 359);
    const auto *fpg0_375 = buffer.data(fpg0 + 375);
    const auto *fpg0_376 = buffer.data(fpg0 + 376);
    const auto *fpg0_377 = buffer.data(fpg0 + 377);
    const auto *fpg0_378 = buffer.data(fpg0 + 378);
    const auto *fpg0_379 = buffer.data(fpg0 + 379);
    const auto *fpg0_380 = buffer.data(fpg0 + 380);
    const auto *fpg0_381 = buffer.data(fpg0 + 381);
    const auto *fpg0_382 = buffer.data(fpg0 + 382);
    const auto *fpg0_383 = buffer.data(fpg0 + 383);
    const auto *fpg0_384 = buffer.data(fpg0 + 384);
    const auto *fpg0_385 = buffer.data(fpg0 + 385);
    const auto *fpg0_386 = buffer.data(fpg0 + 386);
    const auto *fpg0_387 = buffer.data(fpg0 + 387);
    const auto *fpg0_388 = buffer.data(fpg0 + 388);
    const auto *fpg0_389 = buffer.data(fpg0 + 389);

    const auto *fpg1_332 = buffer.data(fpg1 + 332);
    const auto *fpg1_334 = buffer.data(fpg1 + 334);
    const auto *fpg1_335 = buffer.data(fpg1 + 335);
    const auto *fpg1_337 = buffer.data(fpg1 + 337);
    const auto *fpg1_338 = buffer.data(fpg1 + 338);
    const auto *fpg1_339 = buffer.data(fpg1 + 339);
    const auto *fpg1_341 = buffer.data(fpg1 + 341);
    const auto *fpg1_342 = buffer.data(fpg1 + 342);
    const auto *fpg1_343 = buffer.data(fpg1 + 343);
    const auto *fpg1_344 = buffer.data(fpg1 + 344);
    const auto *fpg1_345 = buffer.data(fpg1 + 345);
    const auto *fpg1_346 = buffer.data(fpg1 + 346);
    const auto *fpg1_347 = buffer.data(fpg1 + 347);
    const auto *fpg1_348 = buffer.data(fpg1 + 348);
    const auto *fpg1_349 = buffer.data(fpg1 + 349);
    const auto *fpg1_350 = buffer.data(fpg1 + 350);
    const auto *fpg1_351 = buffer.data(fpg1 + 351);
    const auto *fpg1_352 = buffer.data(fpg1 + 352);
    const auto *fpg1_353 = buffer.data(fpg1 + 353);
    const auto *fpg1_354 = buffer.data(fpg1 + 354);
    const auto *fpg1_355 = buffer.data(fpg1 + 355);
    const auto *fpg1_356 = buffer.data(fpg1 + 356);
    const auto *fpg1_357 = buffer.data(fpg1 + 357);
    const auto *fpg1_358 = buffer.data(fpg1 + 358);
    const auto *fpg1_359 = buffer.data(fpg1 + 359);
    const auto *fpg1_375 = buffer.data(fpg1 + 375);
    const auto *fpg1_376 = buffer.data(fpg1 + 376);
    const auto *fpg1_377 = buffer.data(fpg1 + 377);
    const auto *fpg1_378 = buffer.data(fpg1 + 378);
    const auto *fpg1_379 = buffer.data(fpg1 + 379);
    const auto *fpg1_380 = buffer.data(fpg1 + 380);
    const auto *fpg1_381 = buffer.data(fpg1 + 381);
    const auto *fpg1_382 = buffer.data(fpg1 + 382);
    const auto *fpg1_383 = buffer.data(fpg1 + 383);
    const auto *fpg1_384 = buffer.data(fpg1 + 384);
    const auto *fpg1_385 = buffer.data(fpg1 + 385);
    const auto *fpg1_386 = buffer.data(fpg1 + 386);
    const auto *fpg1_387 = buffer.data(fpg1 + 387);
    const auto *fpg1_388 = buffer.data(fpg1 + 388);
    const auto *fpg1_389 = buffer.data(fpg1 + 389);

    const auto *fph_456 = buffer.data(fph + 456);
    const auto *fph_457 = buffer.data(fph + 457);
    const auto *fph_458 = buffer.data(fph + 458);
    const auto *fph_459 = buffer.data(fph + 459);
    const auto *fph_460 = buffer.data(fph + 460);
    const auto *fph_461 = buffer.data(fph + 461);
    const auto *fph_464 = buffer.data(fph + 464);
    const auto *fph_466 = buffer.data(fph + 466);
    const auto *fph_467 = buffer.data(fph + 467);
    const auto *fph_469 = buffer.data(fph + 469);
    const auto *fph_470 = buffer.data(fph + 470);
    const auto *fph_471 = buffer.data(fph + 471);
    const auto *fph_473 = buffer.data(fph + 473);
    const auto *fph_474 = buffer.data(fph + 474);
    const auto *fph_475 = buffer.data(fph + 475);
    const auto *fph_476 = buffer.data(fph + 476);
    const auto *fph_477 = buffer.data(fph + 477);
    const auto *fph_478 = buffer.data(fph + 478);
    const auto *fph_479 = buffer.data(fph + 479);
    const auto *fph_480 = buffer.data(fph + 480);
    const auto *fph_481 = buffer.data(fph + 481);
    const auto *fph_482 = buffer.data(fph + 482);
    const auto *fph_483 = buffer.data(fph + 483);
    const auto *fph_484 = buffer.data(fph + 484);
    const auto *fph_485 = buffer.data(fph + 485);
    const auto *fph_486 = buffer.data(fph + 486);
    const auto *fph_487 = buffer.data(fph + 487);
    const auto *fph_488 = buffer.data(fph + 488);
    const auto *fph_489 = buffer.data(fph + 489);
    const auto *fph_490 = buffer.data(fph + 490);
    const auto *fph_491 = buffer.data(fph + 491);
    const auto *fph_492 = buffer.data(fph + 492);
    const auto *fph_493 = buffer.data(fph + 493);
    const auto *fph_494 = buffer.data(fph + 494);
    const auto *fph_495 = buffer.data(fph + 495);
    const auto *fph_496 = buffer.data(fph + 496);
    const auto *fph_497 = buffer.data(fph + 497);
    const auto *fph_498 = buffer.data(fph + 498);
    const auto *fph_499 = buffer.data(fph + 499);
    const auto *fph_500 = buffer.data(fph + 500);
    const auto *fph_501 = buffer.data(fph + 501);
    const auto *fph_502 = buffer.data(fph + 502);
    const auto *fph_503 = buffer.data(fph + 503);
    const auto *fph_519 = buffer.data(fph + 519);
    const auto *fph_520 = buffer.data(fph + 520);
    const auto *fph_521 = buffer.data(fph + 521);
    const auto *fph_522 = buffer.data(fph + 522);
    const auto *fph_523 = buffer.data(fph + 523);
    const auto *fph_524 = buffer.data(fph + 524);
    const auto *fph_525 = buffer.data(fph + 525);
    const auto *fph_526 = buffer.data(fph + 526);
    const auto *fph_527 = buffer.data(fph + 527);
    const auto *fph_528 = buffer.data(fph + 528);
    const auto *fph_529 = buffer.data(fph + 529);
    const auto *fph_530 = buffer.data(fph + 530);
    const auto *fph_531 = buffer.data(fph + 531);
    const auto *fph_532 = buffer.data(fph + 532);
    const auto *fph_533 = buffer.data(fph + 533);
    const auto *fph_534 = buffer.data(fph + 534);
    const auto *fph_535 = buffer.data(fph + 535);
    const auto *fph_536 = buffer.data(fph + 536);
    const auto *fph_537 = buffer.data(fph + 537);
    const auto *fph_538 = buffer.data(fph + 538);
    const auto *fph_539 = buffer.data(fph + 539);
    const auto *fph_540 = buffer.data(fph + 540);
    const auto *fph_541 = buffer.data(fph + 541);
    const auto *fph_542 = buffer.data(fph + 542);
    const auto *fph_543 = buffer.data(fph + 543);
    const auto *fph_544 = buffer.data(fph + 544);
    const auto *fph_545 = buffer.data(fph + 545);

#pragma omp simd aligned(t_603, t_604, t_605, t_606, t_607, pc_x, fsh_162, fsh_163, fsh_164, \
                         fsh_165, fsh_166, fph_456, fph_457, fph_458, fph_459, \
                         fph_460 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_603[k] = f_1 * fsh_162[k]
                   + f_4 * pc_x[k] * fph_456[k];

        t_604[k] = f_1 * fsh_163[k]
                   + f_4 * pc_x[k] * fph_457[k];

        t_605[k] = f_1 * fsh_164[k]
                   + f_4 * pc_x[k] * fph_458[k];

        t_606[k] = f_1 * fsh_165[k]
                   + f_4 * pc_x[k] * fph_459[k];

        t_607[k] = f_1 * fsh_166[k]
                   + f_4 * pc_x[k] * fph_460[k];
    }

#pragma omp simd aligned(t_608, t_609, t_610, t_611, pa_z, pb_x, pc_x, pc_z, dpi0_273, \
                         dph_204, dpi1_273, fsi0_219, fsh_167, fsi1_219, fph_456, \
                         fph_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_608[k] = f_1 * fsh_167[k]
                   + f_4 * pc_x[k] * fph_461[k];

        t_609[k] = pa_z[k] * dpi0_273[k]
                   - f_11 * pc_z[k] * dpi1_273[k];

        t_610[k] = f_1 * dph_204[k]
                   + f_4 * pc_z[k] * fph_456[k];

        t_611[k] = pb_x[k] * fsi0_219[k]
                   - f_11 * pc_x[k] * fsi1_219[k];
    }

#pragma omp simd aligned(t_612, t_613, t_614, t_615, pb_x, pc_x, pc_y, dph_272, fsi0_220, \
                         fsi0_221, fsi0_223, fsi1_220, fsi1_221, fsi1_223, \
                         fph_461 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_612[k] = pb_x[k] * fsi0_220[k]
                   - f_11 * pc_x[k] * fsi1_220[k];

        t_613[k] = pb_x[k] * fsi0_221[k]
                   - f_11 * pc_x[k] * fsi1_221[k];

        t_614[k] = f_12 * dph_272[k]
                   + f_4 * pc_y[k] * fph_461[k];

        t_615[k] = pb_x[k] * fsi0_223[k]
                   - f_11 * pc_x[k] * fsi1_223[k];
    }

#pragma omp simd aligned(t_616, t_617, t_618, t_619, pa_z, pc_x, pc_z, dpi0_280, dpi0_281, \
                         dpi0_283, dpi1_280, dpi1_281, dpi1_283, fpg0_332, fpg1_332, \
                         fph_464 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_616[k] = pa_z[k] * dpi0_280[k]
                   - f_11 * pc_z[k] * dpi1_280[k];

        t_617[k] = pa_z[k] * dpi0_281[k]
                   - f_11 * pc_z[k] * dpi1_281[k];

        t_618[k] = f_17 * fpg0_332[k]
                   - f_18 * fpg1_332[k]
                   + f_4 * pc_x[k] * fph_464[k];

        t_619[k] = pa_z[k] * dpi0_283[k]
                   - f_11 * pc_z[k] * dpi1_283[k];
    }

#pragma omp simd aligned(t_620, t_621, t_622, pa_z, pc_x, pc_z, dpi0_286, dpi1_286, fpg0_334, \
                         fpg0_335, fpg1_334, fpg1_335, fph_466, \
                         fph_467 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_620[k] = f_9 * fpg0_334[k]
                   - f_10 * fpg1_334[k]
                   + f_4 * pc_x[k] * fph_466[k];

        t_621[k] = f_9 * fpg0_335[k]
                   - f_10 * fpg1_335[k]
                   + f_4 * pc_x[k] * fph_467[k];

        t_622[k] = pa_z[k] * dpi0_286[k]
                   - f_11 * pc_z[k] * dpi1_286[k];
    }

#pragma omp simd aligned(t_623, t_624, t_625, pc_x, fpg0_337, fpg0_338, fpg0_339, fpg1_337, \
                         fpg1_338, fpg1_339, fph_469, fph_470, \
                         fph_471 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_623[k] = f_7 * fpg0_337[k]
                   - f_8 * fpg1_337[k]
                   + f_4 * pc_x[k] * fph_469[k];

        t_624[k] = f_7 * fpg0_338[k]
                   - f_8 * fpg1_338[k]
                   + f_4 * pc_x[k] * fph_470[k];

        t_625[k] = f_7 * fpg0_339[k]
                   - f_8 * fpg1_339[k]
                   + f_4 * pc_x[k] * fph_471[k];
    }

#pragma omp simd aligned(t_626, t_627, t_628, pa_z, pc_x, pc_z, dpi0_290, dpi1_290, fpg0_341, \
                         fpg0_342, fpg1_341, fpg1_342, fph_473, \
                         fph_474 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_626[k] = pa_z[k] * dpi0_290[k]
                   - f_11 * pc_z[k] * dpi1_290[k];

        t_627[k] = f_5 * fpg0_341[k]
                   - f_6 * fpg1_341[k]
                   + f_4 * pc_x[k] * fph_473[k];

        t_628[k] = f_5 * fpg0_342[k]
                   - f_6 * fpg1_342[k]
                   + f_4 * pc_x[k] * fph_474[k];
    }

#pragma omp simd aligned(t_629, t_630, t_631, t_632, t_633, pc_x, fpg0_343, fpg0_344, \
                         fpg1_343, fpg1_344, fph_475, fph_476, fph_477, fph_478, \
                         fph_479 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_629[k] = f_5 * fpg0_343[k]
                   - f_6 * fpg1_343[k]
                   + f_4 * pc_x[k] * fph_475[k];

        t_630[k] = f_5 * fpg0_344[k]
                   - f_6 * fpg1_344[k]
                   + f_4 * pc_x[k] * fph_476[k];

        t_631[k] = f_4 * pc_x[k] * fph_477[k];

        t_632[k] = f_4 * pc_x[k] * fph_478[k];

        t_633[k] = f_4 * pc_x[k] * fph_479[k];
    }

#pragma omp simd aligned(t_634, t_635, t_636, t_637, t_638, pa_z, pc_x, pc_z, dpi0_301, \
                         dph_225, dpi1_301, fph_477, fph_480, fph_481, \
                         fph_482 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_634[k] = f_4 * pc_x[k] * fph_480[k];

        t_635[k] = f_4 * pc_x[k] * fph_481[k];

        t_636[k] = f_4 * pc_x[k] * fph_482[k];

        t_637[k] = pa_z[k] * dpi0_301[k]
                   - f_11 * pc_z[k] * dpi1_301[k];

        t_638[k] = f_1 * dph_225[k]
                   + f_4 * pc_z[k] * fph_477[k];
    }

#pragma omp simd aligned(t_639, t_640, t_641, pa_z, pc_z, dpi0_303, dpi0_304, dpi0_305, \
                         dph_226, dph_227, dph_228, dpi1_303, dpi1_304, \
                         dpi1_305 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_639[k] = pa_z[k] * dpi0_303[k]
                   + f_12 * dph_226[k]
                   - f_11 * pc_z[k] * dpi1_303[k];

        t_640[k] = pa_z[k] * dpi0_304[k]
                   + f_0 * dph_227[k]
                   - f_11 * pc_z[k] * dpi1_304[k];

        t_641[k] = pa_z[k] * dpi0_305[k]
                   + f_13 * dph_228[k]
                   - f_11 * pc_z[k] * dpi1_305[k];
    }

#pragma omp simd aligned(t_642, t_643, t_644, pc_x, pc_y, pc_z, dph_230, dph_293, fsh_167, \
                         fpg0_344, fpg0_345, fpg1_344, fpg1_345, fph_482, \
                         fph_483 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_642[k] = f_12 * dph_293[k]
                   + f_1 * fsh_167[k]
                   + f_4 * pc_y[k] * fph_482[k];

        t_643[k] = f_1 * dph_230[k]
                   + f_2 * fpg0_344[k]
                   - f_3 * fpg1_344[k]
                   + f_4 * pc_z[k] * fph_482[k];

        t_644[k] = f_2 * fpg0_345[k]
                   - f_3 * fpg1_345[k]
                   + f_4 * pc_x[k] * fph_483[k];
    }

#pragma omp simd aligned(t_645, t_646, t_647, pc_x, fpg0_346, fpg0_347, fpg0_348, fpg1_346, \
                         fpg1_347, fpg1_348, fph_484, fph_485, \
                         fph_486 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_645[k] = f_17 * fpg0_346[k]
                   - f_18 * fpg1_346[k]
                   + f_4 * pc_x[k] * fph_484[k];

        t_646[k] = f_17 * fpg0_347[k]
                   - f_18 * fpg1_347[k]
                   + f_4 * pc_x[k] * fph_485[k];

        t_647[k] = f_9 * fpg0_348[k]
                   - f_10 * fpg1_348[k]
                   + f_4 * pc_x[k] * fph_486[k];
    }

#pragma omp simd aligned(t_648, t_649, t_650, pc_x, fpg0_349, fpg0_350, fpg0_351, fpg1_349, \
                         fpg1_350, fpg1_351, fph_487, fph_488, \
                         fph_489 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_648[k] = f_9 * fpg0_349[k]
                   - f_10 * fpg1_349[k]
                   + f_4 * pc_x[k] * fph_487[k];

        t_649[k] = f_9 * fpg0_350[k]
                   - f_10 * fpg1_350[k]
                   + f_4 * pc_x[k] * fph_488[k];

        t_650[k] = f_7 * fpg0_351[k]
                   - f_8 * fpg1_351[k]
                   + f_4 * pc_x[k] * fph_489[k];
    }

#pragma omp simd aligned(t_651, t_652, t_653, pc_x, fpg0_352, fpg0_353, fpg0_354, fpg1_352, \
                         fpg1_353, fpg1_354, fph_490, fph_491, \
                         fph_492 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_651[k] = f_7 * fpg0_352[k]
                   - f_8 * fpg1_352[k]
                   + f_4 * pc_x[k] * fph_490[k];

        t_652[k] = f_7 * fpg0_353[k]
                   - f_8 * fpg1_353[k]
                   + f_4 * pc_x[k] * fph_491[k];

        t_653[k] = f_7 * fpg0_354[k]
                   - f_8 * fpg1_354[k]
                   + f_4 * pc_x[k] * fph_492[k];
    }

#pragma omp simd aligned(t_654, t_655, t_656, pc_x, fpg0_355, fpg0_356, fpg0_357, fpg1_355, \
                         fpg1_356, fpg1_357, fph_493, fph_494, \
                         fph_495 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_654[k] = f_5 * fpg0_355[k]
                   - f_6 * fpg1_355[k]
                   + f_4 * pc_x[k] * fph_493[k];

        t_655[k] = f_5 * fpg0_356[k]
                   - f_6 * fpg1_356[k]
                   + f_4 * pc_x[k] * fph_494[k];

        t_656[k] = f_5 * fpg0_357[k]
                   - f_6 * fpg1_357[k]
                   + f_4 * pc_x[k] * fph_495[k];
    }

#pragma omp simd aligned(t_657, t_658, t_659, t_660, t_661, pc_x, fpg0_358, fpg0_359, \
                         fpg1_358, fpg1_359, fph_496, fph_497, fph_498, fph_499, \
                         fph_500 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_657[k] = f_5 * fpg0_358[k]
                   - f_6 * fpg1_358[k]
                   + f_4 * pc_x[k] * fph_496[k];

        t_658[k] = f_5 * fpg0_359[k]
                   - f_6 * fpg1_359[k]
                   + f_4 * pc_x[k] * fph_497[k];

        t_659[k] = f_4 * pc_x[k] * fph_498[k];

        t_660[k] = f_4 * pc_x[k] * fph_499[k];

        t_661[k] = f_4 * pc_x[k] * fph_500[k];
    }

#pragma omp simd aligned(t_662, t_663, t_664, t_665, pc_x, pc_y, dph_309, fpg0_355, fpg1_355, \
                         fph_498, fph_501, fph_502, fph_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_662[k] = f_4 * pc_x[k] * fph_501[k];

        t_663[k] = f_4 * pc_x[k] * fph_502[k];

        t_664[k] = f_4 * pc_x[k] * fph_503[k];

        t_665[k] = f_12 * dph_309[k]
                   + f_2 * fpg0_355[k]
                   - f_3 * fpg1_355[k]
                   + f_4 * pc_y[k] * fph_498[k];
    }

#pragma omp simd aligned(t_666, t_667, t_668, pc_y, pc_z, dph_246, dph_311, dph_312, fsh_162, \
                         fpg0_357, fpg0_358, fpg1_357, fpg1_358, fph_498, fph_500, \
                         fph_501 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_666[k] = f_1 * dph_246[k]
                   + f_1 * fsh_162[k]
                   + f_4 * pc_z[k] * fph_498[k];

        t_667[k] = f_12 * dph_311[k]
                   + f_9 * fpg0_357[k]
                   - f_10 * fpg1_357[k]
                   + f_4 * pc_y[k] * fph_500[k];

        t_668[k] = f_12 * dph_312[k]
                   + f_7 * fpg0_358[k]
                   - f_8 * fpg1_358[k]
                   + f_4 * pc_y[k] * fph_501[k];
    }

#pragma omp simd aligned(t_669, t_670, t_671, pa_y, pc_y, ppi0_251, ppi1_251, dpi0_419, \
                         dph_313, dph_314, dpi1_419, fpg0_359, fpg1_359, fph_502, \
                         fph_503 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_669[k] = f_12 * dph_313[k]
                   + f_5 * fpg0_359[k]
                   - f_6 * fpg1_359[k]
                   + f_4 * pc_y[k] * fph_502[k];

        t_670[k] = f_12 * dph_314[k]
                   + f_4 * pc_y[k] * fph_503[k];

        t_671[k] = f_15 * ppi0_251[k]
                   - f_16 * ppi1_251[k]
                   + pa_y[k] * dpi0_419[k]
                   - f_11 * pc_y[k] * dpi1_419[k];
    }

#pragma omp simd aligned(t_672, t_673, t_674, t_675, pa_y, pc_y, dpi0_420, dpi0_421, dpi0_422, \
                         dpi0_423, dph_315, dph_316, dpi1_420, dpi1_421, dpi1_422, \
                         dpi1_423 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_672[k] = pa_y[k] * dpi0_420[k]
                   - f_11 * pc_y[k] * dpi1_420[k];

        t_673[k] = pa_y[k] * dpi0_421[k]
                   + f_1 * dph_315[k]
                   - f_11 * pc_y[k] * dpi1_421[k];

        t_674[k] = pa_y[k] * dpi0_422[k]
                   - f_11 * pc_y[k] * dpi1_422[k];

        t_675[k] = pa_y[k] * dpi0_423[k]
                   + f_12 * dph_316[k]
                   - f_11 * pc_y[k] * dpi1_423[k];
    }

#pragma omp simd aligned(t_676, t_677, t_678, pa_y, pc_y, dpi0_424, dpi0_425, dpi0_426, \
                         dph_317, dph_318, dpi1_424, dpi1_425, \
                         dpi1_426 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_676[k] = pa_y[k] * dpi0_424[k]
                   + f_1 * dph_317[k]
                   - f_11 * pc_y[k] * dpi1_424[k];

        t_677[k] = pa_y[k] * dpi0_425[k]
                   - f_11 * pc_y[k] * dpi1_425[k];

        t_678[k] = pa_y[k] * dpi0_426[k]
                   + f_0 * dph_318[k]
                   - f_11 * pc_y[k] * dpi1_426[k];
    }

#pragma omp simd aligned(t_679, t_680, t_681, pa_y, pc_y, dpi0_427, dpi0_428, dpi0_429, \
                         dph_319, dph_320, dpi1_427, dpi1_428, \
                         dpi1_429 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_679[k] = pa_y[k] * dpi0_427[k]
                   + f_12 * dph_319[k]
                   - f_11 * pc_y[k] * dpi1_427[k];

        t_680[k] = pa_y[k] * dpi0_428[k]
                   + f_1 * dph_320[k]
                   - f_11 * pc_y[k] * dpi1_428[k];

        t_681[k] = pa_y[k] * dpi0_429[k]
                   - f_11 * pc_y[k] * dpi1_429[k];
    }

#pragma omp simd aligned(t_682, t_683, t_684, pa_y, pc_y, dpi0_430, dpi0_431, dpi0_432, \
                         dph_321, dph_322, dph_323, dpi1_430, dpi1_431, \
                         dpi1_432 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_682[k] = pa_y[k] * dpi0_430[k]
                   + f_13 * dph_321[k]
                   - f_11 * pc_y[k] * dpi1_430[k];

        t_683[k] = pa_y[k] * dpi0_431[k]
                   + f_0 * dph_322[k]
                   - f_11 * pc_y[k] * dpi1_431[k];

        t_684[k] = pa_y[k] * dpi0_432[k]
                   + f_12 * dph_323[k]
                   - f_11 * pc_y[k] * dpi1_432[k];
    }

#pragma omp simd aligned(t_685, t_686, t_687, t_688, pa_y, pc_x, pc_y, dpi0_433, dpi0_434, \
                         dph_324, dpi1_433, dpi1_434, fsh_183, fsh_184, fph_519, \
                         fph_520 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_685[k] = pa_y[k] * dpi0_433[k]
                   + f_1 * dph_324[k]
                   - f_11 * pc_y[k] * dpi1_433[k];

        t_686[k] = pa_y[k] * dpi0_434[k]
                   - f_11 * pc_y[k] * dpi1_434[k];

        t_687[k] = f_1 * fsh_183[k]
                   + f_4 * pc_x[k] * fph_519[k];

        t_688[k] = f_1 * fsh_184[k]
                   + f_4 * pc_x[k] * fph_520[k];
    }

#pragma omp simd aligned(t_689, t_690, t_691, t_692, pc_x, fsh_185, fsh_186, fsh_187, fsh_188, \
                         fph_521, fph_522, fph_523, fph_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_689[k] = f_1 * fsh_185[k]
                   + f_4 * pc_x[k] * fph_521[k];

        t_690[k] = f_1 * fsh_186[k]
                   + f_4 * pc_x[k] * fph_522[k];

        t_691[k] = f_1 * fsh_187[k]
                   + f_4 * pc_x[k] * fph_523[k];

        t_692[k] = f_1 * fsh_188[k]
                   + f_4 * pc_x[k] * fph_524[k];
    }

#pragma omp simd aligned(t_693, t_694, t_695, t_696, pb_x, pc_x, pc_z, dph_267, fsi0_245, \
                         fsi0_247, fsi0_248, fsi1_245, fsi1_247, fsi1_248, \
                         fph_519 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_693[k] = pb_x[k] * fsi0_245[k]
                   - f_11 * pc_x[k] * fsi1_245[k];

        t_694[k] = f_12 * dph_267[k]
                   + f_4 * pc_z[k] * fph_519[k];

        t_695[k] = pb_x[k] * fsi0_247[k]
                   - f_11 * pc_x[k] * fsi1_247[k];

        t_696[k] = pb_x[k] * fsi0_248[k]
                   - f_11 * pc_x[k] * fsi1_248[k];
    }

#pragma omp simd aligned(t_697, t_698, t_699, pa_y, pb_x, pc_x, pc_y, dpi0_447, dph_335, \
                         dpi1_447, fsi0_249, fsi1_249, fph_524 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_697[k] = pb_x[k] * fsi0_249[k]
                   - f_11 * pc_x[k] * fsi1_249[k];

        t_698[k] = f_1 * dph_335[k]
                   + f_4 * pc_y[k] * fph_524[k];

        t_699[k] = pa_y[k] * dpi0_447[k]
                   - f_11 * pc_y[k] * dpi1_447[k];
    }

#pragma omp simd aligned(t_700, t_701, t_702, pc_x, fpg0_375, fpg0_376, fpg0_377, fpg1_375, \
                         fpg1_376, fpg1_377, fph_525, fph_526, \
                         fph_527 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_700[k] = f_2 * fpg0_375[k]
                   - f_3 * fpg1_375[k]
                   + f_4 * pc_x[k] * fph_525[k];

        t_701[k] = f_17 * fpg0_376[k]
                   - f_18 * fpg1_376[k]
                   + f_4 * pc_x[k] * fph_526[k];

        t_702[k] = f_17 * fpg0_377[k]
                   - f_18 * fpg1_377[k]
                   + f_4 * pc_x[k] * fph_527[k];
    }

#pragma omp simd aligned(t_703, t_704, t_705, pc_x, fpg0_378, fpg0_379, fpg0_380, fpg1_378, \
                         fpg1_379, fpg1_380, fph_528, fph_529, \
                         fph_530 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_703[k] = f_9 * fpg0_378[k]
                   - f_10 * fpg1_378[k]
                   + f_4 * pc_x[k] * fph_528[k];

        t_704[k] = f_9 * fpg0_379[k]
                   - f_10 * fpg1_379[k]
                   + f_4 * pc_x[k] * fph_529[k];

        t_705[k] = f_9 * fpg0_380[k]
                   - f_10 * fpg1_380[k]
                   + f_4 * pc_x[k] * fph_530[k];
    }

#pragma omp simd aligned(t_706, t_707, t_708, pc_x, fpg0_381, fpg0_382, fpg0_383, fpg1_381, \
                         fpg1_382, fpg1_383, fph_531, fph_532, \
                         fph_533 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_706[k] = f_7 * fpg0_381[k]
                   - f_8 * fpg1_381[k]
                   + f_4 * pc_x[k] * fph_531[k];

        t_707[k] = f_7 * fpg0_382[k]
                   - f_8 * fpg1_382[k]
                   + f_4 * pc_x[k] * fph_532[k];

        t_708[k] = f_7 * fpg0_383[k]
                   - f_8 * fpg1_383[k]
                   + f_4 * pc_x[k] * fph_533[k];
    }

#pragma omp simd aligned(t_709, t_710, t_711, pc_x, fpg0_384, fpg0_385, fpg0_386, fpg1_384, \
                         fpg1_385, fpg1_386, fph_534, fph_535, \
                         fph_536 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_709[k] = f_7 * fpg0_384[k]
                   - f_8 * fpg1_384[k]
                   + f_4 * pc_x[k] * fph_534[k];

        t_710[k] = f_5 * fpg0_385[k]
                   - f_6 * fpg1_385[k]
                   + f_4 * pc_x[k] * fph_535[k];

        t_711[k] = f_5 * fpg0_386[k]
                   - f_6 * fpg1_386[k]
                   + f_4 * pc_x[k] * fph_536[k];
    }

#pragma omp simd aligned(t_712, t_713, t_714, t_715, pc_x, fpg0_387, fpg0_388, fpg0_389, \
                         fpg1_387, fpg1_388, fpg1_389, fph_537, fph_538, fph_539, \
                         fph_540 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_712[k] = f_5 * fpg0_387[k]
                   - f_6 * fpg1_387[k]
                   + f_4 * pc_x[k] * fph_537[k];

        t_713[k] = f_5 * fpg0_388[k]
                   - f_6 * fpg1_388[k]
                   + f_4 * pc_x[k] * fph_538[k];

        t_714[k] = f_5 * fpg0_389[k]
                   - f_6 * fpg1_389[k]
                   + f_4 * pc_x[k] * fph_539[k];

        t_715[k] = f_4 * pc_x[k] * fph_540[k];
    }

#pragma omp simd aligned(t_716, t_717, t_718, t_719, t_720, pc_x, fph_541, fph_542, fph_543, \
                         fph_544, fph_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_716[k] = f_4 * pc_x[k] * fph_541[k];

        t_717[k] = f_4 * pc_x[k] * fph_542[k];

        t_718[k] = f_4 * pc_x[k] * fph_543[k];

        t_719[k] = f_4 * pc_x[k] * fph_544[k];

        t_720[k] = f_4 * pc_x[k] * fph_545[k];
    }
}

static auto
compute_prim_fpi_three_center_electron_repulsion_0_piece6(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t dpi0, const size_t dph,
                                                          const size_t dpi1, const size_t fsi0,
                                                          const size_t fsh, const size_t fsi1,
                                                          const size_t fpg0, const size_t fpg1,
                                                          const size_t fph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.5 / q;
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
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_17 = 2.0 / gamma;
    const auto f_18 = 2.0 * p / (gamma * q);
    const auto f_19 = 2.5 / q;

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
    auto *t_836 = buffer.data(target + 836);
    auto *t_837 = buffer.data(target + 837);
    auto *t_838 = buffer.data(target + 838);
    auto *t_839 = buffer.data(target + 839);

    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *dpi0_476 = buffer.data(dpi0 + 476);
    const auto *dpi0_478 = buffer.data(dpi0 + 478);
    const auto *dpi0_481 = buffer.data(dpi0 + 481);
    const auto *dpi0_485 = buffer.data(dpi0 + 485);
    const auto *dpi0_490 = buffer.data(dpi0 + 490);
    const auto *dpi0_497 = buffer.data(dpi0 + 497);
    const auto *dpi0_499 = buffer.data(dpi0 + 499);
    const auto *dpi0_500 = buffer.data(dpi0 + 500);
    const auto *dpi0_501 = buffer.data(dpi0 + 501);
    const auto *dpi0_503 = buffer.data(dpi0 + 503);

    const auto *dph_288 = buffer.data(dph + 288);
    const auto *dph_293 = buffer.data(dph + 293);
    const auto *dph_309 = buffer.data(dph + 309);
    const auto *dph_351 = buffer.data(dph + 351);
    const auto *dph_353 = buffer.data(dph + 353);
    const auto *dph_354 = buffer.data(dph + 354);
    const auto *dph_355 = buffer.data(dph + 355);
    const auto *dph_356 = buffer.data(dph + 356);
    const auto *dph_372 = buffer.data(dph + 372);
    const auto *dph_374 = buffer.data(dph + 374);
    const auto *dph_375 = buffer.data(dph + 375);
    const auto *dph_376 = buffer.data(dph + 376);
    const auto *dph_377 = buffer.data(dph + 377);

    const auto *dpi1_476 = buffer.data(dpi1 + 476);
    const auto *dpi1_478 = buffer.data(dpi1 + 478);
    const auto *dpi1_481 = buffer.data(dpi1 + 481);
    const auto *dpi1_485 = buffer.data(dpi1 + 485);
    const auto *dpi1_490 = buffer.data(dpi1 + 490);
    const auto *dpi1_497 = buffer.data(dpi1 + 497);
    const auto *dpi1_499 = buffer.data(dpi1 + 499);
    const auto *dpi1_500 = buffer.data(dpi1 + 500);
    const auto *dpi1_501 = buffer.data(dpi1 + 501);
    const auto *dpi1_503 = buffer.data(dpi1 + 503);

    const auto *fsi0_252 = buffer.data(fsi0 + 252);
    const auto *fsi0_254 = buffer.data(fsi0 + 254);
    const auto *fsi0_255 = buffer.data(fsi0 + 255);
    const auto *fsi0_257 = buffer.data(fsi0 + 257);
    const auto *fsi0_258 = buffer.data(fsi0 + 258);
    const auto *fsi0_259 = buffer.data(fsi0 + 259);
    const auto *fsi0_261 = buffer.data(fsi0 + 261);
    const auto *fsi0_262 = buffer.data(fsi0 + 262);
    const auto *fsi0_263 = buffer.data(fsi0 + 263);
    const auto *fsi0_264 = buffer.data(fsi0 + 264);
    const auto *fsi0_266 = buffer.data(fsi0 + 266);
    const auto *fsi0_273 = buffer.data(fsi0 + 273);
    const auto *fsi0_274 = buffer.data(fsi0 + 274);
    const auto *fsi0_275 = buffer.data(fsi0 + 275);
    const auto *fsi0_276 = buffer.data(fsi0 + 276);
    const auto *fsi0_277 = buffer.data(fsi0 + 277);
    const auto *fsi0_279 = buffer.data(fsi0 + 279);

    const auto *fsh_183 = buffer.data(fsh + 183);
    const auto *fsh_185 = buffer.data(fsh + 185);
    const auto *fsh_186 = buffer.data(fsh + 186);
    const auto *fsh_187 = buffer.data(fsh + 187);
    const auto *fsh_188 = buffer.data(fsh + 188);
    const auto *fsh_189 = buffer.data(fsh + 189);
    const auto *fsh_191 = buffer.data(fsh + 191);
    const auto *fsh_192 = buffer.data(fsh + 192);
    const auto *fsh_194 = buffer.data(fsh + 194);
    const auto *fsh_195 = buffer.data(fsh + 195);
    const auto *fsh_196 = buffer.data(fsh + 196);
    const auto *fsh_198 = buffer.data(fsh + 198);
    const auto *fsh_199 = buffer.data(fsh + 199);
    const auto *fsh_200 = buffer.data(fsh + 200);
    const auto *fsh_201 = buffer.data(fsh + 201);
    const auto *fsh_203 = buffer.data(fsh + 203);
    const auto *fsh_204 = buffer.data(fsh + 204);
    const auto *fsh_205 = buffer.data(fsh + 205);
    const auto *fsh_206 = buffer.data(fsh + 206);
    const auto *fsh_207 = buffer.data(fsh + 207);
    const auto *fsh_208 = buffer.data(fsh + 208);
    const auto *fsh_209 = buffer.data(fsh + 209);

    const auto *fsi1_252 = buffer.data(fsi1 + 252);
    const auto *fsi1_254 = buffer.data(fsi1 + 254);
    const auto *fsi1_255 = buffer.data(fsi1 + 255);
    const auto *fsi1_257 = buffer.data(fsi1 + 257);
    const auto *fsi1_258 = buffer.data(fsi1 + 258);
    const auto *fsi1_259 = buffer.data(fsi1 + 259);
    const auto *fsi1_261 = buffer.data(fsi1 + 261);
    const auto *fsi1_262 = buffer.data(fsi1 + 262);
    const auto *fsi1_263 = buffer.data(fsi1 + 263);
    const auto *fsi1_264 = buffer.data(fsi1 + 264);
    const auto *fsi1_266 = buffer.data(fsi1 + 266);
    const auto *fsi1_273 = buffer.data(fsi1 + 273);
    const auto *fsi1_274 = buffer.data(fsi1 + 274);
    const auto *fsi1_275 = buffer.data(fsi1 + 275);
    const auto *fsi1_276 = buffer.data(fsi1 + 276);
    const auto *fsi1_277 = buffer.data(fsi1 + 277);
    const auto *fsi1_279 = buffer.data(fsi1 + 279);

    const auto *fpg0_385 = buffer.data(fpg0 + 385);
    const auto *fpg0_387 = buffer.data(fpg0 + 387);
    const auto *fpg0_388 = buffer.data(fpg0 + 388);
    const auto *fpg0_389 = buffer.data(fpg0 + 389);
    const auto *fpg0_391 = buffer.data(fpg0 + 391);
    const auto *fpg0_393 = buffer.data(fpg0 + 393);
    const auto *fpg0_394 = buffer.data(fpg0 + 394);
    const auto *fpg0_396 = buffer.data(fpg0 + 396);
    const auto *fpg0_397 = buffer.data(fpg0 + 397);
    const auto *fpg0_398 = buffer.data(fpg0 + 398);
    const auto *fpg0_400 = buffer.data(fpg0 + 400);
    const auto *fpg0_401 = buffer.data(fpg0 + 401);
    const auto *fpg0_402 = buffer.data(fpg0 + 402);
    const auto *fpg0_403 = buffer.data(fpg0 + 403);
    const auto *fpg0_423 = buffer.data(fpg0 + 423);
    const auto *fpg0_426 = buffer.data(fpg0 + 426);
    const auto *fpg0_427 = buffer.data(fpg0 + 427);
    const auto *fpg0_430 = buffer.data(fpg0 + 430);
    const auto *fpg0_431 = buffer.data(fpg0 + 431);
    const auto *fpg0_432 = buffer.data(fpg0 + 432);
    const auto *fpg0_435 = buffer.data(fpg0 + 435);
    const auto *fpg0_437 = buffer.data(fpg0 + 437);
    const auto *fpg0_438 = buffer.data(fpg0 + 438);
    const auto *fpg0_440 = buffer.data(fpg0 + 440);
    const auto *fpg0_441 = buffer.data(fpg0 + 441);
    const auto *fpg0_442 = buffer.data(fpg0 + 442);
    const auto *fpg0_444 = buffer.data(fpg0 + 444);
    const auto *fpg0_445 = buffer.data(fpg0 + 445);
    const auto *fpg0_446 = buffer.data(fpg0 + 446);
    const auto *fpg0_447 = buffer.data(fpg0 + 447);
    const auto *fpg0_448 = buffer.data(fpg0 + 448);
    const auto *fpg0_449 = buffer.data(fpg0 + 449);

    const auto *fpg1_385 = buffer.data(fpg1 + 385);
    const auto *fpg1_387 = buffer.data(fpg1 + 387);
    const auto *fpg1_388 = buffer.data(fpg1 + 388);
    const auto *fpg1_389 = buffer.data(fpg1 + 389);
    const auto *fpg1_391 = buffer.data(fpg1 + 391);
    const auto *fpg1_393 = buffer.data(fpg1 + 393);
    const auto *fpg1_394 = buffer.data(fpg1 + 394);
    const auto *fpg1_396 = buffer.data(fpg1 + 396);
    const auto *fpg1_397 = buffer.data(fpg1 + 397);
    const auto *fpg1_398 = buffer.data(fpg1 + 398);
    const auto *fpg1_400 = buffer.data(fpg1 + 400);
    const auto *fpg1_401 = buffer.data(fpg1 + 401);
    const auto *fpg1_402 = buffer.data(fpg1 + 402);
    const auto *fpg1_403 = buffer.data(fpg1 + 403);
    const auto *fpg1_423 = buffer.data(fpg1 + 423);
    const auto *fpg1_426 = buffer.data(fpg1 + 426);
    const auto *fpg1_427 = buffer.data(fpg1 + 427);
    const auto *fpg1_430 = buffer.data(fpg1 + 430);
    const auto *fpg1_431 = buffer.data(fpg1 + 431);
    const auto *fpg1_432 = buffer.data(fpg1 + 432);
    const auto *fpg1_435 = buffer.data(fpg1 + 435);
    const auto *fpg1_437 = buffer.data(fpg1 + 437);
    const auto *fpg1_438 = buffer.data(fpg1 + 438);
    const auto *fpg1_440 = buffer.data(fpg1 + 440);
    const auto *fpg1_441 = buffer.data(fpg1 + 441);
    const auto *fpg1_442 = buffer.data(fpg1 + 442);
    const auto *fpg1_444 = buffer.data(fpg1 + 444);
    const auto *fpg1_445 = buffer.data(fpg1 + 445);
    const auto *fpg1_446 = buffer.data(fpg1 + 446);
    const auto *fpg1_447 = buffer.data(fpg1 + 447);
    const auto *fpg1_448 = buffer.data(fpg1 + 448);
    const auto *fpg1_449 = buffer.data(fpg1 + 449);

    const auto *fph_540 = buffer.data(fph + 540);
    const auto *fph_542 = buffer.data(fph + 542);
    const auto *fph_543 = buffer.data(fph + 543);
    const auto *fph_544 = buffer.data(fph + 544);
    const auto *fph_545 = buffer.data(fph + 545);
    const auto *fph_547 = buffer.data(fph + 547);
    const auto *fph_549 = buffer.data(fph + 549);
    const auto *fph_550 = buffer.data(fph + 550);
    const auto *fph_552 = buffer.data(fph + 552);
    const auto *fph_553 = buffer.data(fph + 553);
    const auto *fph_554 = buffer.data(fph + 554);
    const auto *fph_556 = buffer.data(fph + 556);
    const auto *fph_557 = buffer.data(fph + 557);
    const auto *fph_558 = buffer.data(fph + 558);
    const auto *fph_559 = buffer.data(fph + 559);
    const auto *fph_561 = buffer.data(fph + 561);
    const auto *fph_562 = buffer.data(fph + 562);
    const auto *fph_563 = buffer.data(fph + 563);
    const auto *fph_564 = buffer.data(fph + 564);
    const auto *fph_565 = buffer.data(fph + 565);
    const auto *fph_566 = buffer.data(fph + 566);
    const auto *fph_567 = buffer.data(fph + 567);
    const auto *fph_569 = buffer.data(fph + 569);
    const auto *fph_572 = buffer.data(fph + 572);
    const auto *fph_576 = buffer.data(fph + 576);
    const auto *fph_582 = buffer.data(fph + 582);
    const auto *fph_583 = buffer.data(fph + 583);
    const auto *fph_584 = buffer.data(fph + 584);
    const auto *fph_585 = buffer.data(fph + 585);
    const auto *fph_586 = buffer.data(fph + 586);
    const auto *fph_587 = buffer.data(fph + 587);
    const auto *fph_588 = buffer.data(fph + 588);
    const auto *fph_590 = buffer.data(fph + 590);
    const auto *fph_591 = buffer.data(fph + 591);
    const auto *fph_593 = buffer.data(fph + 593);
    const auto *fph_594 = buffer.data(fph + 594);
    const auto *fph_595 = buffer.data(fph + 595);
    const auto *fph_597 = buffer.data(fph + 597);
    const auto *fph_598 = buffer.data(fph + 598);
    const auto *fph_599 = buffer.data(fph + 599);
    const auto *fph_600 = buffer.data(fph + 600);
    const auto *fph_603 = buffer.data(fph + 603);
    const auto *fph_604 = buffer.data(fph + 604);
    const auto *fph_605 = buffer.data(fph + 605);
    const auto *fph_606 = buffer.data(fph + 606);
    const auto *fph_607 = buffer.data(fph + 607);
    const auto *fph_608 = buffer.data(fph + 608);
    const auto *fph_609 = buffer.data(fph + 609);
    const auto *fph_611 = buffer.data(fph + 611);
    const auto *fph_612 = buffer.data(fph + 612);
    const auto *fph_614 = buffer.data(fph + 614);
    const auto *fph_615 = buffer.data(fph + 615);
    const auto *fph_616 = buffer.data(fph + 616);
    const auto *fph_618 = buffer.data(fph + 618);
    const auto *fph_619 = buffer.data(fph + 619);
    const auto *fph_620 = buffer.data(fph + 620);
    const auto *fph_621 = buffer.data(fph + 621);
    const auto *fph_623 = buffer.data(fph + 623);
    const auto *fph_624 = buffer.data(fph + 624);
    const auto *fph_625 = buffer.data(fph + 625);
    const auto *fph_626 = buffer.data(fph + 626);
    const auto *fph_627 = buffer.data(fph + 627);
    const auto *fph_628 = buffer.data(fph + 628);
    const auto *fph_629 = buffer.data(fph + 629);

#pragma omp simd aligned(t_721, t_722, t_723, pc_y, pc_z, dph_288, dph_351, dph_353, fsh_183, \
                         fsh_185, fpg0_385, fpg0_387, fpg1_385, fpg1_387, fph_540, \
                         fph_542 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_721[k] = f_1 * dph_351[k]
                   + f_1 * fsh_183[k]
                   + f_2 * fpg0_385[k]
                   - f_3 * fpg1_385[k]
                   + f_4 * pc_y[k] * fph_540[k];

        t_722[k] = f_12 * dph_288[k]
                   + f_4 * pc_z[k] * fph_540[k];

        t_723[k] = f_1 * dph_353[k]
                   + f_1 * fsh_185[k]
                   + f_9 * fpg0_387[k]
                   - f_10 * fpg1_387[k]
                   + f_4 * pc_y[k] * fph_542[k];
    }

#pragma omp simd aligned(t_724, t_725, pc_y, dph_354, dph_355, fsh_186, fsh_187, fpg0_388, \
                         fpg0_389, fpg1_388, fpg1_389, fph_543, \
                         fph_544 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_724[k] = f_1 * dph_354[k]
                   + f_1 * fsh_186[k]
                   + f_7 * fpg0_388[k]
                   - f_8 * fpg1_388[k]
                   + f_4 * pc_y[k] * fph_543[k];

        t_725[k] = f_1 * dph_355[k]
                   + f_1 * fsh_187[k]
                   + f_5 * fpg0_389[k]
                   - f_6 * fpg1_389[k]
                   + f_4 * pc_y[k] * fph_544[k];
    }

#pragma omp simd aligned(t_726, t_727, t_728, pa_y, pc_y, pc_z, dpi0_476, dph_293, dph_356, \
                         dpi1_476, fsh_188, fpg0_389, fpg1_389, \
                         fph_545 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_726[k] = f_1 * dph_356[k]
                   + f_1 * fsh_188[k]
                   + f_4 * pc_y[k] * fph_545[k];

        t_727[k] = f_12 * dph_293[k]
                   + f_2 * fpg0_389[k]
                   - f_3 * fpg1_389[k]
                   + f_4 * pc_z[k] * fph_545[k];

        t_728[k] = pa_y[k] * dpi0_476[k]
                   - f_11 * pc_y[k] * dpi1_476[k];
    }

#pragma omp simd aligned(t_729, t_730, t_731, pa_y, pc_x, pc_y, dpi0_478, dpi1_478, fpg0_391, \
                         fpg0_393, fpg1_391, fpg1_393, fph_547, \
                         fph_549 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_729[k] = f_17 * fpg0_391[k]
                   - f_18 * fpg1_391[k]
                   + f_4 * pc_x[k] * fph_547[k];

        t_730[k] = pa_y[k] * dpi0_478[k]
                   - f_11 * pc_y[k] * dpi1_478[k];

        t_731[k] = f_9 * fpg0_393[k]
                   - f_10 * fpg1_393[k]
                   + f_4 * pc_x[k] * fph_549[k];
    }

#pragma omp simd aligned(t_732, t_733, t_734, pa_y, pc_x, pc_y, dpi0_481, dpi1_481, fpg0_394, \
                         fpg0_396, fpg1_394, fpg1_396, fph_550, \
                         fph_552 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_732[k] = f_9 * fpg0_394[k]
                   - f_10 * fpg1_394[k]
                   + f_4 * pc_x[k] * fph_550[k];

        t_733[k] = pa_y[k] * dpi0_481[k]
                   - f_11 * pc_y[k] * dpi1_481[k];

        t_734[k] = f_7 * fpg0_396[k]
                   - f_8 * fpg1_396[k]
                   + f_4 * pc_x[k] * fph_552[k];
    }

#pragma omp simd aligned(t_735, t_736, t_737, pa_y, pc_x, pc_y, dpi0_485, dpi1_485, fpg0_397, \
                         fpg0_398, fpg1_397, fpg1_398, fph_553, \
                         fph_554 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_735[k] = f_7 * fpg0_397[k]
                   - f_8 * fpg1_397[k]
                   + f_4 * pc_x[k] * fph_553[k];

        t_736[k] = f_7 * fpg0_398[k]
                   - f_8 * fpg1_398[k]
                   + f_4 * pc_x[k] * fph_554[k];

        t_737[k] = pa_y[k] * dpi0_485[k]
                   - f_11 * pc_y[k] * dpi1_485[k];
    }

#pragma omp simd aligned(t_738, t_739, t_740, pc_x, fpg0_400, fpg0_401, fpg0_402, fpg1_400, \
                         fpg1_401, fpg1_402, fph_556, fph_557, \
                         fph_558 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_738[k] = f_5 * fpg0_400[k]
                   - f_6 * fpg1_400[k]
                   + f_4 * pc_x[k] * fph_556[k];

        t_739[k] = f_5 * fpg0_401[k]
                   - f_6 * fpg1_401[k]
                   + f_4 * pc_x[k] * fph_557[k];

        t_740[k] = f_5 * fpg0_402[k]
                   - f_6 * fpg1_402[k]
                   + f_4 * pc_x[k] * fph_558[k];
    }

#pragma omp simd aligned(t_741, t_742, t_743, t_744, t_745, pa_y, pc_x, pc_y, dpi0_490, \
                         dpi1_490, fpg0_403, fpg1_403, fph_559, fph_561, fph_562, \
                         fph_563 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_741[k] = f_5 * fpg0_403[k]
                   - f_6 * fpg1_403[k]
                   + f_4 * pc_x[k] * fph_559[k];

        t_742[k] = pa_y[k] * dpi0_490[k]
                   - f_11 * pc_y[k] * dpi1_490[k];

        t_743[k] = f_4 * pc_x[k] * fph_561[k];

        t_744[k] = f_4 * pc_x[k] * fph_562[k];

        t_745[k] = f_4 * pc_x[k] * fph_563[k];
    }

#pragma omp simd aligned(t_746, t_747, t_748, t_749, pa_y, pc_x, pc_y, dpi0_497, dph_372, \
                         dpi1_497, fph_564, fph_565, fph_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_746[k] = f_4 * pc_x[k] * fph_564[k];

        t_747[k] = f_4 * pc_x[k] * fph_565[k];

        t_748[k] = f_4 * pc_x[k] * fph_566[k];

        t_749[k] = pa_y[k] * dpi0_497[k]
                   + f_14 * dph_372[k]
                   - f_11 * pc_y[k] * dpi1_497[k];
    }

#pragma omp simd aligned(t_750, t_751, t_752, pa_y, pc_y, pc_z, dpi0_499, dpi0_500, dph_309, \
                         dph_374, dph_375, dpi1_499, dpi1_500, fsh_183, \
                         fph_561 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_750[k] = f_12 * dph_309[k]
                   + f_1 * fsh_183[k]
                   + f_4 * pc_z[k] * fph_561[k];

        t_751[k] = pa_y[k] * dpi0_499[k]
                   + f_13 * dph_374[k]
                   - f_11 * pc_y[k] * dpi1_499[k];

        t_752[k] = pa_y[k] * dpi0_500[k]
                   + f_0 * dph_375[k]
                   - f_11 * pc_y[k] * dpi1_500[k];
    }

#pragma omp simd aligned(t_753, t_754, t_755, pa_y, pc_y, dpi0_501, dpi0_503, dph_376, \
                         dph_377, dpi1_501, dpi1_503, fph_566 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_753[k] = pa_y[k] * dpi0_501[k]
                   + f_12 * dph_376[k]
                   - f_11 * pc_y[k] * dpi1_501[k];

        t_754[k] = f_1 * dph_377[k]
                   + f_4 * pc_y[k] * fph_566[k];

        t_755[k] = pa_y[k] * dpi0_503[k]
                   - f_11 * pc_y[k] * dpi1_503[k];
    }

#pragma omp simd aligned(t_756, t_757, t_758, pb_x, pc_x, pc_y, fsi0_252, fsi0_254, fsh_189, \
                         fsh_191, fsi1_252, fsi1_254, fph_567 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_756[k] = pb_x[k] * fsi0_252[k]
                   + f_14 * fsh_189[k]
                   - f_11 * pc_x[k] * fsi1_252[k];

        t_757[k] = f_4 * pc_y[k] * fph_567[k];

        t_758[k] = pb_x[k] * fsi0_254[k]
                   + f_19 * fsh_191[k]
                   - f_11 * pc_x[k] * fsi1_254[k];
    }

#pragma omp simd aligned(t_759, t_760, t_761, pb_x, pc_x, pc_y, fsi0_255, fsi0_257, fsh_192, \
                         fsh_194, fsi1_255, fsi1_257, fph_569 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_759[k] = pb_x[k] * fsi0_255[k]
                   + f_13 * fsh_192[k]
                   - f_11 * pc_x[k] * fsi1_255[k];

        t_760[k] = f_4 * pc_y[k] * fph_569[k];

        t_761[k] = pb_x[k] * fsi0_257[k]
                   + f_13 * fsh_194[k]
                   - f_11 * pc_x[k] * fsi1_257[k];
    }

#pragma omp simd aligned(t_762, t_763, t_764, pb_x, pc_x, pc_y, fsi0_258, fsi0_259, fsh_195, \
                         fsh_196, fsi1_258, fsi1_259, fph_572 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_762[k] = pb_x[k] * fsi0_258[k]
                   + f_0 * fsh_195[k]
                   - f_11 * pc_x[k] * fsi1_258[k];

        t_763[k] = pb_x[k] * fsi0_259[k]
                   + f_0 * fsh_196[k]
                   - f_11 * pc_x[k] * fsi1_259[k];

        t_764[k] = f_4 * pc_y[k] * fph_572[k];
    }

#pragma omp simd aligned(t_765, t_766, t_767, pb_x, pc_x, fsi0_261, fsi0_262, fsi0_263, \
                         fsh_198, fsh_199, fsh_200, fsi1_261, fsi1_262, \
                         fsi1_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_765[k] = pb_x[k] * fsi0_261[k]
                   + f_0 * fsh_198[k]
                   - f_11 * pc_x[k] * fsi1_261[k];

        t_766[k] = pb_x[k] * fsi0_262[k]
                   + f_12 * fsh_199[k]
                   - f_11 * pc_x[k] * fsi1_262[k];

        t_767[k] = pb_x[k] * fsi0_263[k]
                   + f_12 * fsh_200[k]
                   - f_11 * pc_x[k] * fsi1_263[k];
    }

#pragma omp simd aligned(t_768, t_769, t_770, t_771, pb_x, pc_x, pc_y, fsi0_264, fsi0_266, \
                         fsh_201, fsh_203, fsh_204, fsi1_264, fsi1_266, fph_576, \
                         fph_582 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_768[k] = pb_x[k] * fsi0_264[k]
                   + f_12 * fsh_201[k]
                   - f_11 * pc_x[k] * fsi1_264[k];

        t_769[k] = f_4 * pc_y[k] * fph_576[k];

        t_770[k] = pb_x[k] * fsi0_266[k]
                   + f_12 * fsh_203[k]
                   - f_11 * pc_x[k] * fsi1_266[k];

        t_771[k] = f_1 * fsh_204[k]
                   + f_4 * pc_x[k] * fph_582[k];
    }

#pragma omp simd aligned(t_772, t_773, t_774, t_775, t_776, pc_x, fsh_205, fsh_206, fsh_207, \
                         fsh_208, fsh_209, fph_583, fph_584, fph_585, fph_586, \
                         fph_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_772[k] = f_1 * fsh_205[k]
                   + f_4 * pc_x[k] * fph_583[k];

        t_773[k] = f_1 * fsh_206[k]
                   + f_4 * pc_x[k] * fph_584[k];

        t_774[k] = f_1 * fsh_207[k]
                   + f_4 * pc_x[k] * fph_585[k];

        t_775[k] = f_1 * fsh_208[k]
                   + f_4 * pc_x[k] * fph_586[k];

        t_776[k] = f_1 * fsh_209[k]
                   + f_4 * pc_x[k] * fph_587[k];
    }

#pragma omp simd aligned(t_777, t_778, t_779, t_780, pb_x, pc_x, fsi0_273, fsi0_274, fsi0_275, \
                         fsi0_276, fsi1_273, fsi1_274, fsi1_275, \
                         fsi1_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_777[k] = pb_x[k] * fsi0_273[k]
                   - f_11 * pc_x[k] * fsi1_273[k];

        t_778[k] = pb_x[k] * fsi0_274[k]
                   - f_11 * pc_x[k] * fsi1_274[k];

        t_779[k] = pb_x[k] * fsi0_275[k]
                   - f_11 * pc_x[k] * fsi1_275[k];

        t_780[k] = pb_x[k] * fsi0_276[k]
                   - f_11 * pc_x[k] * fsi1_276[k];
    }

#pragma omp simd aligned(t_781, t_782, t_783, t_784, pb_x, pb_y, pc_x, pc_y, fsi0_252, \
                         fsi0_277, fsi0_279, fsi1_252, fsi1_277, fsi1_279, \
                         fph_587 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_781[k] = pb_x[k] * fsi0_277[k]
                   - f_11 * pc_x[k] * fsi1_277[k];

        t_782[k] = f_4 * pc_y[k] * fph_587[k];

        t_783[k] = pb_x[k] * fsi0_279[k]
                   - f_11 * pc_x[k] * fsi1_279[k];

        t_784[k] = pb_y[k] * fsi0_252[k]
                   - f_11 * pc_y[k] * fsi1_252[k];
    }

#pragma omp simd aligned(t_785, t_786, t_787, t_788, pb_y, pc_x, pc_y, fsi0_254, fsh_189, \
                         fsh_191, fsi1_254, fpg0_423, fpg1_423, fph_588, fph_590, \
                         fph_591 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_785[k] = f_1 * fsh_189[k]
                   + f_4 * pc_y[k] * fph_588[k];

        t_786[k] = pb_y[k] * fsi0_254[k]
                   - f_11 * pc_y[k] * fsi1_254[k];

        t_787[k] = f_9 * fpg0_423[k]
                   - f_10 * fpg1_423[k]
                   + f_4 * pc_x[k] * fph_591[k];

        t_788[k] = f_1 * fsh_191[k]
                   + f_4 * pc_y[k] * fph_590[k];
    }

#pragma omp simd aligned(t_789, t_790, t_791, pb_y, pc_x, pc_y, fsi0_257, fsi1_257, fpg0_426, \
                         fpg0_427, fpg1_426, fpg1_427, fph_594, \
                         fph_595 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_789[k] = pb_y[k] * fsi0_257[k]
                   - f_11 * pc_y[k] * fsi1_257[k];

        t_790[k] = f_7 * fpg0_426[k]
                   - f_8 * fpg1_426[k]
                   + f_4 * pc_x[k] * fph_594[k];

        t_791[k] = f_7 * fpg0_427[k]
                   - f_8 * fpg1_427[k]
                   + f_4 * pc_x[k] * fph_595[k];
    }

#pragma omp simd aligned(t_792, t_793, t_794, pb_y, pc_x, pc_y, fsi0_261, fsh_194, fsi1_261, \
                         fpg0_430, fpg1_430, fph_593, fph_598 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_792[k] = f_1 * fsh_194[k]
                   + f_4 * pc_y[k] * fph_593[k];

        t_793[k] = pb_y[k] * fsi0_261[k]
                   - f_11 * pc_y[k] * fsi1_261[k];

        t_794[k] = f_5 * fpg0_430[k]
                   - f_6 * fpg1_430[k]
                   + f_4 * pc_x[k] * fph_598[k];
    }

#pragma omp simd aligned(t_795, t_796, t_797, pc_x, pc_y, fsh_198, fpg0_431, fpg0_432, \
                         fpg1_431, fpg1_432, fph_597, fph_599, \
                         fph_600 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_795[k] = f_5 * fpg0_431[k]
                   - f_6 * fpg1_431[k]
                   + f_4 * pc_x[k] * fph_599[k];

        t_796[k] = f_5 * fpg0_432[k]
                   - f_6 * fpg1_432[k]
                   + f_4 * pc_x[k] * fph_600[k];

        t_797[k] = f_1 * fsh_198[k]
                   + f_4 * pc_y[k] * fph_597[k];
    }

#pragma omp simd aligned(t_798, t_799, t_800, t_801, t_802, t_803, pb_y, pc_x, pc_y, fsi0_266, \
                         fsi1_266, fph_603, fph_604, fph_605, fph_606, \
                         fph_607 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_798[k] = pb_y[k] * fsi0_266[k]
                   - f_11 * pc_y[k] * fsi1_266[k];

        t_799[k] = f_4 * pc_x[k] * fph_603[k];

        t_800[k] = f_4 * pc_x[k] * fph_604[k];

        t_801[k] = f_4 * pc_x[k] * fph_605[k];

        t_802[k] = f_4 * pc_x[k] * fph_606[k];

        t_803[k] = f_4 * pc_x[k] * fph_607[k];
    }

#pragma omp simd aligned(t_804, t_805, t_806, pb_y, pc_x, pc_y, fsi0_273, fsi0_274, fsh_204, \
                         fsh_205, fsi1_273, fsi1_274, fph_608 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_804[k] = f_4 * pc_x[k] * fph_608[k];

        t_805[k] = pb_y[k] * fsi0_273[k]
                   + f_14 * fsh_204[k]
                   - f_11 * pc_y[k] * fsi1_273[k];

        t_806[k] = pb_y[k] * fsi0_274[k]
                   + f_19 * fsh_205[k]
                   - f_11 * pc_y[k] * fsi1_274[k];
    }

#pragma omp simd aligned(t_807, t_808, t_809, pb_y, pc_y, fsi0_275, fsi0_276, fsi0_277, \
                         fsh_206, fsh_207, fsh_208, fsi1_275, fsi1_276, \
                         fsi1_277 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_807[k] = pb_y[k] * fsi0_275[k]
                   + f_13 * fsh_206[k]
                   - f_11 * pc_y[k] * fsi1_275[k];

        t_808[k] = pb_y[k] * fsi0_276[k]
                   + f_0 * fsh_207[k]
                   - f_11 * pc_y[k] * fsi1_276[k];

        t_809[k] = pb_y[k] * fsi0_277[k]
                   + f_12 * fsh_208[k]
                   - f_11 * pc_y[k] * fsi1_277[k];
    }

#pragma omp simd aligned(t_810, t_811, t_812, t_813, pb_y, pc_x, pc_y, fsi0_279, fsh_209, \
                         fsi1_279, fpg0_435, fpg1_435, fph_608, \
                         fph_609 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_810[k] = f_1 * fsh_209[k]
                   + f_4 * pc_y[k] * fph_608[k];

        t_811[k] = pb_y[k] * fsi0_279[k]
                   - f_11 * pc_y[k] * fsi1_279[k];

        t_812[k] = f_2 * fpg0_435[k]
                   - f_3 * fpg1_435[k]
                   + f_4 * pc_x[k] * fph_609[k];

        t_813[k] = f_4 * pc_y[k] * fph_609[k];
    }

#pragma omp simd aligned(t_814, t_815, t_816, t_817, pc_x, pc_y, fpg0_437, fpg0_438, fpg0_440, \
                         fpg1_437, fpg1_438, fpg1_440, fph_611, fph_612, \
                         fph_614 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_814[k] = f_17 * fpg0_437[k]
                   - f_18 * fpg1_437[k]
                   + f_4 * pc_x[k] * fph_611[k];

        t_815[k] = f_9 * fpg0_438[k]
                   - f_10 * fpg1_438[k]
                   + f_4 * pc_x[k] * fph_612[k];

        t_816[k] = f_4 * pc_y[k] * fph_611[k];

        t_817[k] = f_9 * fpg0_440[k]
                   - f_10 * fpg1_440[k]
                   + f_4 * pc_x[k] * fph_614[k];
    }

#pragma omp simd aligned(t_818, t_819, t_820, t_821, pc_x, pc_y, fpg0_441, fpg0_442, fpg0_444, \
                         fpg1_441, fpg1_442, fpg1_444, fph_614, fph_615, fph_616, \
                         fph_618 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_818[k] = f_7 * fpg0_441[k]
                   - f_8 * fpg1_441[k]
                   + f_4 * pc_x[k] * fph_615[k];

        t_819[k] = f_7 * fpg0_442[k]
                   - f_8 * fpg1_442[k]
                   + f_4 * pc_x[k] * fph_616[k];

        t_820[k] = f_4 * pc_y[k] * fph_614[k];

        t_821[k] = f_7 * fpg0_444[k]
                   - f_8 * fpg1_444[k]
                   + f_4 * pc_x[k] * fph_618[k];
    }

#pragma omp simd aligned(t_822, t_823, t_824, t_825, pc_x, pc_y, fpg0_445, fpg0_446, fpg0_447, \
                         fpg1_445, fpg1_446, fpg1_447, fph_618, fph_619, fph_620, \
                         fph_621 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_822[k] = f_5 * fpg0_445[k]
                   - f_6 * fpg1_445[k]
                   + f_4 * pc_x[k] * fph_619[k];

        t_823[k] = f_5 * fpg0_446[k]
                   - f_6 * fpg1_446[k]
                   + f_4 * pc_x[k] * fph_620[k];

        t_824[k] = f_5 * fpg0_447[k]
                   - f_6 * fpg1_447[k]
                   + f_4 * pc_x[k] * fph_621[k];

        t_825[k] = f_4 * pc_y[k] * fph_618[k];
    }

#pragma omp simd aligned(t_826, t_827, t_828, t_829, t_830, t_831, pc_x, fpg0_449, fpg1_449, \
                         fph_623, fph_624, fph_625, fph_626, fph_627, \
                         fph_628 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_826[k] = f_5 * fpg0_449[k]
                   - f_6 * fpg1_449[k]
                   + f_4 * pc_x[k] * fph_623[k];

        t_827[k] = f_4 * pc_x[k] * fph_624[k];

        t_828[k] = f_4 * pc_x[k] * fph_625[k];

        t_829[k] = f_4 * pc_x[k] * fph_626[k];

        t_830[k] = f_4 * pc_x[k] * fph_627[k];

        t_831[k] = f_4 * pc_x[k] * fph_628[k];
    }

#pragma omp simd aligned(t_832, t_833, t_834, t_835, pc_x, pc_y, fpg0_445, fpg0_446, fpg0_447, \
                         fpg1_445, fpg1_446, fpg1_447, fph_624, fph_625, fph_626, \
                         fph_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_832[k] = f_4 * pc_x[k] * fph_629[k];

        t_833[k] = f_2 * fpg0_445[k]
                   - f_3 * fpg1_445[k]
                   + f_4 * pc_y[k] * fph_624[k];

        t_834[k] = f_17 * fpg0_446[k]
                   - f_18 * fpg1_446[k]
                   + f_4 * pc_y[k] * fph_625[k];

        t_835[k] = f_9 * fpg0_447[k]
                   - f_10 * fpg1_447[k]
                   + f_4 * pc_y[k] * fph_626[k];
    }

#pragma omp simd aligned(t_836, t_837, t_838, t_839, pc_y, pc_z, dph_377, fsh_209, fpg0_448, \
                         fpg0_449, fpg1_448, fpg1_449, fph_627, fph_628, \
                         fph_629 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_836[k] = f_7 * fpg0_448[k]
                   - f_8 * fpg1_448[k]
                   + f_4 * pc_y[k] * fph_627[k];

        t_837[k] = f_5 * fpg0_449[k]
                   - f_6 * fpg1_449[k]
                   + f_4 * pc_y[k] * fph_628[k];

        t_838[k] = f_4 * pc_y[k] * fph_629[k];

        t_839[k] = f_0 * dph_377[k]
                   + f_1 * fsh_209[k]
                   + f_2 * fpg0_449[k]
                   - f_3 * fpg1_449[k]
                   + f_4 * pc_z[k] * fph_629[k];
    }
}

auto
compute_prim_fpi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppi0,
                                                   const size_t ppi1, const size_t dpi0,
                                                   const size_t dph, const size_t dpi1,
                                                   const size_t fsi0, const size_t fsh,
                                                   const size_t fsi1, const size_t fpg0,
                                                   const size_t fpg1, const size_t fph,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_fpi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, dpi0,
                                                              dph, dpi1, fsi0, fsh, fsi1, fpg0,
                                                              fpg1, fph, ncols, gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppi0,
                                                              ppi1, dpi0, dph, dpi1, fsi0, fsh,
                                                              fsi1, fpg0, fpg1, fph, ncols,
                                                              gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, ppi0,
                                                              ppi1, dpi0, dph, dpi1, fsi0, fsh,
                                                              fsi1, fpg0, fpg1, fph, ncols,
                                                              gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, dpi0,
                                                              dph, dpi1, fsi0, fsh, fsi1, fpg0,
                                                              fpg1, fph, ncols, gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece4(buffer, target, pa, pb, pc, dpi0,
                                                              dph, dpi1, fsi0, fsh, fsi1, fpg0,
                                                              fpg1, fph, ncols, gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece5(buffer, target, pa, pb, pc, ppi0,
                                                              ppi1, dpi0, dph, dpi1, fsi0, fsh,
                                                              fsi1, fpg0, fpg1, fph, ncols,
                                                              gamma, p, q);

    compute_prim_fpi_three_center_electron_repulsion_0_piece6(buffer, target, pa, pb, pc, dpi0,
                                                              dph, dpi1, fsi0, fsh, fsi1, fpg0,
                                                              fpg1, fph, ncols, gamma, p, q);
}

}  // namespace simdt3ceri
