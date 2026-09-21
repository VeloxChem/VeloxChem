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


#include "SimdThreeCenterElectronRepulsionVrrRecDPI.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_dpi_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t pph,
                                                          const size_t ppi1, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t dph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_12 = 1.5 / q;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_0 = buffer.data(ppi0 + 0);
    const auto *ppi0_5 = buffer.data(ppi0 + 5);
    const auto *ppi0_9 = buffer.data(ppi0 + 9);
    const auto *ppi0_14 = buffer.data(ppi0 + 14);
    const auto *ppi0_20 = buffer.data(ppi0 + 20);
    const auto *ppi0_115 = buffer.data(ppi0 + 115);
    const auto *ppi0_118 = buffer.data(ppi0 + 118);
    const auto *ppi0_122 = buffer.data(ppi0 + 122);

    const auto *pph_0 = buffer.data(pph + 0);
    const auto *pph_5 = buffer.data(pph + 5);
    const auto *pph_9 = buffer.data(pph + 9);
    const auto *pph_15 = buffer.data(pph + 15);
    const auto *pph_17 = buffer.data(pph + 17);
    const auto *pph_18 = buffer.data(pph + 18);
    const auto *pph_20 = buffer.data(pph + 20);
    const auto *pph_21 = buffer.data(pph + 21);
    const auto *pph_26 = buffer.data(pph + 26);
    const auto *pph_36 = buffer.data(pph + 36);
    const auto *pph_38 = buffer.data(pph + 38);
    const auto *pph_39 = buffer.data(pph + 39);
    const auto *pph_41 = buffer.data(pph + 41);
    const auto *pph_57 = buffer.data(pph + 57);
    const auto *pph_59 = buffer.data(pph + 59);
    const auto *pph_60 = buffer.data(pph + 60);
    const auto *pph_62 = buffer.data(pph + 62);
    const auto *pph_66 = buffer.data(pph + 66);
    const auto *pph_69 = buffer.data(pph + 69);
    const auto *pph_73 = buffer.data(pph + 73);
    const auto *pph_78 = buffer.data(pph + 78);
    const auto *pph_80 = buffer.data(pph + 80);
    const auto *pph_81 = buffer.data(pph + 81);
    const auto *pph_82 = buffer.data(pph + 82);
    const auto *pph_84 = buffer.data(pph + 84);
    const auto *pph_87 = buffer.data(pph + 87);
    const auto *pph_90 = buffer.data(pph + 90);
    const auto *pph_94 = buffer.data(pph + 94);

    const auto *ppi1_0 = buffer.data(ppi1 + 0);
    const auto *ppi1_5 = buffer.data(ppi1 + 5);
    const auto *ppi1_9 = buffer.data(ppi1 + 9);
    const auto *ppi1_14 = buffer.data(ppi1 + 14);
    const auto *ppi1_20 = buffer.data(ppi1 + 20);
    const auto *ppi1_115 = buffer.data(ppi1 + 115);
    const auto *ppi1_118 = buffer.data(ppi1 + 118);
    const auto *ppi1_122 = buffer.data(ppi1 + 122);

    const auto *dsi0_0 = buffer.data(dsi0 + 0);
    const auto *dsi0_3 = buffer.data(dsi0 + 3);
    const auto *dsi0_5 = buffer.data(dsi0 + 5);
    const auto *dsi0_6 = buffer.data(dsi0 + 6);
    const auto *dsi0_9 = buffer.data(dsi0 + 9);
    const auto *dsi0_10 = buffer.data(dsi0 + 10);
    const auto *dsi0_12 = buffer.data(dsi0 + 12);
    const auto *dsi0_14 = buffer.data(dsi0 + 14);
    const auto *dsi0_21 = buffer.data(dsi0 + 21);
    const auto *dsi0_23 = buffer.data(dsi0 + 23);
    const auto *dsi0_24 = buffer.data(dsi0 + 24);
    const auto *dsi0_25 = buffer.data(dsi0 + 25);
    const auto *dsi0_27 = buffer.data(dsi0 + 27);

    const auto *dsh_0 = buffer.data(dsh + 0);
    const auto *dsh_1 = buffer.data(dsh + 1);
    const auto *dsh_2 = buffer.data(dsh + 2);
    const auto *dsh_3 = buffer.data(dsh + 3);
    const auto *dsh_5 = buffer.data(dsh + 5);
    const auto *dsh_6 = buffer.data(dsh + 6);
    const auto *dsh_8 = buffer.data(dsh + 8);
    const auto *dsh_9 = buffer.data(dsh + 9);
    const auto *dsh_10 = buffer.data(dsh + 10);
    const auto *dsh_14 = buffer.data(dsh + 14);
    const auto *dsh_15 = buffer.data(dsh + 15);
    const auto *dsh_17 = buffer.data(dsh + 17);
    const auto *dsh_18 = buffer.data(dsh + 18);
    const auto *dsh_19 = buffer.data(dsh + 19);
    const auto *dsh_20 = buffer.data(dsh + 20);
    const auto *dsh_21 = buffer.data(dsh + 21);
    const auto *dsh_24 = buffer.data(dsh + 24);
    const auto *dsh_26 = buffer.data(dsh + 26);
    const auto *dsh_27 = buffer.data(dsh + 27);
    const auto *dsh_31 = buffer.data(dsh + 31);
    const auto *dsh_36 = buffer.data(dsh + 36);
    const auto *dsh_38 = buffer.data(dsh + 38);
    const auto *dsh_39 = buffer.data(dsh + 39);
    const auto *dsh_40 = buffer.data(dsh + 40);

    const auto *dsi1_0 = buffer.data(dsi1 + 0);
    const auto *dsi1_3 = buffer.data(dsi1 + 3);
    const auto *dsi1_5 = buffer.data(dsi1 + 5);
    const auto *dsi1_6 = buffer.data(dsi1 + 6);
    const auto *dsi1_9 = buffer.data(dsi1 + 9);
    const auto *dsi1_10 = buffer.data(dsi1 + 10);
    const auto *dsi1_12 = buffer.data(dsi1 + 12);
    const auto *dsi1_14 = buffer.data(dsi1 + 14);
    const auto *dsi1_21 = buffer.data(dsi1 + 21);
    const auto *dsi1_23 = buffer.data(dsi1 + 23);
    const auto *dsi1_24 = buffer.data(dsi1 + 24);
    const auto *dsi1_25 = buffer.data(dsi1 + 25);
    const auto *dsi1_27 = buffer.data(dsi1 + 27);

    const auto *dpg0_0 = buffer.data(dpg0 + 0);
    const auto *dpg0_1 = buffer.data(dpg0 + 1);
    const auto *dpg0_2 = buffer.data(dpg0 + 2);
    const auto *dpg0_3 = buffer.data(dpg0 + 3);
    const auto *dpg0_5 = buffer.data(dpg0 + 5);
    const auto *dpg0_10 = buffer.data(dpg0 + 10);
    const auto *dpg0_12 = buffer.data(dpg0 + 12);
    const auto *dpg0_13 = buffer.data(dpg0 + 13);
    const auto *dpg0_14 = buffer.data(dpg0 + 14);
    const auto *dpg0_35 = buffer.data(dpg0 + 35);
    const auto *dpg0_42 = buffer.data(dpg0 + 42);
    const auto *dpg0_43 = buffer.data(dpg0 + 43);
    const auto *dpg0_44 = buffer.data(dpg0 + 44);
    const auto *dpg0_48 = buffer.data(dpg0 + 48);
    const auto *dpg0_51 = buffer.data(dpg0 + 51);
    const auto *dpg0_55 = buffer.data(dpg0 + 55);
    const auto *dpg0_56 = buffer.data(dpg0 + 56);
    const auto *dpg0_57 = buffer.data(dpg0 + 57);
    const auto *dpg0_59 = buffer.data(dpg0 + 59);
    const auto *dpg0_60 = buffer.data(dpg0 + 60);
    const auto *dpg0_62 = buffer.data(dpg0 + 62);

    const auto *dpg1_0 = buffer.data(dpg1 + 0);
    const auto *dpg1_1 = buffer.data(dpg1 + 1);
    const auto *dpg1_2 = buffer.data(dpg1 + 2);
    const auto *dpg1_3 = buffer.data(dpg1 + 3);
    const auto *dpg1_5 = buffer.data(dpg1 + 5);
    const auto *dpg1_10 = buffer.data(dpg1 + 10);
    const auto *dpg1_12 = buffer.data(dpg1 + 12);
    const auto *dpg1_13 = buffer.data(dpg1 + 13);
    const auto *dpg1_14 = buffer.data(dpg1 + 14);
    const auto *dpg1_35 = buffer.data(dpg1 + 35);
    const auto *dpg1_42 = buffer.data(dpg1 + 42);
    const auto *dpg1_43 = buffer.data(dpg1 + 43);
    const auto *dpg1_44 = buffer.data(dpg1 + 44);
    const auto *dpg1_48 = buffer.data(dpg1 + 48);
    const auto *dpg1_51 = buffer.data(dpg1 + 51);
    const auto *dpg1_55 = buffer.data(dpg1 + 55);
    const auto *dpg1_56 = buffer.data(dpg1 + 56);
    const auto *dpg1_57 = buffer.data(dpg1 + 57);
    const auto *dpg1_59 = buffer.data(dpg1 + 59);
    const auto *dpg1_60 = buffer.data(dpg1 + 60);
    const auto *dpg1_62 = buffer.data(dpg1 + 62);

    const auto *dph_0 = buffer.data(dph + 0);
    const auto *dph_1 = buffer.data(dph + 1);
    const auto *dph_2 = buffer.data(dph + 2);
    const auto *dph_3 = buffer.data(dph + 3);
    const auto *dph_5 = buffer.data(dph + 5);
    const auto *dph_6 = buffer.data(dph + 6);
    const auto *dph_8 = buffer.data(dph + 8);
    const auto *dph_9 = buffer.data(dph + 9);
    const auto *dph_10 = buffer.data(dph + 10);
    const auto *dph_14 = buffer.data(dph + 14);
    const auto *dph_15 = buffer.data(dph + 15);
    const auto *dph_17 = buffer.data(dph + 17);
    const auto *dph_18 = buffer.data(dph + 18);
    const auto *dph_19 = buffer.data(dph + 19);
    const auto *dph_20 = buffer.data(dph + 20);
    const auto *dph_21 = buffer.data(dph + 21);
    const auto *dph_23 = buffer.data(dph + 23);
    const auto *dph_24 = buffer.data(dph + 24);
    const auto *dph_26 = buffer.data(dph + 26);
    const auto *dph_27 = buffer.data(dph + 27);
    const auto *dph_30 = buffer.data(dph + 30);
    const auto *dph_31 = buffer.data(dph + 31);
    const auto *dph_35 = buffer.data(dph + 35);
    const auto *dph_36 = buffer.data(dph + 36);
    const auto *dph_38 = buffer.data(dph + 38);
    const auto *dph_39 = buffer.data(dph + 39);
    const auto *dph_41 = buffer.data(dph + 41);
    const auto *dph_42 = buffer.data(dph + 42);
    const auto *dph_44 = buffer.data(dph + 44);
    const auto *dph_45 = buffer.data(dph + 45);
    const auto *dph_47 = buffer.data(dph + 47);
    const auto *dph_48 = buffer.data(dph + 48);
    const auto *dph_50 = buffer.data(dph + 50);
    const auto *dph_51 = buffer.data(dph + 51);
    const auto *dph_52 = buffer.data(dph + 52);
    const auto *dph_56 = buffer.data(dph + 56);
    const auto *dph_57 = buffer.data(dph + 57);
    const auto *dph_59 = buffer.data(dph + 59);
    const auto *dph_60 = buffer.data(dph + 60);
    const auto *dph_61 = buffer.data(dph + 61);
    const auto *dph_62 = buffer.data(dph + 62);
    const auto *dph_63 = buffer.data(dph + 63);
    const auto *dph_64 = buffer.data(dph + 64);
    const auto *dph_66 = buffer.data(dph + 66);
    const auto *dph_68 = buffer.data(dph + 68);
    const auto *dph_69 = buffer.data(dph + 69);
    const auto *dph_70 = buffer.data(dph + 70);
    const auto *dph_72 = buffer.data(dph + 72);
    const auto *dph_73 = buffer.data(dph + 73);
    const auto *dph_78 = buffer.data(dph + 78);
    const auto *dph_79 = buffer.data(dph + 79);
    const auto *dph_80 = buffer.data(dph + 80);
    const auto *dph_81 = buffer.data(dph + 81);
    const auto *dph_82 = buffer.data(dph + 82);
    const auto *dph_83 = buffer.data(dph + 83);
    const auto *dph_84 = buffer.data(dph + 84);
    const auto *dph_85 = buffer.data(dph + 85);
    const auto *dph_86 = buffer.data(dph + 86);
    const auto *dph_87 = buffer.data(dph + 87);
    const auto *dph_89 = buffer.data(dph + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, pph_0, dsh_0, dpg0_0, \
                         dpg1_0, dph_0, dph_1, dph_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * pph_0[k]
                 + f_1 * dsh_0[k]
                 + f_2 * dpg0_0[k]
                 - f_3 * dpg1_0[k]
                 + f_4 * pc_x[k] * dph_0[k];

        t_1[k] = f_4 * pc_y[k] * dph_0[k];

        t_2[k] = f_4 * pc_z[k] * dph_0[k];

        t_3[k] = f_5 * dpg0_0[k]
                 - f_6 * dpg1_0[k]
                 + f_4 * pc_y[k] * dph_1[k];

        t_4[k] = f_4 * pc_y[k] * dph_2[k];

        t_5[k] = f_5 * dpg0_0[k]
                 - f_6 * dpg1_0[k]
                 + f_4 * pc_z[k] * dph_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pc_y, pc_z, dpg0_1, dpg0_2, dpg0_3, dpg1_1, \
                         dpg1_2, dpg1_3, dph_3, dph_5, dph_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_7 * dpg0_1[k]
                 - f_8 * dpg1_1[k]
                 + f_4 * pc_y[k] * dph_3[k];

        t_7[k] = f_4 * pc_z[k] * dph_3[k];

        t_8[k] = f_4 * pc_y[k] * dph_5[k];

        t_9[k] = f_7 * dpg0_2[k]
                 - f_8 * dpg1_2[k]
                 + f_4 * pc_z[k] * dph_5[k];

        t_10[k] = f_9 * dpg0_3[k]
                  - f_10 * dpg1_3[k]
                  + f_4 * pc_y[k] * dph_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pc_x, pc_y, pc_z, pph_15, dsh_15, \
                         dpg0_5, dpg1_5, dph_6, dph_8, dph_9, dph_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_4 * pc_z[k] * dph_6[k];

        t_12[k] = f_5 * dpg0_5[k]
                  - f_6 * dpg1_5[k]
                  + f_4 * pc_y[k] * dph_8[k];

        t_13[k] = f_4 * pc_y[k] * dph_9[k];

        t_14[k] = f_9 * dpg0_5[k]
                  - f_10 * dpg1_5[k]
                  + f_4 * pc_z[k] * dph_9[k];

        t_15[k] = f_0 * pph_15[k]
                  + f_1 * dsh_15[k]
                  + f_4 * pc_x[k] * dph_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pc_x, pc_y, pc_z, pph_17, pph_18, dsh_17, \
                         dsh_18, dph_10, dph_14, dph_17, dph_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_4 * pc_z[k] * dph_10[k];

        t_17[k] = f_0 * pph_17[k]
                  + f_1 * dsh_17[k]
                  + f_4 * pc_x[k] * dph_17[k];

        t_18[k] = f_0 * pph_18[k]
                  + f_1 * dsh_18[k]
                  + f_4 * pc_x[k] * dph_18[k];

        t_19[k] = f_4 * pc_y[k] * dph_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, t_23, pc_x, pc_y, pc_z, pph_20, dsh_20, dpg0_10, \
                         dpg0_12, dpg1_10, dpg1_12, dph_15, dph_17, \
                         dph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = f_0 * pph_20[k]
                  + f_1 * dsh_20[k]
                  + f_4 * pc_x[k] * dph_20[k];

        t_21[k] = f_2 * dpg0_10[k]
                  - f_3 * dpg1_10[k]
                  + f_4 * pc_y[k] * dph_15[k];

        t_22[k] = f_4 * pc_z[k] * dph_15[k];

        t_23[k] = f_9 * dpg0_12[k]
                  - f_10 * dpg1_12[k]
                  + f_4 * pc_y[k] * dph_17[k];
    }

#pragma omp simd aligned(t_24, t_25, t_26, t_27, pc_y, pc_z, dpg0_13, dpg0_14, dpg1_13, \
                         dpg1_14, dph_18, dph_19, dph_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_24[k] = f_7 * dpg0_13[k]
                  - f_8 * dpg1_13[k]
                  + f_4 * pc_y[k] * dph_18[k];

        t_25[k] = f_5 * dpg0_14[k]
                  - f_6 * dpg1_14[k]
                  + f_4 * pc_y[k] * dph_19[k];

        t_26[k] = f_4 * pc_y[k] * dph_20[k];

        t_27[k] = f_2 * dpg0_14[k]
                  - f_3 * dpg1_14[k]
                  + f_4 * pc_z[k] * dph_20[k];
    }

#pragma omp simd aligned(t_28, t_29, t_30, t_31, pb_y, pc_y, pc_z, dsi0_0, dsi0_3, dsh_0, \
                         dsh_1, dsi1_0, dsi1_3, dph_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_28[k] = pb_y[k] * dsi0_0[k]
                  - f_11 * pc_y[k] * dsi1_0[k];

        t_29[k] = f_1 * dsh_0[k]
                  + f_4 * pc_y[k] * dph_21[k];

        t_30[k] = f_4 * pc_z[k] * dph_21[k];

        t_31[k] = pb_y[k] * dsi0_3[k]
                  + f_0 * dsh_1[k]
                  - f_11 * pc_y[k] * dsi1_3[k];
    }

#pragma omp simd aligned(t_32, t_33, t_34, t_35, pb_y, pc_y, pc_z, dsi0_5, dsi0_6, dsh_2, \
                         dsh_3, dsi1_5, dsi1_6, dph_23, dph_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_32[k] = f_1 * dsh_2[k]
                  + f_4 * pc_y[k] * dph_23[k];

        t_33[k] = pb_y[k] * dsi0_5[k]
                  - f_11 * pc_y[k] * dsi1_5[k];

        t_34[k] = pb_y[k] * dsi0_6[k]
                  + f_12 * dsh_3[k]
                  - f_11 * pc_y[k] * dsi1_6[k];

        t_35[k] = f_4 * pc_z[k] * dph_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pb_y, pc_y, pc_z, dsi0_9, dsi0_10, dsh_5, \
                         dsh_6, dsi1_9, dsi1_10, dph_26, dph_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_1 * dsh_5[k]
                  + f_4 * pc_y[k] * dph_26[k];

        t_37[k] = pb_y[k] * dsi0_9[k]
                  - f_11 * pc_y[k] * dsi1_9[k];

        t_38[k] = pb_y[k] * dsi0_10[k]
                  + f_13 * dsh_6[k]
                  - f_11 * pc_y[k] * dsi1_10[k];

        t_39[k] = f_4 * pc_z[k] * dph_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pb_y, pc_x, pc_y, pph_36, dsi0_12, dsi0_14, \
                         dsh_8, dsh_9, dsi1_12, dsi1_14, dph_30, \
                         dph_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pb_y[k] * dsi0_12[k]
                  + f_0 * dsh_8[k]
                  - f_11 * pc_y[k] * dsi1_12[k];

        t_41[k] = f_1 * dsh_9[k]
                  + f_4 * pc_y[k] * dph_30[k];

        t_42[k] = pb_y[k] * dsi0_14[k]
                  - f_11 * pc_y[k] * dsi1_14[k];

        t_43[k] = f_0 * pph_36[k]
                  + f_4 * pc_x[k] * dph_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, pc_x, pc_y, pc_z, pph_38, pph_39, dsh_14, \
                         dph_31, dph_35, dph_38, dph_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = f_4 * pc_z[k] * dph_31[k];

        t_45[k] = f_0 * pph_38[k]
                  + f_4 * pc_x[k] * dph_38[k];

        t_46[k] = f_0 * pph_39[k]
                  + f_4 * pc_x[k] * dph_39[k];

        t_47[k] = f_1 * dsh_14[k]
                  + f_4 * pc_y[k] * dph_35[k];
    }

#pragma omp simd aligned(t_48, t_49, t_50, pb_y, pc_x, pc_y, pc_z, pph_41, dsi0_21, dsh_15, \
                         dsi1_21, dph_36, dph_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_48[k] = f_0 * pph_41[k]
                  + f_4 * pc_x[k] * dph_41[k];

        t_49[k] = pb_y[k] * dsi0_21[k]
                  + f_14 * dsh_15[k]
                  - f_11 * pc_y[k] * dsi1_21[k];

        t_50[k] = f_4 * pc_z[k] * dph_36[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, pb_y, pc_y, dsi0_23, dsi0_24, dsi0_25, dsh_17, \
                         dsh_18, dsh_19, dsi1_23, dsi1_24, dsi1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = pb_y[k] * dsi0_23[k]
                  + f_13 * dsh_17[k]
                  - f_11 * pc_y[k] * dsi1_23[k];

        t_52[k] = pb_y[k] * dsi0_24[k]
                  + f_12 * dsh_18[k]
                  - f_11 * pc_y[k] * dsi1_24[k];

        t_53[k] = pb_y[k] * dsi0_25[k]
                  + f_0 * dsh_19[k]
                  - f_11 * pc_y[k] * dsi1_25[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pb_y, pb_z, pc_y, pc_z, dsi0_0, dsi0_27, \
                         dsh_20, dsi1_0, dsi1_27, dph_41, dph_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_1 * dsh_20[k]
                  + f_4 * pc_y[k] * dph_41[k];

        t_55[k] = pb_y[k] * dsi0_27[k]
                  - f_11 * pc_y[k] * dsi1_27[k];

        t_56[k] = pb_z[k] * dsi0_0[k]
                  - f_11 * pc_z[k] * dsi1_0[k];

        t_57[k] = f_4 * pc_y[k] * dph_42[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, t_61, pb_z, pc_y, pc_z, dsi0_3, dsi0_5, dsh_0, \
                         dsh_2, dsi1_3, dsi1_5, dph_42, dph_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_1 * dsh_0[k]
                  + f_4 * pc_z[k] * dph_42[k];

        t_59[k] = pb_z[k] * dsi0_3[k]
                  - f_11 * pc_z[k] * dsi1_3[k];

        t_60[k] = f_4 * pc_y[k] * dph_44[k];

        t_61[k] = pb_z[k] * dsi0_5[k]
                  + f_0 * dsh_2[k]
                  - f_11 * pc_z[k] * dsi1_5[k];
    }

#pragma omp simd aligned(t_62, t_63, t_64, t_65, pb_z, pc_y, pc_z, dsi0_6, dsi0_9, dsh_3, \
                         dsh_5, dsi1_6, dsi1_9, dph_45, dph_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = pb_z[k] * dsi0_6[k]
                  - f_11 * pc_z[k] * dsi1_6[k];

        t_63[k] = f_1 * dsh_3[k]
                  + f_4 * pc_z[k] * dph_45[k];

        t_64[k] = f_4 * pc_y[k] * dph_47[k];

        t_65[k] = pb_z[k] * dsi0_9[k]
                  + f_12 * dsh_5[k]
                  - f_11 * pc_z[k] * dsi1_9[k];
    }

#pragma omp simd aligned(t_66, t_67, t_68, t_69, pb_z, pc_y, pc_z, dsi0_10, dsh_6, dsi1_10, \
                         dpg0_35, dpg1_35, dph_48, dph_50, dph_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_66[k] = pb_z[k] * dsi0_10[k]
                  - f_11 * pc_z[k] * dsi1_10[k];

        t_67[k] = f_1 * dsh_6[k]
                  + f_4 * pc_z[k] * dph_48[k];

        t_68[k] = f_5 * dpg0_35[k]
                  - f_6 * dpg1_35[k]
                  + f_4 * pc_y[k] * dph_50[k];

        t_69[k] = f_4 * pc_y[k] * dph_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, t_73, pb_z, pc_x, pc_z, pph_57, pph_59, dsi0_14, \
                         dsh_9, dsh_10, dsi1_14, dph_52, dph_57, \
                         dph_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = pb_z[k] * dsi0_14[k]
                  + f_13 * dsh_9[k]
                  - f_11 * pc_z[k] * dsi1_14[k];

        t_71[k] = f_0 * pph_57[k]
                  + f_4 * pc_x[k] * dph_57[k];

        t_72[k] = f_1 * dsh_10[k]
                  + f_4 * pc_z[k] * dph_52[k];

        t_73[k] = f_0 * pph_59[k]
                  + f_4 * pc_x[k] * dph_59[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, t_77, pb_z, pc_x, pc_y, pc_z, pph_60, pph_62, \
                         dsi0_21, dsi1_21, dph_56, dph_60, dph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_0 * pph_60[k]
                  + f_4 * pc_x[k] * dph_60[k];

        t_75[k] = f_4 * pc_y[k] * dph_56[k];

        t_76[k] = f_0 * pph_62[k]
                  + f_4 * pc_x[k] * dph_62[k];

        t_77[k] = pb_z[k] * dsi0_21[k]
                  - f_11 * pc_z[k] * dsi1_21[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, pc_y, pc_z, dsh_15, dpg0_42, dpg0_43, dpg1_42, \
                         dpg1_43, dph_57, dph_59, dph_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = f_1 * dsh_15[k]
                  + f_4 * pc_z[k] * dph_57[k];

        t_79[k] = f_9 * dpg0_42[k]
                  - f_10 * dpg1_42[k]
                  + f_4 * pc_y[k] * dph_59[k];

        t_80[k] = f_7 * dpg0_43[k]
                  - f_8 * dpg1_43[k]
                  + f_4 * pc_y[k] * dph_60[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, pb_z, pc_y, pc_z, dsi0_27, dsh_20, dsi1_27, \
                         dpg0_44, dpg1_44, dph_61, dph_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_5 * dpg0_44[k]
                  - f_6 * dpg1_44[k]
                  + f_4 * pc_y[k] * dph_61[k];

        t_82[k] = f_4 * pc_y[k] * dph_62[k];

        t_83[k] = pb_z[k] * dsi0_27[k]
                  + f_14 * dsh_20[k]
                  - f_11 * pc_z[k] * dsi1_27[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, pa_y, pc_y, pc_z, ppi0_0, pph_0, ppi1_0, \
                         dph_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = pa_y[k] * ppi0_0[k]
                  - f_11 * pc_y[k] * ppi1_0[k];

        t_85[k] = f_1 * pph_0[k]
                  + f_4 * pc_y[k] * dph_63[k];

        t_86[k] = f_4 * pc_z[k] * dph_63[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, pa_y, pc_x, pc_y, pc_z, ppi0_5, pph_66, ppi1_5, \
                         dsh_24, dpg0_48, dpg1_48, dph_64, dph_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_1 * pph_66[k]
                  + f_1 * dsh_24[k]
                  + f_9 * dpg0_48[k]
                  - f_10 * dpg1_48[k]
                  + f_4 * pc_x[k] * dph_66[k];

        t_88[k] = f_4 * pc_z[k] * dph_64[k];

        t_89[k] = pa_y[k] * ppi0_5[k]
                  - f_11 * pc_y[k] * ppi1_5[k];
    }

#pragma omp simd aligned(t_90, t_91, t_92, pc_x, pc_y, pc_z, pph_5, pph_69, dsh_27, dpg0_51, \
                         dpg1_51, dph_66, dph_68, dph_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_90[k] = f_1 * pph_69[k]
                  + f_1 * dsh_27[k]
                  + f_7 * dpg0_51[k]
                  - f_8 * dpg1_51[k]
                  + f_4 * pc_x[k] * dph_69[k];

        t_91[k] = f_4 * pc_z[k] * dph_66[k];

        t_92[k] = f_1 * pph_5[k]
                  + f_4 * pc_y[k] * dph_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, pa_y, pc_x, pc_y, pc_z, ppi0_9, pph_73, ppi1_9, \
                         dsh_31, dpg0_55, dpg1_55, dph_69, dph_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = pa_y[k] * ppi0_9[k]
                  - f_11 * pc_y[k] * ppi1_9[k];

        t_94[k] = f_1 * pph_73[k]
                  + f_1 * dsh_31[k]
                  + f_5 * dpg0_55[k]
                  - f_6 * dpg1_55[k]
                  + f_4 * pc_x[k] * dph_73[k];

        t_95[k] = f_4 * pc_z[k] * dph_69[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pa_y, pc_y, pc_z, ppi0_14, pph_9, ppi1_14, dpg0_48, \
                         dpg1_48, dph_70, dph_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * dpg0_48[k]
                  - f_6 * dpg1_48[k]
                  + f_4 * pc_z[k] * dph_70[k];

        t_97[k] = f_1 * pph_9[k]
                  + f_4 * pc_y[k] * dph_72[k];

        t_98[k] = pa_y[k] * ppi0_14[k]
                  - f_11 * pc_y[k] * ppi1_14[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pc_x, pc_z, pph_78, pph_80, pph_81, \
                         dsh_36, dsh_38, dsh_39, dph_73, dph_78, dph_80, \
                         dph_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_1 * pph_78[k]
                  + f_1 * dsh_36[k]
                  + f_4 * pc_x[k] * dph_78[k];

        t_100[k] = f_4 * pc_z[k] * dph_73[k];

        t_101[k] = f_1 * pph_80[k]
                   + f_1 * dsh_38[k]
                   + f_4 * pc_x[k] * dph_80[k];

        t_102[k] = f_1 * pph_81[k]
                   + f_1 * dsh_39[k]
                   + f_4 * pc_x[k] * dph_81[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_y, pc_x, pc_y, ppi0_20, pph_15, pph_82, \
                         ppi1_20, dsh_40, dpg0_55, dpg1_55, dph_78, \
                         dph_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * pph_82[k]
                   + f_1 * dsh_40[k]
                   + f_4 * pc_x[k] * dph_82[k];

        t_104[k] = pa_y[k] * ppi0_20[k]
                   - f_11 * pc_y[k] * ppi1_20[k];

        t_105[k] = f_1 * pph_15[k]
                   + f_2 * dpg0_55[k]
                   - f_3 * dpg1_55[k]
                   + f_4 * pc_y[k] * dph_78[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pc_z, dpg0_55, dpg0_56, dpg0_57, dpg1_55, \
                         dpg1_56, dpg1_57, dph_78, dph_79, dph_80, \
                         dph_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_4 * pc_z[k] * dph_78[k];

        t_107[k] = f_5 * dpg0_55[k]
                   - f_6 * dpg1_55[k]
                   + f_4 * pc_z[k] * dph_79[k];

        t_108[k] = f_7 * dpg0_56[k]
                   - f_8 * dpg1_56[k]
                   + f_4 * pc_z[k] * dph_80[k];

        t_109[k] = f_9 * dpg0_57[k]
                   - f_10 * dpg1_57[k]
                   + f_4 * pc_z[k] * dph_81[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pc_x, pc_y, pc_z, pph_20, pph_84, dpg0_59, \
                         dpg0_60, dpg1_59, dpg1_60, dph_83, dph_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_1 * pph_20[k]
                   + f_4 * pc_y[k] * dph_83[k];

        t_111[k] = f_2 * dpg0_59[k]
                   - f_3 * dpg1_59[k]
                   + f_4 * pc_z[k] * dph_83[k];

        t_112[k] = f_1 * pph_84[k]
                   + f_2 * dpg0_60[k]
                   - f_3 * dpg1_60[k]
                   + f_4 * pc_x[k] * dph_84[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pc_x, pc_y, pc_z, ppi0_115, pph_21, \
                         pph_87, ppi1_115, dsh_21, dph_84, dph_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_1 * pph_21[k]
                   + f_1 * dsh_21[k]
                   + f_4 * pc_y[k] * dph_84[k];

        t_114[k] = f_4 * pc_z[k] * dph_84[k];

        t_115[k] = pa_x[k] * ppi0_115[k]
                   + f_13 * pph_87[k]
                   - f_11 * pc_x[k] * ppi1_115[k];

        t_116[k] = f_4 * pc_z[k] * dph_85[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pc_x, pc_z, ppi0_118, pph_90, ppi1_118, \
                         dpg0_60, dpg1_60, dph_86, dph_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_5 * dpg0_60[k]
                   - f_6 * dpg1_60[k]
                   + f_4 * pc_z[k] * dph_86[k];

        t_118[k] = pa_x[k] * ppi0_118[k]
                   + f_12 * pph_90[k]
                   - f_11 * pc_x[k] * ppi1_118[k];

        t_119[k] = f_4 * pc_z[k] * dph_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_x, pc_x, pc_y, pc_z, ppi0_122, pph_26, \
                         pph_94, ppi1_122, dsh_26, dpg0_62, dpg1_62, \
                         dph_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_1 * pph_26[k]
                   + f_1 * dsh_26[k]
                   + f_4 * pc_y[k] * dph_89[k];

        t_121[k] = f_7 * dpg0_62[k]
                   - f_8 * dpg1_62[k]
                   + f_4 * pc_z[k] * dph_89[k];

        t_122[k] = pa_x[k] * ppi0_122[k]
                   + f_0 * pph_94[k]
                   - f_11 * pc_x[k] * ppi1_122[k];
    }
}

static auto
compute_prim_dpi_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t pph,
                                                          const size_t ppi1, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t dph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_0 = buffer.data(ppi0 + 0);
    const auto *ppi0_3 = buffer.data(ppi0 + 3);
    const auto *ppi0_6 = buffer.data(ppi0 + 6);
    const auto *ppi0_10 = buffer.data(ppi0 + 10);
    const auto *ppi0_15 = buffer.data(ppi0 + 15);
    const auto *ppi0_28 = buffer.data(ppi0 + 28);
    const auto *ppi0_31 = buffer.data(ppi0 + 31);
    const auto *ppi0_34 = buffer.data(ppi0 + 34);
    const auto *ppi0_38 = buffer.data(ppi0 + 38);
    const auto *ppi0_56 = buffer.data(ppi0 + 56);
    const auto *ppi0_61 = buffer.data(ppi0 + 61);
    const auto *ppi0_65 = buffer.data(ppi0 + 65);
    const auto *ppi0_70 = buffer.data(ppi0 + 70);
    const auto *ppi0_133 = buffer.data(ppi0 + 133);
    const auto *ppi0_135 = buffer.data(ppi0 + 135);
    const auto *ppi0_136 = buffer.data(ppi0 + 136);
    const auto *ppi0_137 = buffer.data(ppi0 + 137);
    const auto *ppi0_138 = buffer.data(ppi0 + 138);
    const auto *ppi0_139 = buffer.data(ppi0 + 139);
    const auto *ppi0_152 = buffer.data(ppi0 + 152);
    const auto *ppi0_161 = buffer.data(ppi0 + 161);
    const auto *ppi0_163 = buffer.data(ppi0 + 163);
    const auto *ppi0_164 = buffer.data(ppi0 + 164);
    const auto *ppi0_165 = buffer.data(ppi0 + 165);
    const auto *ppi0_167 = buffer.data(ppi0 + 167);
    const auto *ppi0_203 = buffer.data(ppi0 + 203);
    const auto *ppi0_207 = buffer.data(ppi0 + 207);
    const auto *ppi0_208 = buffer.data(ppi0 + 208);
    const auto *ppi0_217 = buffer.data(ppi0 + 217);
    const auto *ppi0_218 = buffer.data(ppi0 + 218);
    const auto *ppi0_219 = buffer.data(ppi0 + 219);
    const auto *ppi0_220 = buffer.data(ppi0 + 220);
    const auto *ppi0_221 = buffer.data(ppi0 + 221);
    const auto *ppi0_223 = buffer.data(ppi0 + 223);
    const auto *ppi0_229 = buffer.data(ppi0 + 229);
    const auto *ppi0_233 = buffer.data(ppi0 + 233);
    const auto *ppi0_238 = buffer.data(ppi0 + 238);
    const auto *ppi0_245 = buffer.data(ppi0 + 245);

    const auto *pph_0 = buffer.data(pph + 0);
    const auto *pph_20 = buffer.data(pph + 20);
    const auto *pph_21 = buffer.data(pph + 21);
    const auto *pph_30 = buffer.data(pph + 30);
    const auto *pph_42 = buffer.data(pph + 42);
    const auto *pph_47 = buffer.data(pph + 47);
    const auto *pph_51 = buffer.data(pph + 51);
    const auto *pph_62 = buffer.data(pph + 62);
    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_101 = buffer.data(pph + 101);
    const auto *pph_102 = buffer.data(pph + 102);
    const auto *pph_103 = buffer.data(pph + 103);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_117 = buffer.data(pph + 117);
    const auto *pph_120 = buffer.data(pph + 120);
    const auto *pph_122 = buffer.data(pph + 122);
    const auto *pph_123 = buffer.data(pph + 123);
    const auto *pph_124 = buffer.data(pph + 124);
    const auto *pph_125 = buffer.data(pph + 125);
    const auto *pph_131 = buffer.data(pph + 131);
    const auto *pph_135 = buffer.data(pph + 135);
    const auto *pph_140 = buffer.data(pph + 140);
    const auto *pph_142 = buffer.data(pph + 142);
    const auto *pph_143 = buffer.data(pph + 143);
    const auto *pph_144 = buffer.data(pph + 144);
    const auto *pph_146 = buffer.data(pph + 146);
    const auto *pph_154 = buffer.data(pph + 154);
    const auto *pph_158 = buffer.data(pph + 158);
    const auto *pph_159 = buffer.data(pph + 159);
    const auto *pph_162 = buffer.data(pph + 162);
    const auto *pph_163 = buffer.data(pph + 163);
    const auto *pph_164 = buffer.data(pph + 164);
    const auto *pph_165 = buffer.data(pph + 165);
    const auto *pph_167 = buffer.data(pph + 167);
    const auto *pph_168 = buffer.data(pph + 168);
    const auto *pph_173 = buffer.data(pph + 173);
    const auto *pph_177 = buffer.data(pph + 177);
    const auto *pph_182 = buffer.data(pph + 182);
    const auto *pph_183 = buffer.data(pph + 183);
    const auto *pph_184 = buffer.data(pph + 184);
    const auto *pph_185 = buffer.data(pph + 185);
    const auto *pph_186 = buffer.data(pph + 186);
    const auto *pph_188 = buffer.data(pph + 188);

    const auto *ppi1_0 = buffer.data(ppi1 + 0);
    const auto *ppi1_3 = buffer.data(ppi1 + 3);
    const auto *ppi1_6 = buffer.data(ppi1 + 6);
    const auto *ppi1_10 = buffer.data(ppi1 + 10);
    const auto *ppi1_15 = buffer.data(ppi1 + 15);
    const auto *ppi1_28 = buffer.data(ppi1 + 28);
    const auto *ppi1_31 = buffer.data(ppi1 + 31);
    const auto *ppi1_34 = buffer.data(ppi1 + 34);
    const auto *ppi1_38 = buffer.data(ppi1 + 38);
    const auto *ppi1_56 = buffer.data(ppi1 + 56);
    const auto *ppi1_61 = buffer.data(ppi1 + 61);
    const auto *ppi1_65 = buffer.data(ppi1 + 65);
    const auto *ppi1_70 = buffer.data(ppi1 + 70);
    const auto *ppi1_133 = buffer.data(ppi1 + 133);
    const auto *ppi1_135 = buffer.data(ppi1 + 135);
    const auto *ppi1_136 = buffer.data(ppi1 + 136);
    const auto *ppi1_137 = buffer.data(ppi1 + 137);
    const auto *ppi1_138 = buffer.data(ppi1 + 138);
    const auto *ppi1_139 = buffer.data(ppi1 + 139);
    const auto *ppi1_152 = buffer.data(ppi1 + 152);
    const auto *ppi1_161 = buffer.data(ppi1 + 161);
    const auto *ppi1_163 = buffer.data(ppi1 + 163);
    const auto *ppi1_164 = buffer.data(ppi1 + 164);
    const auto *ppi1_165 = buffer.data(ppi1 + 165);
    const auto *ppi1_167 = buffer.data(ppi1 + 167);
    const auto *ppi1_203 = buffer.data(ppi1 + 203);
    const auto *ppi1_207 = buffer.data(ppi1 + 207);
    const auto *ppi1_208 = buffer.data(ppi1 + 208);
    const auto *ppi1_217 = buffer.data(ppi1 + 217);
    const auto *ppi1_218 = buffer.data(ppi1 + 218);
    const auto *ppi1_219 = buffer.data(ppi1 + 219);
    const auto *ppi1_220 = buffer.data(ppi1 + 220);
    const auto *ppi1_221 = buffer.data(ppi1 + 221);
    const auto *ppi1_223 = buffer.data(ppi1 + 223);
    const auto *ppi1_229 = buffer.data(ppi1 + 229);
    const auto *ppi1_233 = buffer.data(ppi1 + 233);
    const auto *ppi1_238 = buffer.data(ppi1 + 238);
    const auto *ppi1_245 = buffer.data(ppi1 + 245);

    const auto *dsi0_31 = buffer.data(dsi0 + 31);
    const auto *dsi0_34 = buffer.data(dsi0 + 34);
    const auto *dsi0_38 = buffer.data(dsi0 + 38);
    const auto *dsi0_61 = buffer.data(dsi0 + 61);
    const auto *dsi0_65 = buffer.data(dsi0 + 65);
    const auto *dsi0_70 = buffer.data(dsi0 + 70);

    const auto *dsh_21 = buffer.data(dsh + 21);
    const auto *dsh_22 = buffer.data(dsh + 22);
    const auto *dsh_24 = buffer.data(dsh + 24);
    const auto *dsh_27 = buffer.data(dsh + 27);
    const auto *dsh_30 = buffer.data(dsh + 30);
    const auto *dsh_31 = buffer.data(dsh + 31);
    const auto *dsh_36 = buffer.data(dsh + 36);
    const auto *dsh_42 = buffer.data(dsh + 42);
    const auto *dsh_44 = buffer.data(dsh + 44);
    const auto *dsh_47 = buffer.data(dsh + 47);
    const auto *dsh_51 = buffer.data(dsh + 51);
    const auto *dsh_56 = buffer.data(dsh + 56);
    const auto *dsh_58 = buffer.data(dsh + 58);
    const auto *dsh_59 = buffer.data(dsh + 59);
    const auto *dsh_60 = buffer.data(dsh + 60);
    const auto *dsh_62 = buffer.data(dsh + 62);

    const auto *dsi1_31 = buffer.data(dsi1 + 31);
    const auto *dsi1_34 = buffer.data(dsi1 + 34);
    const auto *dsi1_38 = buffer.data(dsi1 + 38);
    const auto *dsi1_61 = buffer.data(dsi1 + 61);
    const auto *dsi1_65 = buffer.data(dsi1 + 65);
    const auto *dsi1_70 = buffer.data(dsi1 + 70);

    const auto *dpg0_63 = buffer.data(dpg0 + 63);
    const auto *dpg0_65 = buffer.data(dpg0 + 65);
    const auto *dpg0_92 = buffer.data(dpg0 + 92);
    const auto *dpg0_94 = buffer.data(dpg0 + 94);
    const auto *dpg0_95 = buffer.data(dpg0 + 95);
    const auto *dpg0_99 = buffer.data(dpg0 + 99);
    const auto *dpg0_100 = buffer.data(dpg0 + 100);
    const auto *dpg0_101 = buffer.data(dpg0 + 101);
    const auto *dpg0_102 = buffer.data(dpg0 + 102);
    const auto *dpg0_103 = buffer.data(dpg0 + 103);
    const auto *dpg0_104 = buffer.data(dpg0 + 104);
    const auto *dpg0_120 = buffer.data(dpg0 + 120);
    const auto *dpg0_121 = buffer.data(dpg0 + 121);
    const auto *dpg0_122 = buffer.data(dpg0 + 122);
    const auto *dpg0_123 = buffer.data(dpg0 + 123);
    const auto *dpg0_124 = buffer.data(dpg0 + 124);
    const auto *dpg0_125 = buffer.data(dpg0 + 125);

    const auto *dpg1_63 = buffer.data(dpg1 + 63);
    const auto *dpg1_65 = buffer.data(dpg1 + 65);
    const auto *dpg1_92 = buffer.data(dpg1 + 92);
    const auto *dpg1_94 = buffer.data(dpg1 + 94);
    const auto *dpg1_95 = buffer.data(dpg1 + 95);
    const auto *dpg1_99 = buffer.data(dpg1 + 99);
    const auto *dpg1_100 = buffer.data(dpg1 + 100);
    const auto *dpg1_101 = buffer.data(dpg1 + 101);
    const auto *dpg1_102 = buffer.data(dpg1 + 102);
    const auto *dpg1_103 = buffer.data(dpg1 + 103);
    const auto *dpg1_104 = buffer.data(dpg1 + 104);
    const auto *dpg1_120 = buffer.data(dpg1 + 120);
    const auto *dpg1_121 = buffer.data(dpg1 + 121);
    const auto *dpg1_122 = buffer.data(dpg1 + 122);
    const auto *dpg1_123 = buffer.data(dpg1 + 123);
    const auto *dpg1_124 = buffer.data(dpg1 + 124);
    const auto *dpg1_125 = buffer.data(dpg1 + 125);

    const auto *dph_90 = buffer.data(dph + 90);
    const auto *dph_91 = buffer.data(dph + 91);
    const auto *dph_93 = buffer.data(dph + 93);
    const auto *dph_94 = buffer.data(dph + 94);
    const auto *dph_99 = buffer.data(dph + 99);
    const auto *dph_101 = buffer.data(dph + 101);
    const auto *dph_102 = buffer.data(dph + 102);
    const auto *dph_103 = buffer.data(dph + 103);
    const auto *dph_104 = buffer.data(dph + 104);
    const auto *dph_105 = buffer.data(dph + 105);
    const auto *dph_106 = buffer.data(dph + 106);
    const auto *dph_108 = buffer.data(dph + 108);
    const auto *dph_110 = buffer.data(dph + 110);
    const auto *dph_111 = buffer.data(dph + 111);
    const auto *dph_114 = buffer.data(dph + 114);
    const auto *dph_115 = buffer.data(dph + 115);
    const auto *dph_120 = buffer.data(dph + 120);
    const auto *dph_122 = buffer.data(dph + 122);
    const auto *dph_123 = buffer.data(dph + 123);
    const auto *dph_124 = buffer.data(dph + 124);
    const auto *dph_125 = buffer.data(dph + 125);
    const auto *dph_126 = buffer.data(dph + 126);
    const auto *dph_128 = buffer.data(dph + 128);
    const auto *dph_130 = buffer.data(dph + 130);
    const auto *dph_131 = buffer.data(dph + 131);
    const auto *dph_133 = buffer.data(dph + 133);
    const auto *dph_134 = buffer.data(dph + 134);
    const auto *dph_135 = buffer.data(dph + 135);
    const auto *dph_140 = buffer.data(dph + 140);
    const auto *dph_141 = buffer.data(dph + 141);
    const auto *dph_142 = buffer.data(dph + 142);
    const auto *dph_143 = buffer.data(dph + 143);
    const auto *dph_144 = buffer.data(dph + 144);
    const auto *dph_145 = buffer.data(dph + 145);
    const auto *dph_146 = buffer.data(dph + 146);
    const auto *dph_147 = buffer.data(dph + 147);
    const auto *dph_149 = buffer.data(dph + 149);
    const auto *dph_152 = buffer.data(dph + 152);
    const auto *dph_156 = buffer.data(dph + 156);
    const auto *dph_161 = buffer.data(dph + 161);
    const auto *dph_162 = buffer.data(dph + 162);
    const auto *dph_163 = buffer.data(dph + 163);
    const auto *dph_164 = buffer.data(dph + 164);
    const auto *dph_165 = buffer.data(dph + 165);
    const auto *dph_167 = buffer.data(dph + 167);
    const auto *dph_168 = buffer.data(dph + 168);
    const auto *dph_169 = buffer.data(dph + 169);
    const auto *dph_170 = buffer.data(dph + 170);
    const auto *dph_171 = buffer.data(dph + 171);
    const auto *dph_172 = buffer.data(dph + 172);
    const auto *dph_173 = buffer.data(dph + 173);
    const auto *dph_174 = buffer.data(dph + 174);
    const auto *dph_175 = buffer.data(dph + 175);
    const auto *dph_176 = buffer.data(dph + 176);
    const auto *dph_177 = buffer.data(dph + 177);
    const auto *dph_182 = buffer.data(dph + 182);
    const auto *dph_183 = buffer.data(dph + 183);
    const auto *dph_184 = buffer.data(dph + 184);
    const auto *dph_185 = buffer.data(dph + 185);
    const auto *dph_186 = buffer.data(dph + 186);
    const auto *dph_188 = buffer.data(dph + 188);

#pragma omp simd aligned(t_123, t_124, t_125, t_126, pc_y, pc_z, pph_30, dsh_30, dpg0_63, \
                         dpg0_65, dpg1_63, dpg1_65, dph_90, dph_91, \
                         dph_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_4 * pc_z[k] * dph_90[k];

        t_124[k] = f_5 * dpg0_63[k]
                   - f_6 * dpg1_63[k]
                   + f_4 * pc_z[k] * dph_91[k];

        t_125[k] = f_1 * pph_30[k]
                   + f_1 * dsh_30[k]
                   + f_4 * pc_y[k] * dph_93[k];

        t_126[k] = f_9 * dpg0_65[k]
                   - f_10 * dpg1_65[k]
                   + f_4 * pc_z[k] * dph_93[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, t_131, pc_x, pc_z, pph_99, pph_101, \
                         pph_102, pph_103, dph_94, dph_99, dph_101, dph_102, \
                         dph_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_1 * pph_99[k]
                   + f_4 * pc_x[k] * dph_99[k];

        t_128[k] = f_4 * pc_z[k] * dph_94[k];

        t_129[k] = f_1 * pph_101[k]
                   + f_4 * pc_x[k] * dph_101[k];

        t_130[k] = f_1 * pph_102[k]
                   + f_4 * pc_x[k] * dph_102[k];

        t_131[k] = f_1 * pph_103[k]
                   + f_4 * pc_x[k] * dph_103[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_x, pc_x, pc_z, ppi0_133, ppi0_135, \
                         pph_104, ppi1_133, ppi1_135, dph_99, dph_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_1 * pph_104[k]
                   + f_4 * pc_x[k] * dph_104[k];

        t_133[k] = pa_x[k] * ppi0_133[k]
                   - f_11 * pc_x[k] * ppi1_133[k];

        t_134[k] = f_4 * pc_z[k] * dph_99[k];

        t_135[k] = pa_x[k] * ppi0_135[k]
                   - f_11 * pc_x[k] * ppi1_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, pc_x, ppi0_136, ppi0_137, ppi0_138, \
                         ppi0_139, ppi1_136, ppi1_137, ppi1_138, \
                         ppi1_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * ppi0_136[k]
                   - f_11 * pc_x[k] * ppi1_136[k];

        t_137[k] = pa_x[k] * ppi0_137[k]
                   - f_11 * pc_x[k] * ppi1_137[k];

        t_138[k] = pa_x[k] * ppi0_138[k]
                   - f_11 * pc_x[k] * ppi1_138[k];

        t_139[k] = pa_x[k] * ppi0_139[k]
                   - f_11 * pc_x[k] * ppi1_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_y, pb_z, pc_y, pc_z, ppi0_56, pph_42, \
                         ppi1_56, dsi0_31, dsh_21, dsi1_31, dph_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pa_y[k] * ppi0_56[k]
                   - f_11 * pc_y[k] * ppi1_56[k];

        t_141[k] = f_1 * pph_42[k]
                   + f_4 * pc_y[k] * dph_105[k];

        t_142[k] = f_1 * dsh_21[k]
                   + f_4 * pc_z[k] * dph_105[k];

        t_143[k] = pb_z[k] * dsi0_31[k]
                   - f_11 * pc_z[k] * dsi1_31[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_y, pb_z, pc_y, pc_z, ppi0_61, ppi1_61, \
                         dsi0_34, dsh_22, dsh_24, dsi1_34, dph_106, \
                         dph_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_1 * dsh_22[k]
                   + f_4 * pc_z[k] * dph_106[k];

        t_145[k] = pa_y[k] * ppi0_61[k]
                   - f_11 * pc_y[k] * ppi1_61[k];

        t_146[k] = pb_z[k] * dsi0_34[k]
                   - f_11 * pc_z[k] * dsi1_34[k];

        t_147[k] = f_1 * dsh_24[k]
                   + f_4 * pc_z[k] * dph_108[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_y, pb_z, pc_y, pc_z, ppi0_65, pph_47, \
                         ppi1_65, dsi0_38, dsh_27, dsi1_38, dph_110, \
                         dph_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_1 * pph_47[k]
                   + f_4 * pc_y[k] * dph_110[k];

        t_149[k] = pa_y[k] * ppi0_65[k]
                   - f_11 * pc_y[k] * ppi1_65[k];

        t_150[k] = pb_z[k] * dsi0_38[k]
                   - f_11 * pc_z[k] * dsi1_38[k];

        t_151[k] = f_1 * dsh_27[k]
                   + f_4 * pc_z[k] * dph_111[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, pa_x, pa_y, pc_x, pc_y, ppi0_70, ppi0_152, \
                         pph_51, pph_117, ppi1_70, ppi1_152, dph_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = pa_x[k] * ppi0_152[k]
                   + f_0 * pph_117[k]
                   - f_11 * pc_x[k] * ppi1_152[k];

        t_153[k] = f_1 * pph_51[k]
                   + f_4 * pc_y[k] * dph_114[k];

        t_154[k] = pa_y[k] * ppi0_70[k]
                   - f_11 * pc_y[k] * ppi1_70[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pc_x, pc_z, pph_120, pph_122, pph_123, \
                         dsh_31, dph_115, dph_120, dph_122, dph_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_1 * pph_120[k]
                   + f_4 * pc_x[k] * dph_120[k];

        t_156[k] = f_1 * dsh_31[k]
                   + f_4 * pc_z[k] * dph_115[k];

        t_157[k] = f_1 * pph_122[k]
                   + f_4 * pc_x[k] * dph_122[k];

        t_158[k] = f_1 * pph_123[k]
                   + f_4 * pc_x[k] * dph_123[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pc_x, pc_z, ppi0_161, pph_124, \
                         pph_125, ppi1_161, dsh_36, dph_120, dph_124, \
                         dph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_1 * pph_124[k]
                   + f_4 * pc_x[k] * dph_124[k];

        t_160[k] = f_1 * pph_125[k]
                   + f_4 * pc_x[k] * dph_125[k];

        t_161[k] = pa_x[k] * ppi0_161[k]
                   - f_11 * pc_x[k] * ppi1_161[k];

        t_162[k] = f_1 * dsh_36[k]
                   + f_4 * pc_z[k] * dph_120[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_x, pc_x, pc_y, ppi0_163, ppi0_164, \
                         ppi0_165, pph_62, ppi1_163, ppi1_164, ppi1_165, \
                         dph_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * ppi0_163[k]
                   - f_11 * pc_x[k] * ppi1_163[k];

        t_164[k] = pa_x[k] * ppi0_164[k]
                   - f_11 * pc_x[k] * ppi1_164[k];

        t_165[k] = pa_x[k] * ppi0_165[k]
                   - f_11 * pc_x[k] * ppi1_165[k];

        t_166[k] = f_1 * pph_62[k]
                   + f_4 * pc_y[k] * dph_125[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pa_x, pa_z, pc_x, pc_y, pc_z, ppi0_0, \
                         ppi0_167, pph_0, ppi1_0, ppi1_167, dph_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = pa_x[k] * ppi0_167[k]
                   - f_11 * pc_x[k] * ppi1_167[k];

        t_168[k] = pa_z[k] * ppi0_0[k]
                   - f_11 * pc_z[k] * ppi1_0[k];

        t_169[k] = f_4 * pc_y[k] * dph_126[k];

        t_170[k] = f_1 * pph_0[k]
                   + f_4 * pc_z[k] * dph_126[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_z, pc_x, pc_y, pc_z, ppi0_3, pph_131, ppi1_3, \
                         dsh_47, dpg0_95, dpg1_95, dph_128, dph_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = pa_z[k] * ppi0_3[k]
                   - f_11 * pc_z[k] * ppi1_3[k];

        t_172[k] = f_4 * pc_y[k] * dph_128[k];

        t_173[k] = f_1 * pph_131[k]
                   + f_1 * dsh_47[k]
                   + f_9 * dpg0_95[k]
                   - f_10 * dpg1_95[k]
                   + f_4 * pc_x[k] * dph_131[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_z, pc_y, pc_z, ppi0_6, ppi1_6, dpg0_92, \
                         dpg1_92, dph_130, dph_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = pa_z[k] * ppi0_6[k]
                   - f_11 * pc_z[k] * ppi1_6[k];

        t_175[k] = f_5 * dpg0_92[k]
                   - f_6 * dpg1_92[k]
                   + f_4 * pc_y[k] * dph_130[k];

        t_176[k] = f_4 * pc_y[k] * dph_131[k];
    }

#pragma omp simd aligned(t_177, t_178, pa_z, pc_x, pc_z, ppi0_10, pph_135, ppi1_10, dsh_51, \
                         dpg0_99, dpg1_99, dph_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_1 * pph_135[k]
                   + f_1 * dsh_51[k]
                   + f_7 * dpg0_99[k]
                   - f_8 * dpg1_99[k]
                   + f_4 * pc_x[k] * dph_135[k];

        t_178[k] = pa_z[k] * ppi0_10[k]
                   - f_11 * pc_z[k] * ppi1_10[k];
    }

#pragma omp simd aligned(t_179, t_180, t_181, pc_y, dpg0_94, dpg0_95, dpg1_94, dpg1_95, \
                         dph_133, dph_134, dph_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_179[k] = f_7 * dpg0_94[k]
                   - f_8 * dpg1_94[k]
                   + f_4 * pc_y[k] * dph_133[k];

        t_180[k] = f_5 * dpg0_95[k]
                   - f_6 * dpg1_95[k]
                   + f_4 * pc_y[k] * dph_134[k];

        t_181[k] = f_4 * pc_y[k] * dph_135[k];
    }

#pragma omp simd aligned(t_182, t_183, t_184, pa_z, pc_x, pc_z, ppi0_15, pph_140, pph_142, \
                         ppi1_15, dsh_56, dsh_58, dpg0_104, dpg1_104, dph_140, \
                         dph_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_182[k] = f_1 * pph_140[k]
                   + f_1 * dsh_56[k]
                   + f_5 * dpg0_104[k]
                   - f_6 * dpg1_104[k]
                   + f_4 * pc_x[k] * dph_140[k];

        t_183[k] = pa_z[k] * ppi0_15[k]
                   - f_11 * pc_z[k] * ppi1_15[k];

        t_184[k] = f_1 * pph_142[k]
                   + f_1 * dsh_58[k]
                   + f_4 * pc_x[k] * dph_142[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pc_x, pc_y, pph_143, pph_144, pph_146, \
                         dsh_59, dsh_60, dsh_62, dph_140, dph_143, dph_144, \
                         dph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_1 * pph_143[k]
                   + f_1 * dsh_59[k]
                   + f_4 * pc_x[k] * dph_143[k];

        t_186[k] = f_1 * pph_144[k]
                   + f_1 * dsh_60[k]
                   + f_4 * pc_x[k] * dph_144[k];

        t_187[k] = f_4 * pc_y[k] * dph_140[k];

        t_188[k] = f_1 * pph_146[k]
                   + f_1 * dsh_62[k]
                   + f_4 * pc_x[k] * dph_146[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pc_y, dpg0_100, dpg0_101, dpg0_102, dpg1_100, \
                         dpg1_101, dpg1_102, dph_141, dph_142, \
                         dph_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_2 * dpg0_100[k]
                   - f_3 * dpg1_100[k]
                   + f_4 * pc_y[k] * dph_141[k];

        t_190[k] = f_15 * dpg0_101[k]
                   - f_16 * dpg1_101[k]
                   + f_4 * pc_y[k] * dph_142[k];

        t_191[k] = f_9 * dpg0_102[k]
                   - f_10 * dpg1_102[k]
                   + f_4 * pc_y[k] * dph_143[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, t_195, pc_y, pc_z, pph_20, dpg0_103, dpg0_104, \
                         dpg1_103, dpg1_104, dph_144, dph_145, \
                         dph_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_7 * dpg0_103[k]
                   - f_8 * dpg1_103[k]
                   + f_4 * pc_y[k] * dph_144[k];

        t_193[k] = f_5 * dpg0_104[k]
                   - f_6 * dpg1_104[k]
                   + f_4 * pc_y[k] * dph_145[k];

        t_194[k] = f_4 * pc_y[k] * dph_146[k];

        t_195[k] = f_1 * pph_20[k]
                   + f_2 * dpg0_104[k]
                   - f_3 * dpg1_104[k]
                   + f_4 * pc_z[k] * dph_146[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_z, pc_y, pc_z, ppi0_28, ppi0_31, \
                         pph_21, ppi1_28, ppi1_31, dsh_42, dph_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_z[k] * ppi0_28[k]
                   - f_11 * pc_z[k] * ppi1_28[k];

        t_197[k] = f_1 * dsh_42[k]
                   + f_4 * pc_y[k] * dph_147[k];

        t_198[k] = f_1 * pph_21[k]
                   + f_4 * pc_z[k] * dph_147[k];

        t_199[k] = pa_z[k] * ppi0_31[k]
                   - f_11 * pc_z[k] * ppi1_31[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pa_z, pb_y, pc_y, pc_z, ppi0_34, ppi1_34, \
                         dsi0_61, dsh_44, dsi1_61, dph_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_1 * dsh_44[k]
                   + f_4 * pc_y[k] * dph_149[k];

        t_201[k] = pb_y[k] * dsi0_61[k]
                   - f_11 * pc_y[k] * dsi1_61[k];

        t_202[k] = pa_z[k] * ppi0_34[k]
                   - f_11 * pc_z[k] * ppi1_34[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pa_x, pb_y, pc_x, pc_y, ppi0_203, pph_154, \
                         ppi1_203, dsi0_65, dsh_47, dsi1_65, dph_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = pa_x[k] * ppi0_203[k]
                   + f_12 * pph_154[k]
                   - f_11 * pc_x[k] * ppi1_203[k];

        t_204[k] = f_1 * dsh_47[k]
                   + f_4 * pc_y[k] * dph_152[k];

        t_205[k] = pb_y[k] * dsi0_65[k]
                   - f_11 * pc_y[k] * dsi1_65[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_x, pa_z, pc_x, pc_z, ppi0_38, ppi0_207, \
                         ppi0_208, pph_158, pph_159, ppi1_38, ppi1_207, \
                         ppi1_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_z[k] * ppi0_38[k]
                   - f_11 * pc_z[k] * ppi1_38[k];

        t_207[k] = pa_x[k] * ppi0_207[k]
                   + f_0 * pph_158[k]
                   - f_11 * pc_x[k] * ppi1_207[k];

        t_208[k] = pa_x[k] * ppi0_208[k]
                   + f_0 * pph_159[k]
                   - f_11 * pc_x[k] * ppi1_208[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, t_212, pb_y, pc_x, pc_y, pph_162, pph_163, \
                         dsi0_70, dsh_51, dsi1_70, dph_156, dph_162, \
                         dph_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = f_1 * dsh_51[k]
                   + f_4 * pc_y[k] * dph_156[k];

        t_210[k] = pb_y[k] * dsi0_70[k]
                   - f_11 * pc_y[k] * dsi1_70[k];

        t_211[k] = f_1 * pph_162[k]
                   + f_4 * pc_x[k] * dph_162[k];

        t_212[k] = f_1 * pph_163[k]
                   + f_4 * pc_x[k] * dph_163[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pc_x, pc_y, pph_164, pph_165, pph_167, \
                         dsh_56, dph_161, dph_164, dph_165, dph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_1 * pph_164[k]
                   + f_4 * pc_x[k] * dph_164[k];

        t_214[k] = f_1 * pph_165[k]
                   + f_4 * pc_x[k] * dph_165[k];

        t_215[k] = f_1 * dsh_56[k]
                   + f_4 * pc_y[k] * dph_161[k];

        t_216[k] = f_1 * pph_167[k]
                   + f_4 * pc_x[k] * dph_167[k];
    }

#pragma omp simd aligned(t_217, t_218, t_219, t_220, pa_x, pc_x, ppi0_217, ppi0_218, ppi0_219, \
                         ppi0_220, ppi1_217, ppi1_218, ppi1_219, \
                         ppi1_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = pa_x[k] * ppi0_217[k]
                   - f_11 * pc_x[k] * ppi1_217[k];

        t_218[k] = pa_x[k] * ppi0_218[k]
                   - f_11 * pc_x[k] * ppi1_218[k];

        t_219[k] = pa_x[k] * ppi0_219[k]
                   - f_11 * pc_x[k] * ppi1_219[k];

        t_220[k] = pa_x[k] * ppi0_220[k]
                   - f_11 * pc_x[k] * ppi1_220[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pa_x, pc_x, pc_y, ppi0_221, ppi0_223, ppi1_221, \
                         ppi1_223, dsh_62, dph_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = pa_x[k] * ppi0_221[k]
                   - f_11 * pc_x[k] * ppi1_221[k];

        t_222[k] = f_1 * dsh_62[k]
                   + f_4 * pc_y[k] * dph_167[k];

        t_223[k] = pa_x[k] * ppi0_223[k]
                   - f_11 * pc_x[k] * ppi1_223[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pc_x, pc_y, pc_z, pph_42, pph_168, \
                         dsh_42, dpg0_120, dpg1_120, dph_168, dph_169, \
                         dph_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_1 * pph_168[k]
                   + f_2 * dpg0_120[k]
                   - f_3 * dpg1_120[k]
                   + f_4 * pc_x[k] * dph_168[k];

        t_225[k] = f_4 * pc_y[k] * dph_168[k];

        t_226[k] = f_1 * pph_42[k]
                   + f_1 * dsh_42[k]
                   + f_4 * pc_z[k] * dph_168[k];

        t_227[k] = f_5 * dpg0_120[k]
                   - f_6 * dpg1_120[k]
                   + f_4 * pc_y[k] * dph_169[k];

        t_228[k] = f_4 * pc_y[k] * dph_170[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, pa_x, pc_x, pc_y, ppi0_229, pph_173, ppi1_229, \
                         dpg0_121, dpg0_122, dpg1_121, dpg1_122, dph_171, \
                         dph_172 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = pa_x[k] * ppi0_229[k]
                   + f_13 * pph_173[k]
                   - f_11 * pc_x[k] * ppi1_229[k];

        t_230[k] = f_7 * dpg0_121[k]
                   - f_8 * dpg1_121[k]
                   + f_4 * pc_y[k] * dph_171[k];

        t_231[k] = f_5 * dpg0_122[k]
                   - f_6 * dpg1_122[k]
                   + f_4 * pc_y[k] * dph_172[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pa_x, pc_x, pc_y, ppi0_233, pph_177, ppi1_233, \
                         dpg0_123, dpg1_123, dph_173, dph_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = f_4 * pc_y[k] * dph_173[k];

        t_233[k] = pa_x[k] * ppi0_233[k]
                   + f_12 * pph_177[k]
                   - f_11 * pc_x[k] * ppi1_233[k];

        t_234[k] = f_9 * dpg0_123[k]
                   - f_10 * dpg1_123[k]
                   + f_4 * pc_y[k] * dph_174[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pc_y, dpg0_124, dpg0_125, dpg1_124, dpg1_125, \
                         dph_175, dph_176, dph_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_7 * dpg0_124[k]
                   - f_8 * dpg1_124[k]
                   + f_4 * pc_y[k] * dph_175[k];

        t_236[k] = f_5 * dpg0_125[k]
                   - f_6 * dpg1_125[k]
                   + f_4 * pc_y[k] * dph_176[k];

        t_237[k] = f_4 * pc_y[k] * dph_177[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_x, pc_x, ppi0_238, pph_182, pph_183, \
                         pph_184, pph_185, ppi1_238, dph_183, dph_184, \
                         dph_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_x[k] * ppi0_238[k]
                   + f_0 * pph_182[k]
                   - f_11 * pc_x[k] * ppi1_238[k];

        t_239[k] = f_1 * pph_183[k]
                   + f_4 * pc_x[k] * dph_183[k];

        t_240[k] = f_1 * pph_184[k]
                   + f_4 * pc_x[k] * dph_184[k];

        t_241[k] = f_1 * pph_185[k]
                   + f_4 * pc_x[k] * dph_185[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, t_245, pa_x, pc_x, pc_y, ppi0_245, pph_186, \
                         pph_188, ppi1_245, dph_182, dph_186, dph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_1 * pph_186[k]
                   + f_4 * pc_x[k] * dph_186[k];

        t_243[k] = f_4 * pc_y[k] * dph_182[k];

        t_244[k] = f_1 * pph_188[k]
                   + f_4 * pc_x[k] * dph_188[k];

        t_245[k] = pa_x[k] * ppi0_245[k]
                   - f_11 * pc_x[k] * ppi1_245[k];
    }
}

static auto
compute_prim_dpi_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t pph,
                                                          const size_t ppi1, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t dph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 2.5 / q;

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
    auto *t_364 = buffer.data(target + 364);
    auto *t_365 = buffer.data(target + 365);
    auto *t_366 = buffer.data(target + 366);
    auto *t_367 = buffer.data(target + 367);
    auto *t_368 = buffer.data(target + 368);
    auto *t_369 = buffer.data(target + 369);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_85 = buffer.data(ppi0 + 85);
    const auto *ppi0_87 = buffer.data(ppi0 + 87);
    const auto *ppi0_90 = buffer.data(ppi0 + 90);
    const auto *ppi0_94 = buffer.data(ppi0 + 94);
    const auto *ppi0_105 = buffer.data(ppi0 + 105);
    const auto *ppi0_113 = buffer.data(ppi0 + 113);
    const auto *ppi0_115 = buffer.data(ppi0 + 115);
    const auto *ppi0_168 = buffer.data(ppi0 + 168);
    const auto *ppi0_170 = buffer.data(ppi0 + 170);
    const auto *ppi0_172 = buffer.data(ppi0 + 172);
    const auto *ppi0_173 = buffer.data(ppi0 + 173);
    const auto *ppi0_175 = buffer.data(ppi0 + 175);
    const auto *ppi0_176 = buffer.data(ppi0 + 176);
    const auto *ppi0_177 = buffer.data(ppi0 + 177);
    const auto *ppi0_179 = buffer.data(ppi0 + 179);
    const auto *ppi0_180 = buffer.data(ppi0 + 180);
    const auto *ppi0_181 = buffer.data(ppi0 + 181);
    const auto *ppi0_182 = buffer.data(ppi0 + 182);
    const auto *ppi0_195 = buffer.data(ppi0 + 195);
    const auto *ppi0_246 = buffer.data(ppi0 + 246);
    const auto *ppi0_247 = buffer.data(ppi0 + 247);
    const auto *ppi0_248 = buffer.data(ppi0 + 248);
    const auto *ppi0_249 = buffer.data(ppi0 + 249);
    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *pph_78 = buffer.data(pph + 78);
    const auto *pph_83 = buffer.data(pph + 83);
    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_125 = buffer.data(pph + 125);
    const auto *pph_128 = buffer.data(pph + 128);
    const auto *pph_130 = buffer.data(pph + 130);
    const auto *pph_131 = buffer.data(pph + 131);
    const auto *pph_133 = buffer.data(pph + 133);
    const auto *pph_134 = buffer.data(pph + 134);
    const auto *pph_135 = buffer.data(pph + 135);
    const auto *pph_146 = buffer.data(pph + 146);

    const auto *ppi1_85 = buffer.data(ppi1 + 85);
    const auto *ppi1_87 = buffer.data(ppi1 + 87);
    const auto *ppi1_90 = buffer.data(ppi1 + 90);
    const auto *ppi1_94 = buffer.data(ppi1 + 94);
    const auto *ppi1_105 = buffer.data(ppi1 + 105);
    const auto *ppi1_113 = buffer.data(ppi1 + 113);
    const auto *ppi1_115 = buffer.data(ppi1 + 115);
    const auto *ppi1_168 = buffer.data(ppi1 + 168);
    const auto *ppi1_170 = buffer.data(ppi1 + 170);
    const auto *ppi1_172 = buffer.data(ppi1 + 172);
    const auto *ppi1_173 = buffer.data(ppi1 + 173);
    const auto *ppi1_175 = buffer.data(ppi1 + 175);
    const auto *ppi1_176 = buffer.data(ppi1 + 176);
    const auto *ppi1_177 = buffer.data(ppi1 + 177);
    const auto *ppi1_179 = buffer.data(ppi1 + 179);
    const auto *ppi1_180 = buffer.data(ppi1 + 180);
    const auto *ppi1_181 = buffer.data(ppi1 + 181);
    const auto *ppi1_182 = buffer.data(ppi1 + 182);
    const auto *ppi1_195 = buffer.data(ppi1 + 195);
    const auto *ppi1_246 = buffer.data(ppi1 + 246);
    const auto *ppi1_247 = buffer.data(ppi1 + 247);
    const auto *ppi1_248 = buffer.data(ppi1 + 248);
    const auto *ppi1_249 = buffer.data(ppi1 + 249);
    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *dsi0_84 = buffer.data(dsi0 + 84);
    const auto *dsi0_85 = buffer.data(dsi0 + 85);
    const auto *dsi0_87 = buffer.data(dsi0 + 87);
    const auto *dsi0_89 = buffer.data(dsi0 + 89);
    const auto *dsi0_90 = buffer.data(dsi0 + 90);
    const auto *dsi0_92 = buffer.data(dsi0 + 92);
    const auto *dsi0_93 = buffer.data(dsi0 + 93);
    const auto *dsi0_94 = buffer.data(dsi0 + 94);
    const auto *dsi0_96 = buffer.data(dsi0 + 96);
    const auto *dsi0_97 = buffer.data(dsi0 + 97);
    const auto *dsi0_98 = buffer.data(dsi0 + 98);
    const auto *dsi0_105 = buffer.data(dsi0 + 105);
    const auto *dsi0_107 = buffer.data(dsi0 + 107);
    const auto *dsi0_108 = buffer.data(dsi0 + 108);
    const auto *dsi0_109 = buffer.data(dsi0 + 109);
    const auto *dsi0_111 = buffer.data(dsi0 + 111);
    const auto *dsi0_135 = buffer.data(dsi0 + 135);
    const auto *dsi0_136 = buffer.data(dsi0 + 136);
    const auto *dsi0_137 = buffer.data(dsi0 + 137);

    const auto *dsh_63 = buffer.data(dsh + 63);
    const auto *dsh_64 = buffer.data(dsh + 64);
    const auto *dsh_66 = buffer.data(dsh + 66);
    const auto *dsh_68 = buffer.data(dsh + 68);
    const auto *dsh_69 = buffer.data(dsh + 69);
    const auto *dsh_71 = buffer.data(dsh + 71);
    const auto *dsh_72 = buffer.data(dsh + 72);
    const auto *dsh_73 = buffer.data(dsh + 73);
    const auto *dsh_75 = buffer.data(dsh + 75);
    const auto *dsh_76 = buffer.data(dsh + 76);
    const auto *dsh_77 = buffer.data(dsh + 77);
    const auto *dsh_78 = buffer.data(dsh + 78);
    const auto *dsh_79 = buffer.data(dsh + 79);
    const auto *dsh_80 = buffer.data(dsh + 80);
    const auto *dsh_81 = buffer.data(dsh + 81);
    const auto *dsh_82 = buffer.data(dsh + 82);
    const auto *dsh_83 = buffer.data(dsh + 83);
    const auto *dsh_99 = buffer.data(dsh + 99);
    const auto *dsh_100 = buffer.data(dsh + 100);
    const auto *dsh_101 = buffer.data(dsh + 101);
    const auto *dsh_102 = buffer.data(dsh + 102);
    const auto *dsh_103 = buffer.data(dsh + 103);
    const auto *dsh_104 = buffer.data(dsh + 104);

    const auto *dsi1_84 = buffer.data(dsi1 + 84);
    const auto *dsi1_85 = buffer.data(dsi1 + 85);
    const auto *dsi1_87 = buffer.data(dsi1 + 87);
    const auto *dsi1_89 = buffer.data(dsi1 + 89);
    const auto *dsi1_90 = buffer.data(dsi1 + 90);
    const auto *dsi1_92 = buffer.data(dsi1 + 92);
    const auto *dsi1_93 = buffer.data(dsi1 + 93);
    const auto *dsi1_94 = buffer.data(dsi1 + 94);
    const auto *dsi1_96 = buffer.data(dsi1 + 96);
    const auto *dsi1_97 = buffer.data(dsi1 + 97);
    const auto *dsi1_98 = buffer.data(dsi1 + 98);
    const auto *dsi1_105 = buffer.data(dsi1 + 105);
    const auto *dsi1_107 = buffer.data(dsi1 + 107);
    const auto *dsi1_108 = buffer.data(dsi1 + 108);
    const auto *dsi1_109 = buffer.data(dsi1 + 109);
    const auto *dsi1_111 = buffer.data(dsi1 + 111);
    const auto *dsi1_135 = buffer.data(dsi1 + 135);
    const auto *dsi1_136 = buffer.data(dsi1 + 136);
    const auto *dsi1_137 = buffer.data(dsi1 + 137);

    const auto *dpg0_150 = buffer.data(dpg0 + 150);
    const auto *dpg0_151 = buffer.data(dpg0 + 151);
    const auto *dpg0_153 = buffer.data(dpg0 + 153);
    const auto *dpg0_155 = buffer.data(dpg0 + 155);
    const auto *dpg0_156 = buffer.data(dpg0 + 156);
    const auto *dpg0_158 = buffer.data(dpg0 + 158);
    const auto *dpg0_159 = buffer.data(dpg0 + 159);
    const auto *dpg0_160 = buffer.data(dpg0 + 160);
    const auto *dpg0_161 = buffer.data(dpg0 + 161);
    const auto *dpg0_162 = buffer.data(dpg0 + 162);
    const auto *dpg0_163 = buffer.data(dpg0 + 163);
    const auto *dpg0_164 = buffer.data(dpg0 + 164);
    const auto *dpg0_170 = buffer.data(dpg0 + 170);
    const auto *dpg0_173 = buffer.data(dpg0 + 173);
    const auto *dpg0_174 = buffer.data(dpg0 + 174);
    const auto *dpg0_177 = buffer.data(dpg0 + 177);
    const auto *dpg0_178 = buffer.data(dpg0 + 178);
    const auto *dpg0_179 = buffer.data(dpg0 + 179);
    const auto *dpg0_195 = buffer.data(dpg0 + 195);
    const auto *dpg0_197 = buffer.data(dpg0 + 197);
    const auto *dpg0_199 = buffer.data(dpg0 + 199);
    const auto *dpg0_200 = buffer.data(dpg0 + 200);

    const auto *dpg1_150 = buffer.data(dpg1 + 150);
    const auto *dpg1_151 = buffer.data(dpg1 + 151);
    const auto *dpg1_153 = buffer.data(dpg1 + 153);
    const auto *dpg1_155 = buffer.data(dpg1 + 155);
    const auto *dpg1_156 = buffer.data(dpg1 + 156);
    const auto *dpg1_158 = buffer.data(dpg1 + 158);
    const auto *dpg1_159 = buffer.data(dpg1 + 159);
    const auto *dpg1_160 = buffer.data(dpg1 + 160);
    const auto *dpg1_161 = buffer.data(dpg1 + 161);
    const auto *dpg1_162 = buffer.data(dpg1 + 162);
    const auto *dpg1_163 = buffer.data(dpg1 + 163);
    const auto *dpg1_164 = buffer.data(dpg1 + 164);
    const auto *dpg1_170 = buffer.data(dpg1 + 170);
    const auto *dpg1_173 = buffer.data(dpg1 + 173);
    const auto *dpg1_174 = buffer.data(dpg1 + 174);
    const auto *dpg1_177 = buffer.data(dpg1 + 177);
    const auto *dpg1_178 = buffer.data(dpg1 + 178);
    const auto *dpg1_179 = buffer.data(dpg1 + 179);
    const auto *dpg1_195 = buffer.data(dpg1 + 195);
    const auto *dpg1_197 = buffer.data(dpg1 + 197);
    const auto *dpg1_199 = buffer.data(dpg1 + 199);
    const auto *dpg1_200 = buffer.data(dpg1 + 200);

    const auto *dph_188 = buffer.data(dph + 188);
    const auto *dph_189 = buffer.data(dph + 189);
    const auto *dph_190 = buffer.data(dph + 190);
    const auto *dph_192 = buffer.data(dph + 192);
    const auto *dph_195 = buffer.data(dph + 195);
    const auto *dph_204 = buffer.data(dph + 204);
    const auto *dph_205 = buffer.data(dph + 205);
    const auto *dph_206 = buffer.data(dph + 206);
    const auto *dph_207 = buffer.data(dph + 207);
    const auto *dph_208 = buffer.data(dph + 208);
    const auto *dph_209 = buffer.data(dph + 209);
    const auto *dph_210 = buffer.data(dph + 210);
    const auto *dph_211 = buffer.data(dph + 211);
    const auto *dph_213 = buffer.data(dph + 213);
    const auto *dph_215 = buffer.data(dph + 215);
    const auto *dph_216 = buffer.data(dph + 216);
    const auto *dph_218 = buffer.data(dph + 218);
    const auto *dph_219 = buffer.data(dph + 219);
    const auto *dph_220 = buffer.data(dph + 220);
    const auto *dph_222 = buffer.data(dph + 222);
    const auto *dph_223 = buffer.data(dph + 223);
    const auto *dph_224 = buffer.data(dph + 224);
    const auto *dph_225 = buffer.data(dph + 225);
    const auto *dph_226 = buffer.data(dph + 226);
    const auto *dph_227 = buffer.data(dph + 227);
    const auto *dph_228 = buffer.data(dph + 228);
    const auto *dph_229 = buffer.data(dph + 229);
    const auto *dph_230 = buffer.data(dph + 230);
    const auto *dph_231 = buffer.data(dph + 231);
    const auto *dph_232 = buffer.data(dph + 232);
    const auto *dph_234 = buffer.data(dph + 234);
    const auto *dph_236 = buffer.data(dph + 236);
    const auto *dph_237 = buffer.data(dph + 237);
    const auto *dph_239 = buffer.data(dph + 239);
    const auto *dph_240 = buffer.data(dph + 240);
    const auto *dph_243 = buffer.data(dph + 243);
    const auto *dph_244 = buffer.data(dph + 244);
    const auto *dph_245 = buffer.data(dph + 245);
    const auto *dph_246 = buffer.data(dph + 246);
    const auto *dph_247 = buffer.data(dph + 247);
    const auto *dph_248 = buffer.data(dph + 248);
    const auto *dph_249 = buffer.data(dph + 249);
    const auto *dph_250 = buffer.data(dph + 250);
    const auto *dph_251 = buffer.data(dph + 251);
    const auto *dph_267 = buffer.data(dph + 267);
    const auto *dph_268 = buffer.data(dph + 268);
    const auto *dph_269 = buffer.data(dph + 269);
    const auto *dph_270 = buffer.data(dph + 270);
    const auto *dph_271 = buffer.data(dph + 271);
    const auto *dph_272 = buffer.data(dph + 272);
    const auto *dph_273 = buffer.data(dph + 273);
    const auto *dph_275 = buffer.data(dph + 275);
    const auto *dph_277 = buffer.data(dph + 277);
    const auto *dph_278 = buffer.data(dph + 278);

#pragma omp simd aligned(t_246, t_247, t_248, t_249, pa_x, pc_x, ppi0_246, ppi0_247, ppi0_248, \
                         ppi0_249, ppi1_246, ppi1_247, ppi1_248, \
                         ppi1_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_x[k] * ppi0_246[k]
                   - f_11 * pc_x[k] * ppi1_246[k];

        t_247[k] = pa_x[k] * ppi0_247[k]
                   - f_11 * pc_x[k] * ppi1_247[k];

        t_248[k] = pa_x[k] * ppi0_248[k]
                   - f_11 * pc_x[k] * ppi1_248[k];

        t_249[k] = pa_x[k] * ppi0_249[k]
                   - f_11 * pc_x[k] * ppi1_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pa_x, pb_x, pc_x, pc_y, ppi0_251, ppi1_251, \
                         dsi0_84, dsh_63, dsi1_84, dph_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_4 * pc_y[k] * dph_188[k];

        t_251[k] = pa_x[k] * ppi0_251[k]
                   - f_11 * pc_x[k] * ppi1_251[k];

        t_252[k] = pb_x[k] * dsi0_84[k]
                   + f_14 * dsh_63[k]
                   - f_11 * pc_x[k] * dsi1_84[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pb_x, pc_x, pc_z, dsi0_85, dsi0_87, \
                         dsh_64, dsh_66, dsi1_85, dsi1_87, dph_189, \
                         dph_190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = pb_x[k] * dsi0_85[k]
                   + f_17 * dsh_64[k]
                   - f_11 * pc_x[k] * dsi1_85[k];

        t_254[k] = f_4 * pc_z[k] * dph_189[k];

        t_255[k] = pb_x[k] * dsi0_87[k]
                   + f_13 * dsh_66[k]
                   - f_11 * pc_x[k] * dsi1_87[k];

        t_256[k] = f_4 * pc_z[k] * dph_190[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pb_x, pc_x, pc_z, dsi0_89, dsi0_90, dsh_68, \
                         dsh_69, dsi1_89, dsi1_90, dph_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = pb_x[k] * dsi0_89[k]
                   + f_13 * dsh_68[k]
                   - f_11 * pc_x[k] * dsi1_89[k];

        t_258[k] = pb_x[k] * dsi0_90[k]
                   + f_12 * dsh_69[k]
                   - f_11 * pc_x[k] * dsi1_90[k];

        t_259[k] = f_4 * pc_z[k] * dph_192[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, pb_x, pc_x, dsi0_92, dsi0_93, dsi0_94, dsh_71, \
                         dsh_72, dsh_73, dsi1_92, dsi1_93, dsi1_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pb_x[k] * dsi0_92[k]
                   + f_12 * dsh_71[k]
                   - f_11 * pc_x[k] * dsi1_92[k];

        t_261[k] = pb_x[k] * dsi0_93[k]
                   + f_12 * dsh_72[k]
                   - f_11 * pc_x[k] * dsi1_93[k];

        t_262[k] = pb_x[k] * dsi0_94[k]
                   + f_0 * dsh_73[k]
                   - f_11 * pc_x[k] * dsi1_94[k];
    }

#pragma omp simd aligned(t_263, t_264, t_265, pb_x, pc_x, pc_z, dsi0_96, dsi0_97, dsh_75, \
                         dsh_76, dsi1_96, dsi1_97, dph_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_263[k] = f_4 * pc_z[k] * dph_195[k];

        t_264[k] = pb_x[k] * dsi0_96[k]
                   + f_0 * dsh_75[k]
                   - f_11 * pc_x[k] * dsi1_96[k];

        t_265[k] = pb_x[k] * dsi0_97[k]
                   + f_0 * dsh_76[k]
                   - f_11 * pc_x[k] * dsi1_97[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, t_269, pb_x, pc_x, dsi0_98, dsh_77, dsh_78, \
                         dsh_79, dsh_80, dsi1_98, dph_204, dph_205, \
                         dph_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = pb_x[k] * dsi0_98[k]
                   + f_0 * dsh_77[k]
                   - f_11 * pc_x[k] * dsi1_98[k];

        t_267[k] = f_1 * dsh_78[k]
                   + f_4 * pc_x[k] * dph_204[k];

        t_268[k] = f_1 * dsh_79[k]
                   + f_4 * pc_x[k] * dph_205[k];

        t_269[k] = f_1 * dsh_80[k]
                   + f_4 * pc_x[k] * dph_206[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pb_x, pc_x, dsi0_105, dsh_81, dsh_82, \
                         dsh_83, dsi1_105, dph_207, dph_208, dph_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * dsh_81[k]
                   + f_4 * pc_x[k] * dph_207[k];

        t_271[k] = f_1 * dsh_82[k]
                   + f_4 * pc_x[k] * dph_208[k];

        t_272[k] = f_1 * dsh_83[k]
                   + f_4 * pc_x[k] * dph_209[k];

        t_273[k] = pb_x[k] * dsi0_105[k]
                   - f_11 * pc_x[k] * dsi1_105[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, pb_x, pc_x, pc_z, dsi0_107, dsi0_108, \
                         dsi0_109, dsi1_107, dsi1_108, dsi1_109, \
                         dph_204 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_4 * pc_z[k] * dph_204[k];

        t_275[k] = pb_x[k] * dsi0_107[k]
                   - f_11 * pc_x[k] * dsi1_107[k];

        t_276[k] = pb_x[k] * dsi0_108[k]
                   - f_11 * pc_x[k] * dsi1_108[k];

        t_277[k] = pb_x[k] * dsi0_109[k]
                   - f_11 * pc_x[k] * dsi1_109[k];
    }

#pragma omp simd aligned(t_278, t_279, t_280, pb_x, pc_x, pc_y, pph_83, dsi0_111, dsi1_111, \
                         dpg0_150, dpg1_150, dph_209, dph_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_278[k] = f_0 * pph_83[k]
                   + f_4 * pc_y[k] * dph_209[k];

        t_279[k] = pb_x[k] * dsi0_111[k]
                   - f_11 * pc_x[k] * dsi1_111[k];

        t_280[k] = f_2 * dpg0_150[k]
                   - f_3 * dpg1_150[k]
                   + f_4 * pc_x[k] * dph_210[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, t_284, pc_x, pc_z, dpg0_151, dpg0_153, dpg1_151, \
                         dpg1_153, dph_210, dph_211, dph_213 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_15 * dpg0_151[k]
                   - f_16 * dpg1_151[k]
                   + f_4 * pc_x[k] * dph_211[k];

        t_282[k] = f_4 * pc_z[k] * dph_210[k];

        t_283[k] = f_9 * dpg0_153[k]
                   - f_10 * dpg1_153[k]
                   + f_4 * pc_x[k] * dph_213[k];

        t_284[k] = f_4 * pc_z[k] * dph_211[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, t_288, pc_x, pc_z, dpg0_155, dpg0_156, dpg0_158, \
                         dpg1_155, dpg1_156, dpg1_158, dph_213, dph_215, dph_216, \
                         dph_218 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_9 * dpg0_155[k]
                   - f_10 * dpg1_155[k]
                   + f_4 * pc_x[k] * dph_215[k];

        t_286[k] = f_7 * dpg0_156[k]
                   - f_8 * dpg1_156[k]
                   + f_4 * pc_x[k] * dph_216[k];

        t_287[k] = f_4 * pc_z[k] * dph_213[k];

        t_288[k] = f_7 * dpg0_158[k]
                   - f_8 * dpg1_158[k]
                   + f_4 * pc_x[k] * dph_218[k];
    }

#pragma omp simd aligned(t_289, t_290, t_291, t_292, pc_x, pc_z, dpg0_159, dpg0_160, dpg0_162, \
                         dpg1_159, dpg1_160, dpg1_162, dph_216, dph_219, dph_220, \
                         dph_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_289[k] = f_7 * dpg0_159[k]
                   - f_8 * dpg1_159[k]
                   + f_4 * pc_x[k] * dph_219[k];

        t_290[k] = f_5 * dpg0_160[k]
                   - f_6 * dpg1_160[k]
                   + f_4 * pc_x[k] * dph_220[k];

        t_291[k] = f_4 * pc_z[k] * dph_216[k];

        t_292[k] = f_5 * dpg0_162[k]
                   - f_6 * dpg1_162[k]
                   + f_4 * pc_x[k] * dph_222[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, t_297, pc_x, dpg0_163, dpg0_164, \
                         dpg1_163, dpg1_164, dph_223, dph_224, dph_225, dph_226, \
                         dph_227 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_5 * dpg0_163[k]
                   - f_6 * dpg1_163[k]
                   + f_4 * pc_x[k] * dph_223[k];

        t_294[k] = f_5 * dpg0_164[k]
                   - f_6 * dpg1_164[k]
                   + f_4 * pc_x[k] * dph_224[k];

        t_295[k] = f_4 * pc_x[k] * dph_225[k];

        t_296[k] = f_4 * pc_x[k] * dph_226[k];

        t_297[k] = f_4 * pc_x[k] * dph_227[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, t_302, pc_x, pc_y, pc_z, pph_99, dsh_78, \
                         dpg0_160, dpg1_160, dph_225, dph_228, dph_229, \
                         dph_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_4 * pc_x[k] * dph_228[k];

        t_299[k] = f_4 * pc_x[k] * dph_229[k];

        t_300[k] = f_4 * pc_x[k] * dph_230[k];

        t_301[k] = f_0 * pph_99[k]
                   + f_1 * dsh_78[k]
                   + f_2 * dpg0_160[k]
                   - f_3 * dpg1_160[k]
                   + f_4 * pc_y[k] * dph_225[k];

        t_302[k] = f_4 * pc_z[k] * dph_225[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pc_z, dpg0_160, dpg0_161, dpg0_162, dpg1_160, \
                         dpg1_161, dpg1_162, dph_226, dph_227, \
                         dph_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_5 * dpg0_160[k]
                   - f_6 * dpg1_160[k]
                   + f_4 * pc_z[k] * dph_226[k];

        t_304[k] = f_7 * dpg0_161[k]
                   - f_8 * dpg1_161[k]
                   + f_4 * pc_z[k] * dph_227[k];

        t_305[k] = f_9 * dpg0_162[k]
                   - f_10 * dpg1_162[k]
                   + f_4 * pc_z[k] * dph_228[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pb_z, pc_y, pc_z, pph_104, dsi0_84, \
                         dsi0_85, dsh_83, dsi1_84, dsi1_85, dpg0_164, dpg1_164, \
                         dph_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_0 * pph_104[k]
                   + f_1 * dsh_83[k]
                   + f_4 * pc_y[k] * dph_230[k];

        t_307[k] = f_2 * dpg0_164[k]
                   - f_3 * dpg1_164[k]
                   + f_4 * pc_z[k] * dph_230[k];

        t_308[k] = pb_z[k] * dsi0_84[k]
                   - f_11 * pc_z[k] * dsi1_84[k];

        t_309[k] = pb_z[k] * dsi0_85[k]
                   - f_11 * pc_z[k] * dsi1_85[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pb_z, pc_x, pc_z, dsi0_87, dsh_63, \
                         dsh_64, dsi1_87, dpg0_170, dpg1_170, dph_231, dph_232, \
                         dph_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_1 * dsh_63[k]
                   + f_4 * pc_z[k] * dph_231[k];

        t_311[k] = pb_z[k] * dsi0_87[k]
                   - f_11 * pc_z[k] * dsi1_87[k];

        t_312[k] = f_1 * dsh_64[k]
                   + f_4 * pc_z[k] * dph_232[k];

        t_313[k] = f_9 * dpg0_170[k]
                   - f_10 * dpg1_170[k]
                   + f_4 * pc_x[k] * dph_236[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, pb_z, pc_x, pc_z, dsi0_90, dsh_66, dsi1_90, \
                         dpg0_173, dpg1_173, dph_234, dph_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pb_z[k] * dsi0_90[k]
                   - f_11 * pc_z[k] * dsi1_90[k];

        t_315[k] = f_1 * dsh_66[k]
                   + f_4 * pc_z[k] * dph_234[k];

        t_316[k] = f_7 * dpg0_173[k]
                   - f_8 * dpg1_173[k]
                   + f_4 * pc_x[k] * dph_239[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, pb_z, pc_x, pc_z, dsi0_94, dsh_69, dsi1_94, \
                         dpg0_174, dpg1_174, dph_237, dph_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = f_7 * dpg0_174[k]
                   - f_8 * dpg1_174[k]
                   + f_4 * pc_x[k] * dph_240[k];

        t_318[k] = pb_z[k] * dsi0_94[k]
                   - f_11 * pc_z[k] * dsi1_94[k];

        t_319[k] = f_1 * dsh_69[k]
                   + f_4 * pc_z[k] * dph_237[k];
    }

#pragma omp simd aligned(t_320, t_321, t_322, t_323, pc_x, dpg0_177, dpg0_178, dpg0_179, \
                         dpg1_177, dpg1_178, dpg1_179, dph_243, dph_244, dph_245, \
                         dph_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_320[k] = f_5 * dpg0_177[k]
                   - f_6 * dpg1_177[k]
                   + f_4 * pc_x[k] * dph_243[k];

        t_321[k] = f_5 * dpg0_178[k]
                   - f_6 * dpg1_178[k]
                   + f_4 * pc_x[k] * dph_244[k];

        t_322[k] = f_5 * dpg0_179[k]
                   - f_6 * dpg1_179[k]
                   + f_4 * pc_x[k] * dph_245[k];

        t_323[k] = f_4 * pc_x[k] * dph_246[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pb_z, pc_x, pc_z, dsi0_105, \
                         dsi1_105, dph_247, dph_248, dph_249, dph_250, \
                         dph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_4 * pc_x[k] * dph_247[k];

        t_325[k] = f_4 * pc_x[k] * dph_248[k];

        t_326[k] = f_4 * pc_x[k] * dph_249[k];

        t_327[k] = f_4 * pc_x[k] * dph_250[k];

        t_328[k] = f_4 * pc_x[k] * dph_251[k];

        t_329[k] = pb_z[k] * dsi0_105[k]
                   - f_11 * pc_z[k] * dsi1_105[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pb_z, pc_z, dsi0_107, dsi0_108, dsh_78, dsh_79, \
                         dsh_80, dsi1_107, dsi1_108, dph_246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_1 * dsh_78[k]
                   + f_4 * pc_z[k] * dph_246[k];

        t_331[k] = pb_z[k] * dsi0_107[k]
                   + f_0 * dsh_79[k]
                   - f_11 * pc_z[k] * dsi1_107[k];

        t_332[k] = pb_z[k] * dsi0_108[k]
                   + f_12 * dsh_80[k]
                   - f_11 * pc_z[k] * dsi1_108[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pb_z, pc_y, pc_z, pph_125, dsi0_109, dsi0_111, \
                         dsh_81, dsh_83, dsi1_109, dsi1_111, dph_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = pb_z[k] * dsi0_109[k]
                   + f_13 * dsh_81[k]
                   - f_11 * pc_z[k] * dsi1_109[k];

        t_334[k] = f_0 * pph_125[k]
                   + f_4 * pc_y[k] * dph_251[k];

        t_335[k] = pb_z[k] * dsi0_111[k]
                   + f_14 * dsh_83[k]
                   - f_11 * pc_z[k] * dsi1_111[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_y, pa_z, pc_y, pc_z, ppi0_85, ppi0_87, \
                         ppi0_168, ppi0_170, ppi1_85, ppi1_87, ppi1_168, \
                         ppi1_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_y[k] * ppi0_168[k]
                   - f_11 * pc_y[k] * ppi1_168[k];

        t_337[k] = pa_z[k] * ppi0_85[k]
                   - f_11 * pc_z[k] * ppi1_85[k];

        t_338[k] = pa_y[k] * ppi0_170[k]
                   - f_11 * pc_y[k] * ppi1_170[k];

        t_339[k] = pa_z[k] * ppi0_87[k]
                   - f_11 * pc_z[k] * ppi1_87[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_y, pa_z, pc_y, pc_z, ppi0_90, ppi0_172, \
                         ppi0_173, pph_128, ppi1_90, ppi1_172, \
                         ppi1_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = pa_y[k] * ppi0_172[k]
                   + f_1 * pph_128[k]
                   - f_11 * pc_y[k] * ppi1_172[k];

        t_341[k] = pa_y[k] * ppi0_173[k]
                   - f_11 * pc_y[k] * ppi1_173[k];

        t_342[k] = pa_z[k] * ppi0_90[k]
                   - f_11 * pc_z[k] * ppi1_90[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pa_y, pc_y, ppi0_175, ppi0_176, ppi0_177, \
                         pph_130, pph_131, ppi1_175, ppi1_176, \
                         ppi1_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = pa_y[k] * ppi0_175[k]
                   + f_0 * pph_130[k]
                   - f_11 * pc_y[k] * ppi1_175[k];

        t_344[k] = pa_y[k] * ppi0_176[k]
                   + f_1 * pph_131[k]
                   - f_11 * pc_y[k] * ppi1_176[k];

        t_345[k] = pa_y[k] * ppi0_177[k]
                   - f_11 * pc_y[k] * ppi1_177[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pa_y, pa_z, pc_y, pc_z, ppi0_94, ppi0_179, \
                         ppi0_180, pph_133, pph_134, ppi1_94, ppi1_179, \
                         ppi1_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = pa_z[k] * ppi0_94[k]
                   - f_11 * pc_z[k] * ppi1_94[k];

        t_347[k] = pa_y[k] * ppi0_179[k]
                   + f_12 * pph_133[k]
                   - f_11 * pc_y[k] * ppi1_179[k];

        t_348[k] = pa_y[k] * ppi0_180[k]
                   + f_0 * pph_134[k]
                   - f_11 * pc_y[k] * ppi1_180[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, t_352, pa_y, pc_x, pc_y, ppi0_181, ppi0_182, \
                         pph_135, ppi1_181, ppi1_182, dsh_99, dsh_100, dph_267, \
                         dph_268 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = pa_y[k] * ppi0_181[k]
                   + f_1 * pph_135[k]
                   - f_11 * pc_y[k] * ppi1_181[k];

        t_350[k] = pa_y[k] * ppi0_182[k]
                   - f_11 * pc_y[k] * ppi1_182[k];

        t_351[k] = f_1 * dsh_99[k]
                   + f_4 * pc_x[k] * dph_267[k];

        t_352[k] = f_1 * dsh_100[k]
                   + f_4 * pc_x[k] * dph_268[k];
    }

#pragma omp simd aligned(t_353, t_354, t_355, t_356, pc_x, dsh_101, dsh_102, dsh_103, dsh_104, \
                         dph_269, dph_270, dph_271, dph_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_353[k] = f_1 * dsh_101[k]
                   + f_4 * pc_x[k] * dph_269[k];

        t_354[k] = f_1 * dsh_102[k]
                   + f_4 * pc_x[k] * dph_270[k];

        t_355[k] = f_1 * dsh_103[k]
                   + f_4 * pc_x[k] * dph_271[k];

        t_356[k] = f_1 * dsh_104[k]
                   + f_4 * pc_x[k] * dph_272[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, t_360, pa_z, pb_x, pc_x, pc_z, ppi0_105, pph_78, \
                         ppi1_105, dsi0_135, dsi0_136, dsi1_135, dsi1_136, \
                         dph_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = pa_z[k] * ppi0_105[k]
                   - f_11 * pc_z[k] * ppi1_105[k];

        t_358[k] = f_1 * pph_78[k]
                   + f_4 * pc_z[k] * dph_267[k];

        t_359[k] = pb_x[k] * dsi0_135[k]
                   - f_11 * pc_x[k] * dsi1_135[k];

        t_360[k] = pb_x[k] * dsi0_136[k]
                   - f_11 * pc_x[k] * dsi1_136[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pa_y, pb_x, pc_x, pc_y, ppi0_195, pph_146, \
                         ppi1_195, dsi0_137, dsi1_137, dph_272 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = pb_x[k] * dsi0_137[k]
                   - f_11 * pc_x[k] * dsi1_137[k];

        t_362[k] = f_1 * pph_146[k]
                   + f_4 * pc_y[k] * dph_272[k];

        t_363[k] = pa_y[k] * ppi0_195[k]
                   - f_11 * pc_y[k] * ppi1_195[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pa_z, pc_x, pc_z, ppi0_113, ppi1_113, dpg0_195, \
                         dpg0_197, dpg1_195, dpg1_197, dph_273, \
                         dph_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = f_2 * dpg0_195[k]
                   - f_3 * dpg1_195[k]
                   + f_4 * pc_x[k] * dph_273[k];

        t_365[k] = pa_z[k] * ppi0_113[k]
                   - f_11 * pc_z[k] * ppi1_113[k];

        t_366[k] = f_15 * dpg0_197[k]
                   - f_16 * dpg1_197[k]
                   + f_4 * pc_x[k] * dph_275[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, pa_z, pc_x, pc_z, ppi0_115, ppi1_115, dpg0_199, \
                         dpg0_200, dpg1_199, dpg1_200, dph_277, \
                         dph_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = pa_z[k] * ppi0_115[k]
                   - f_11 * pc_z[k] * ppi1_115[k];

        t_368[k] = f_9 * dpg0_199[k]
                   - f_10 * dpg1_199[k]
                   + f_4 * pc_x[k] * dph_277[k];

        t_369[k] = f_9 * dpg0_200[k]
                   - f_10 * dpg1_200[k]
                   + f_4 * pc_x[k] * dph_278[k];
    }
}

static auto
compute_prim_dpi_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t ppi0, const size_t pph,
                                                          const size_t ppi1, const size_t dsi0,
                                                          const size_t dsh, const size_t dsi1,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t dph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_12 = 1.5 / q;
    const auto f_13 = 2.0 / q;
    const auto f_14 = 3.0 / q;
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);
    const auto f_17 = 2.5 / q;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *ppi0_118 = buffer.data(ppi0 + 118);
    const auto *ppi0_122 = buffer.data(ppi0 + 122);
    const auto *ppi0_133 = buffer.data(ppi0 + 133);
    const auto *ppi0_224 = buffer.data(ppi0 + 224);
    const auto *ppi0_226 = buffer.data(ppi0 + 226);
    const auto *ppi0_229 = buffer.data(ppi0 + 229);
    const auto *ppi0_233 = buffer.data(ppi0 + 233);
    const auto *ppi0_238 = buffer.data(ppi0 + 238);
    const auto *ppi0_247 = buffer.data(ppi0 + 247);
    const auto *ppi0_248 = buffer.data(ppi0 + 248);
    const auto *ppi0_249 = buffer.data(ppi0 + 249);
    const auto *ppi0_251 = buffer.data(ppi0 + 251);

    const auto *pph_99 = buffer.data(pph + 99);
    const auto *pph_104 = buffer.data(pph + 104);
    const auto *pph_120 = buffer.data(pph + 120);
    const auto *pph_164 = buffer.data(pph + 164);
    const auto *pph_165 = buffer.data(pph + 165);
    const auto *pph_166 = buffer.data(pph + 166);
    const auto *pph_167 = buffer.data(pph + 167);
    const auto *pph_183 = buffer.data(pph + 183);
    const auto *pph_185 = buffer.data(pph + 185);
    const auto *pph_186 = buffer.data(pph + 186);
    const auto *pph_187 = buffer.data(pph + 187);
    const auto *pph_188 = buffer.data(pph + 188);

    const auto *ppi1_118 = buffer.data(ppi1 + 118);
    const auto *ppi1_122 = buffer.data(ppi1 + 122);
    const auto *ppi1_133 = buffer.data(ppi1 + 133);
    const auto *ppi1_224 = buffer.data(ppi1 + 224);
    const auto *ppi1_226 = buffer.data(ppi1 + 226);
    const auto *ppi1_229 = buffer.data(ppi1 + 229);
    const auto *ppi1_233 = buffer.data(ppi1 + 233);
    const auto *ppi1_238 = buffer.data(ppi1 + 238);
    const auto *ppi1_247 = buffer.data(ppi1 + 247);
    const auto *ppi1_248 = buffer.data(ppi1 + 248);
    const auto *ppi1_249 = buffer.data(ppi1 + 249);
    const auto *ppi1_251 = buffer.data(ppi1 + 251);

    const auto *dsi0_140 = buffer.data(dsi0 + 140);
    const auto *dsi0_142 = buffer.data(dsi0 + 142);
    const auto *dsi0_143 = buffer.data(dsi0 + 143);
    const auto *dsi0_145 = buffer.data(dsi0 + 145);
    const auto *dsi0_146 = buffer.data(dsi0 + 146);
    const auto *dsi0_147 = buffer.data(dsi0 + 147);
    const auto *dsi0_149 = buffer.data(dsi0 + 149);
    const auto *dsi0_150 = buffer.data(dsi0 + 150);
    const auto *dsi0_151 = buffer.data(dsi0 + 151);
    const auto *dsi0_152 = buffer.data(dsi0 + 152);
    const auto *dsi0_154 = buffer.data(dsi0 + 154);
    const auto *dsi0_161 = buffer.data(dsi0 + 161);
    const auto *dsi0_162 = buffer.data(dsi0 + 162);
    const auto *dsi0_163 = buffer.data(dsi0 + 163);
    const auto *dsi0_164 = buffer.data(dsi0 + 164);
    const auto *dsi0_165 = buffer.data(dsi0 + 165);
    const auto *dsi0_167 = buffer.data(dsi0 + 167);

    const auto *dsh_99 = buffer.data(dsh + 99);
    const auto *dsh_101 = buffer.data(dsh + 101);
    const auto *dsh_102 = buffer.data(dsh + 102);
    const auto *dsh_103 = buffer.data(dsh + 103);
    const auto *dsh_104 = buffer.data(dsh + 104);
    const auto *dsh_105 = buffer.data(dsh + 105);
    const auto *dsh_107 = buffer.data(dsh + 107);
    const auto *dsh_108 = buffer.data(dsh + 108);
    const auto *dsh_110 = buffer.data(dsh + 110);
    const auto *dsh_111 = buffer.data(dsh + 111);
    const auto *dsh_112 = buffer.data(dsh + 112);
    const auto *dsh_114 = buffer.data(dsh + 114);
    const auto *dsh_115 = buffer.data(dsh + 115);
    const auto *dsh_116 = buffer.data(dsh + 116);
    const auto *dsh_117 = buffer.data(dsh + 117);
    const auto *dsh_119 = buffer.data(dsh + 119);
    const auto *dsh_120 = buffer.data(dsh + 120);
    const auto *dsh_121 = buffer.data(dsh + 121);
    const auto *dsh_122 = buffer.data(dsh + 122);
    const auto *dsh_123 = buffer.data(dsh + 123);
    const auto *dsh_124 = buffer.data(dsh + 124);
    const auto *dsh_125 = buffer.data(dsh + 125);

    const auto *dsi1_140 = buffer.data(dsi1 + 140);
    const auto *dsi1_142 = buffer.data(dsi1 + 142);
    const auto *dsi1_143 = buffer.data(dsi1 + 143);
    const auto *dsi1_145 = buffer.data(dsi1 + 145);
    const auto *dsi1_146 = buffer.data(dsi1 + 146);
    const auto *dsi1_147 = buffer.data(dsi1 + 147);
    const auto *dsi1_149 = buffer.data(dsi1 + 149);
    const auto *dsi1_150 = buffer.data(dsi1 + 150);
    const auto *dsi1_151 = buffer.data(dsi1 + 151);
    const auto *dsi1_152 = buffer.data(dsi1 + 152);
    const auto *dsi1_154 = buffer.data(dsi1 + 154);
    const auto *dsi1_161 = buffer.data(dsi1 + 161);
    const auto *dsi1_162 = buffer.data(dsi1 + 162);
    const auto *dsi1_163 = buffer.data(dsi1 + 163);
    const auto *dsi1_164 = buffer.data(dsi1 + 164);
    const auto *dsi1_165 = buffer.data(dsi1 + 165);
    const auto *dsi1_167 = buffer.data(dsi1 + 167);

    const auto *dpg0_202 = buffer.data(dpg0 + 202);
    const auto *dpg0_203 = buffer.data(dpg0 + 203);
    const auto *dpg0_204 = buffer.data(dpg0 + 204);
    const auto *dpg0_206 = buffer.data(dpg0 + 206);
    const auto *dpg0_207 = buffer.data(dpg0 + 207);
    const auto *dpg0_208 = buffer.data(dpg0 + 208);
    const auto *dpg0_209 = buffer.data(dpg0 + 209);
    const auto *dpg0_211 = buffer.data(dpg0 + 211);
    const auto *dpg0_213 = buffer.data(dpg0 + 213);
    const auto *dpg0_214 = buffer.data(dpg0 + 214);
    const auto *dpg0_216 = buffer.data(dpg0 + 216);
    const auto *dpg0_217 = buffer.data(dpg0 + 217);
    const auto *dpg0_218 = buffer.data(dpg0 + 218);
    const auto *dpg0_220 = buffer.data(dpg0 + 220);
    const auto *dpg0_221 = buffer.data(dpg0 + 221);
    const auto *dpg0_222 = buffer.data(dpg0 + 222);
    const auto *dpg0_223 = buffer.data(dpg0 + 223);
    const auto *dpg0_243 = buffer.data(dpg0 + 243);
    const auto *dpg0_246 = buffer.data(dpg0 + 246);
    const auto *dpg0_247 = buffer.data(dpg0 + 247);
    const auto *dpg0_250 = buffer.data(dpg0 + 250);
    const auto *dpg0_251 = buffer.data(dpg0 + 251);
    const auto *dpg0_252 = buffer.data(dpg0 + 252);
    const auto *dpg0_255 = buffer.data(dpg0 + 255);
    const auto *dpg0_257 = buffer.data(dpg0 + 257);
    const auto *dpg0_258 = buffer.data(dpg0 + 258);
    const auto *dpg0_260 = buffer.data(dpg0 + 260);
    const auto *dpg0_261 = buffer.data(dpg0 + 261);
    const auto *dpg0_262 = buffer.data(dpg0 + 262);
    const auto *dpg0_264 = buffer.data(dpg0 + 264);
    const auto *dpg0_265 = buffer.data(dpg0 + 265);
    const auto *dpg0_266 = buffer.data(dpg0 + 266);
    const auto *dpg0_267 = buffer.data(dpg0 + 267);

    const auto *dpg1_202 = buffer.data(dpg1 + 202);
    const auto *dpg1_203 = buffer.data(dpg1 + 203);
    const auto *dpg1_204 = buffer.data(dpg1 + 204);
    const auto *dpg1_206 = buffer.data(dpg1 + 206);
    const auto *dpg1_207 = buffer.data(dpg1 + 207);
    const auto *dpg1_208 = buffer.data(dpg1 + 208);
    const auto *dpg1_209 = buffer.data(dpg1 + 209);
    const auto *dpg1_211 = buffer.data(dpg1 + 211);
    const auto *dpg1_213 = buffer.data(dpg1 + 213);
    const auto *dpg1_214 = buffer.data(dpg1 + 214);
    const auto *dpg1_216 = buffer.data(dpg1 + 216);
    const auto *dpg1_217 = buffer.data(dpg1 + 217);
    const auto *dpg1_218 = buffer.data(dpg1 + 218);
    const auto *dpg1_220 = buffer.data(dpg1 + 220);
    const auto *dpg1_221 = buffer.data(dpg1 + 221);
    const auto *dpg1_222 = buffer.data(dpg1 + 222);
    const auto *dpg1_223 = buffer.data(dpg1 + 223);
    const auto *dpg1_243 = buffer.data(dpg1 + 243);
    const auto *dpg1_246 = buffer.data(dpg1 + 246);
    const auto *dpg1_247 = buffer.data(dpg1 + 247);
    const auto *dpg1_250 = buffer.data(dpg1 + 250);
    const auto *dpg1_251 = buffer.data(dpg1 + 251);
    const auto *dpg1_252 = buffer.data(dpg1 + 252);
    const auto *dpg1_255 = buffer.data(dpg1 + 255);
    const auto *dpg1_257 = buffer.data(dpg1 + 257);
    const auto *dpg1_258 = buffer.data(dpg1 + 258);
    const auto *dpg1_260 = buffer.data(dpg1 + 260);
    const auto *dpg1_261 = buffer.data(dpg1 + 261);
    const auto *dpg1_262 = buffer.data(dpg1 + 262);
    const auto *dpg1_264 = buffer.data(dpg1 + 264);
    const auto *dpg1_265 = buffer.data(dpg1 + 265);
    const auto *dpg1_266 = buffer.data(dpg1 + 266);
    const auto *dpg1_267 = buffer.data(dpg1 + 267);

    const auto *dph_280 = buffer.data(dph + 280);
    const auto *dph_281 = buffer.data(dph + 281);
    const auto *dph_282 = buffer.data(dph + 282);
    const auto *dph_284 = buffer.data(dph + 284);
    const auto *dph_285 = buffer.data(dph + 285);
    const auto *dph_286 = buffer.data(dph + 286);
    const auto *dph_287 = buffer.data(dph + 287);
    const auto *dph_288 = buffer.data(dph + 288);
    const auto *dph_289 = buffer.data(dph + 289);
    const auto *dph_290 = buffer.data(dph + 290);
    const auto *dph_291 = buffer.data(dph + 291);
    const auto *dph_292 = buffer.data(dph + 292);
    const auto *dph_293 = buffer.data(dph + 293);
    const auto *dph_295 = buffer.data(dph + 295);
    const auto *dph_297 = buffer.data(dph + 297);
    const auto *dph_298 = buffer.data(dph + 298);
    const auto *dph_300 = buffer.data(dph + 300);
    const auto *dph_301 = buffer.data(dph + 301);
    const auto *dph_302 = buffer.data(dph + 302);
    const auto *dph_304 = buffer.data(dph + 304);
    const auto *dph_305 = buffer.data(dph + 305);
    const auto *dph_306 = buffer.data(dph + 306);
    const auto *dph_307 = buffer.data(dph + 307);
    const auto *dph_309 = buffer.data(dph + 309);
    const auto *dph_310 = buffer.data(dph + 310);
    const auto *dph_311 = buffer.data(dph + 311);
    const auto *dph_312 = buffer.data(dph + 312);
    const auto *dph_313 = buffer.data(dph + 313);
    const auto *dph_314 = buffer.data(dph + 314);
    const auto *dph_315 = buffer.data(dph + 315);
    const auto *dph_317 = buffer.data(dph + 317);
    const auto *dph_320 = buffer.data(dph + 320);
    const auto *dph_324 = buffer.data(dph + 324);
    const auto *dph_330 = buffer.data(dph + 330);
    const auto *dph_331 = buffer.data(dph + 331);
    const auto *dph_332 = buffer.data(dph + 332);
    const auto *dph_333 = buffer.data(dph + 333);
    const auto *dph_334 = buffer.data(dph + 334);
    const auto *dph_335 = buffer.data(dph + 335);
    const auto *dph_336 = buffer.data(dph + 336);
    const auto *dph_338 = buffer.data(dph + 338);
    const auto *dph_339 = buffer.data(dph + 339);
    const auto *dph_341 = buffer.data(dph + 341);
    const auto *dph_342 = buffer.data(dph + 342);
    const auto *dph_343 = buffer.data(dph + 343);
    const auto *dph_345 = buffer.data(dph + 345);
    const auto *dph_346 = buffer.data(dph + 346);
    const auto *dph_347 = buffer.data(dph + 347);
    const auto *dph_348 = buffer.data(dph + 348);
    const auto *dph_351 = buffer.data(dph + 351);
    const auto *dph_352 = buffer.data(dph + 352);
    const auto *dph_353 = buffer.data(dph + 353);
    const auto *dph_354 = buffer.data(dph + 354);
    const auto *dph_355 = buffer.data(dph + 355);
    const auto *dph_356 = buffer.data(dph + 356);
    const auto *dph_357 = buffer.data(dph + 357);
    const auto *dph_359 = buffer.data(dph + 359);
    const auto *dph_360 = buffer.data(dph + 360);
    const auto *dph_362 = buffer.data(dph + 362);
    const auto *dph_363 = buffer.data(dph + 363);
    const auto *dph_364 = buffer.data(dph + 364);
    const auto *dph_366 = buffer.data(dph + 366);
    const auto *dph_367 = buffer.data(dph + 367);
    const auto *dph_368 = buffer.data(dph + 368);
    const auto *dph_369 = buffer.data(dph + 369);

#pragma omp simd aligned(t_370, t_371, t_372, pa_z, pc_x, pc_z, ppi0_118, ppi1_118, dpg0_202, \
                         dpg0_203, dpg1_202, dpg1_203, dph_280, \
                         dph_281 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pa_z[k] * ppi0_118[k]
                   - f_11 * pc_z[k] * ppi1_118[k];

        t_371[k] = f_7 * dpg0_202[k]
                   - f_8 * dpg1_202[k]
                   + f_4 * pc_x[k] * dph_280[k];

        t_372[k] = f_7 * dpg0_203[k]
                   - f_8 * dpg1_203[k]
                   + f_4 * pc_x[k] * dph_281[k];
    }

#pragma omp simd aligned(t_373, t_374, t_375, pa_z, pc_x, pc_z, ppi0_122, ppi1_122, dpg0_204, \
                         dpg0_206, dpg1_204, dpg1_206, dph_282, \
                         dph_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_7 * dpg0_204[k]
                   - f_8 * dpg1_204[k]
                   + f_4 * pc_x[k] * dph_282[k];

        t_374[k] = pa_z[k] * ppi0_122[k]
                   - f_11 * pc_z[k] * ppi1_122[k];

        t_375[k] = f_5 * dpg0_206[k]
                   - f_6 * dpg1_206[k]
                   + f_4 * pc_x[k] * dph_284[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, t_379, pc_x, dpg0_207, dpg0_208, dpg0_209, \
                         dpg1_207, dpg1_208, dpg1_209, dph_285, dph_286, dph_287, \
                         dph_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = f_5 * dpg0_207[k]
                   - f_6 * dpg1_207[k]
                   + f_4 * pc_x[k] * dph_285[k];

        t_377[k] = f_5 * dpg0_208[k]
                   - f_6 * dpg1_208[k]
                   + f_4 * pc_x[k] * dph_286[k];

        t_378[k] = f_5 * dpg0_209[k]
                   - f_6 * dpg1_209[k]
                   + f_4 * pc_x[k] * dph_287[k];

        t_379[k] = f_4 * pc_x[k] * dph_288[k];
    }

#pragma omp simd aligned(t_380, t_381, t_382, t_383, t_384, t_385, pa_z, pc_x, pc_z, ppi0_133, \
                         ppi1_133, dph_289, dph_290, dph_291, dph_292, \
                         dph_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_380[k] = f_4 * pc_x[k] * dph_289[k];

        t_381[k] = f_4 * pc_x[k] * dph_290[k];

        t_382[k] = f_4 * pc_x[k] * dph_291[k];

        t_383[k] = f_4 * pc_x[k] * dph_292[k];

        t_384[k] = f_4 * pc_x[k] * dph_293[k];

        t_385[k] = pa_z[k] * ppi0_133[k]
                   - f_11 * pc_z[k] * ppi1_133[k];
    }

#pragma omp simd aligned(t_386, t_387, pc_y, pc_z, pph_99, pph_164, dsh_101, dpg0_207, \
                         dpg1_207, dph_288, dph_290 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_1 * pph_99[k]
                   + f_4 * pc_z[k] * dph_288[k];

        t_387[k] = f_1 * pph_164[k]
                   + f_1 * dsh_101[k]
                   + f_9 * dpg0_207[k]
                   - f_10 * dpg1_207[k]
                   + f_4 * pc_y[k] * dph_290[k];
    }

#pragma omp simd aligned(t_388, t_389, pc_y, pph_165, pph_166, dsh_102, dsh_103, dpg0_208, \
                         dpg0_209, dpg1_208, dpg1_209, dph_291, \
                         dph_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_1 * pph_165[k]
                   + f_1 * dsh_102[k]
                   + f_7 * dpg0_208[k]
                   - f_8 * dpg1_208[k]
                   + f_4 * pc_y[k] * dph_291[k];

        t_389[k] = f_1 * pph_166[k]
                   + f_1 * dsh_103[k]
                   + f_5 * dpg0_209[k]
                   - f_6 * dpg1_209[k]
                   + f_4 * pc_y[k] * dph_292[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, pa_y, pc_y, pc_z, ppi0_224, pph_104, pph_167, \
                         ppi1_224, dsh_104, dpg0_209, dpg1_209, \
                         dph_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = f_1 * pph_167[k]
                   + f_1 * dsh_104[k]
                   + f_4 * pc_y[k] * dph_293[k];

        t_391[k] = f_1 * pph_104[k]
                   + f_2 * dpg0_209[k]
                   - f_3 * dpg1_209[k]
                   + f_4 * pc_z[k] * dph_293[k];

        t_392[k] = pa_y[k] * ppi0_224[k]
                   - f_11 * pc_y[k] * ppi1_224[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, pa_y, pc_x, pc_y, ppi0_226, ppi1_226, dpg0_211, \
                         dpg0_213, dpg1_211, dpg1_213, dph_295, \
                         dph_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_15 * dpg0_211[k]
                   - f_16 * dpg1_211[k]
                   + f_4 * pc_x[k] * dph_295[k];

        t_394[k] = pa_y[k] * ppi0_226[k]
                   - f_11 * pc_y[k] * ppi1_226[k];

        t_395[k] = f_9 * dpg0_213[k]
                   - f_10 * dpg1_213[k]
                   + f_4 * pc_x[k] * dph_297[k];
    }

#pragma omp simd aligned(t_396, t_397, t_398, pa_y, pc_x, pc_y, ppi0_229, ppi1_229, dpg0_214, \
                         dpg0_216, dpg1_214, dpg1_216, dph_298, \
                         dph_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_396[k] = f_9 * dpg0_214[k]
                   - f_10 * dpg1_214[k]
                   + f_4 * pc_x[k] * dph_298[k];

        t_397[k] = pa_y[k] * ppi0_229[k]
                   - f_11 * pc_y[k] * ppi1_229[k];

        t_398[k] = f_7 * dpg0_216[k]
                   - f_8 * dpg1_216[k]
                   + f_4 * pc_x[k] * dph_300[k];
    }

#pragma omp simd aligned(t_399, t_400, t_401, pa_y, pc_x, pc_y, ppi0_233, ppi1_233, dpg0_217, \
                         dpg0_218, dpg1_217, dpg1_218, dph_301, \
                         dph_302 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_399[k] = f_7 * dpg0_217[k]
                   - f_8 * dpg1_217[k]
                   + f_4 * pc_x[k] * dph_301[k];

        t_400[k] = f_7 * dpg0_218[k]
                   - f_8 * dpg1_218[k]
                   + f_4 * pc_x[k] * dph_302[k];

        t_401[k] = pa_y[k] * ppi0_233[k]
                   - f_11 * pc_y[k] * ppi1_233[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pc_x, dpg0_220, dpg0_221, dpg0_222, dpg1_220, \
                         dpg1_221, dpg1_222, dph_304, dph_305, \
                         dph_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = f_5 * dpg0_220[k]
                   - f_6 * dpg1_220[k]
                   + f_4 * pc_x[k] * dph_304[k];

        t_403[k] = f_5 * dpg0_221[k]
                   - f_6 * dpg1_221[k]
                   + f_4 * pc_x[k] * dph_305[k];

        t_404[k] = f_5 * dpg0_222[k]
                   - f_6 * dpg1_222[k]
                   + f_4 * pc_x[k] * dph_306[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, t_409, pa_y, pc_x, pc_y, ppi0_238, \
                         ppi1_238, dpg0_223, dpg1_223, dph_307, dph_309, dph_310, \
                         dph_311 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_5 * dpg0_223[k]
                   - f_6 * dpg1_223[k]
                   + f_4 * pc_x[k] * dph_307[k];

        t_406[k] = pa_y[k] * ppi0_238[k]
                   - f_11 * pc_y[k] * ppi1_238[k];

        t_407[k] = f_4 * pc_x[k] * dph_309[k];

        t_408[k] = f_4 * pc_x[k] * dph_310[k];

        t_409[k] = f_4 * pc_x[k] * dph_311[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pc_x, pc_y, pph_183, dpg0_220, dpg1_220, \
                         dph_309, dph_312, dph_313, dph_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = f_4 * pc_x[k] * dph_312[k];

        t_411[k] = f_4 * pc_x[k] * dph_313[k];

        t_412[k] = f_4 * pc_x[k] * dph_314[k];

        t_413[k] = f_1 * pph_183[k]
                   + f_2 * dpg0_220[k]
                   - f_3 * dpg1_220[k]
                   + f_4 * pc_y[k] * dph_309[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pa_y, pc_y, pc_z, ppi0_247, ppi0_248, pph_120, \
                         pph_185, pph_186, ppi1_247, ppi1_248, dsh_99, \
                         dph_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_1 * pph_120[k]
                   + f_1 * dsh_99[k]
                   + f_4 * pc_z[k] * dph_309[k];

        t_415[k] = pa_y[k] * ppi0_247[k]
                   + f_13 * pph_185[k]
                   - f_11 * pc_y[k] * ppi1_247[k];

        t_416[k] = pa_y[k] * ppi0_248[k]
                   + f_12 * pph_186[k]
                   - f_11 * pc_y[k] * ppi1_248[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pa_y, pc_y, ppi0_249, ppi0_251, pph_187, \
                         pph_188, ppi1_249, ppi1_251, dph_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = pa_y[k] * ppi0_249[k]
                   + f_0 * pph_187[k]
                   - f_11 * pc_y[k] * ppi1_249[k];

        t_418[k] = f_1 * pph_188[k]
                   + f_4 * pc_y[k] * dph_314[k];

        t_419[k] = pa_y[k] * ppi0_251[k]
                   - f_11 * pc_y[k] * ppi1_251[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, pb_x, pc_x, pc_y, dsi0_140, dsi0_142, dsh_105, \
                         dsh_107, dsi1_140, dsi1_142, dph_315 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_x[k] * dsi0_140[k]
                   + f_14 * dsh_105[k]
                   - f_11 * pc_x[k] * dsi1_140[k];

        t_421[k] = f_4 * pc_y[k] * dph_315[k];

        t_422[k] = pb_x[k] * dsi0_142[k]
                   + f_17 * dsh_107[k]
                   - f_11 * pc_x[k] * dsi1_142[k];
    }

#pragma omp simd aligned(t_423, t_424, t_425, pb_x, pc_x, pc_y, dsi0_143, dsi0_145, dsh_108, \
                         dsh_110, dsi1_143, dsi1_145, dph_317 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_423[k] = pb_x[k] * dsi0_143[k]
                   + f_13 * dsh_108[k]
                   - f_11 * pc_x[k] * dsi1_143[k];

        t_424[k] = f_4 * pc_y[k] * dph_317[k];

        t_425[k] = pb_x[k] * dsi0_145[k]
                   + f_13 * dsh_110[k]
                   - f_11 * pc_x[k] * dsi1_145[k];
    }

#pragma omp simd aligned(t_426, t_427, t_428, pb_x, pc_x, pc_y, dsi0_146, dsi0_147, dsh_111, \
                         dsh_112, dsi1_146, dsi1_147, dph_320 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_426[k] = pb_x[k] * dsi0_146[k]
                   + f_12 * dsh_111[k]
                   - f_11 * pc_x[k] * dsi1_146[k];

        t_427[k] = pb_x[k] * dsi0_147[k]
                   + f_12 * dsh_112[k]
                   - f_11 * pc_x[k] * dsi1_147[k];

        t_428[k] = f_4 * pc_y[k] * dph_320[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pb_x, pc_x, dsi0_149, dsi0_150, dsi0_151, \
                         dsh_114, dsh_115, dsh_116, dsi1_149, dsi1_150, \
                         dsi1_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = pb_x[k] * dsi0_149[k]
                   + f_12 * dsh_114[k]
                   - f_11 * pc_x[k] * dsi1_149[k];

        t_430[k] = pb_x[k] * dsi0_150[k]
                   + f_0 * dsh_115[k]
                   - f_11 * pc_x[k] * dsi1_150[k];

        t_431[k] = pb_x[k] * dsi0_151[k]
                   + f_0 * dsh_116[k]
                   - f_11 * pc_x[k] * dsi1_151[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, t_435, pb_x, pc_x, pc_y, dsi0_152, dsi0_154, \
                         dsh_117, dsh_119, dsh_120, dsi1_152, dsi1_154, dph_324, \
                         dph_330 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = pb_x[k] * dsi0_152[k]
                   + f_0 * dsh_117[k]
                   - f_11 * pc_x[k] * dsi1_152[k];

        t_433[k] = f_4 * pc_y[k] * dph_324[k];

        t_434[k] = pb_x[k] * dsi0_154[k]
                   + f_0 * dsh_119[k]
                   - f_11 * pc_x[k] * dsi1_154[k];

        t_435[k] = f_1 * dsh_120[k]
                   + f_4 * pc_x[k] * dph_330[k];
    }

#pragma omp simd aligned(t_436, t_437, t_438, t_439, t_440, pc_x, dsh_121, dsh_122, dsh_123, \
                         dsh_124, dsh_125, dph_331, dph_332, dph_333, dph_334, \
                         dph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_436[k] = f_1 * dsh_121[k]
                   + f_4 * pc_x[k] * dph_331[k];

        t_437[k] = f_1 * dsh_122[k]
                   + f_4 * pc_x[k] * dph_332[k];

        t_438[k] = f_1 * dsh_123[k]
                   + f_4 * pc_x[k] * dph_333[k];

        t_439[k] = f_1 * dsh_124[k]
                   + f_4 * pc_x[k] * dph_334[k];

        t_440[k] = f_1 * dsh_125[k]
                   + f_4 * pc_x[k] * dph_335[k];
    }

#pragma omp simd aligned(t_441, t_442, t_443, t_444, pb_x, pc_x, dsi0_161, dsi0_162, dsi0_163, \
                         dsi0_164, dsi1_161, dsi1_162, dsi1_163, \
                         dsi1_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_441[k] = pb_x[k] * dsi0_161[k]
                   - f_11 * pc_x[k] * dsi1_161[k];

        t_442[k] = pb_x[k] * dsi0_162[k]
                   - f_11 * pc_x[k] * dsi1_162[k];

        t_443[k] = pb_x[k] * dsi0_163[k]
                   - f_11 * pc_x[k] * dsi1_163[k];

        t_444[k] = pb_x[k] * dsi0_164[k]
                   - f_11 * pc_x[k] * dsi1_164[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pb_x, pb_y, pc_x, pc_y, dsi0_140, \
                         dsi0_165, dsi0_167, dsi1_140, dsi1_165, dsi1_167, \
                         dph_335 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = pb_x[k] * dsi0_165[k]
                   - f_11 * pc_x[k] * dsi1_165[k];

        t_446[k] = f_4 * pc_y[k] * dph_335[k];

        t_447[k] = pb_x[k] * dsi0_167[k]
                   - f_11 * pc_x[k] * dsi1_167[k];

        t_448[k] = pb_y[k] * dsi0_140[k]
                   - f_11 * pc_y[k] * dsi1_140[k];
    }

#pragma omp simd aligned(t_449, t_450, t_451, t_452, pb_y, pc_x, pc_y, dsi0_142, dsh_105, \
                         dsh_107, dsi1_142, dpg0_243, dpg1_243, dph_336, dph_338, \
                         dph_339 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_1 * dsh_105[k]
                   + f_4 * pc_y[k] * dph_336[k];

        t_450[k] = pb_y[k] * dsi0_142[k]
                   - f_11 * pc_y[k] * dsi1_142[k];

        t_451[k] = f_9 * dpg0_243[k]
                   - f_10 * dpg1_243[k]
                   + f_4 * pc_x[k] * dph_339[k];

        t_452[k] = f_1 * dsh_107[k]
                   + f_4 * pc_y[k] * dph_338[k];
    }

#pragma omp simd aligned(t_453, t_454, t_455, pb_y, pc_x, pc_y, dsi0_145, dsi1_145, dpg0_246, \
                         dpg0_247, dpg1_246, dpg1_247, dph_342, \
                         dph_343 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_453[k] = pb_y[k] * dsi0_145[k]
                   - f_11 * pc_y[k] * dsi1_145[k];

        t_454[k] = f_7 * dpg0_246[k]
                   - f_8 * dpg1_246[k]
                   + f_4 * pc_x[k] * dph_342[k];

        t_455[k] = f_7 * dpg0_247[k]
                   - f_8 * dpg1_247[k]
                   + f_4 * pc_x[k] * dph_343[k];
    }

#pragma omp simd aligned(t_456, t_457, t_458, pb_y, pc_x, pc_y, dsi0_149, dsh_110, dsi1_149, \
                         dpg0_250, dpg1_250, dph_341, dph_346 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_456[k] = f_1 * dsh_110[k]
                   + f_4 * pc_y[k] * dph_341[k];

        t_457[k] = pb_y[k] * dsi0_149[k]
                   - f_11 * pc_y[k] * dsi1_149[k];

        t_458[k] = f_5 * dpg0_250[k]
                   - f_6 * dpg1_250[k]
                   + f_4 * pc_x[k] * dph_346[k];
    }

#pragma omp simd aligned(t_459, t_460, t_461, pc_x, pc_y, dsh_114, dpg0_251, dpg0_252, \
                         dpg1_251, dpg1_252, dph_345, dph_347, \
                         dph_348 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_459[k] = f_5 * dpg0_251[k]
                   - f_6 * dpg1_251[k]
                   + f_4 * pc_x[k] * dph_347[k];

        t_460[k] = f_5 * dpg0_252[k]
                   - f_6 * dpg1_252[k]
                   + f_4 * pc_x[k] * dph_348[k];

        t_461[k] = f_1 * dsh_114[k]
                   + f_4 * pc_y[k] * dph_345[k];
    }

#pragma omp simd aligned(t_462, t_463, t_464, t_465, t_466, t_467, pb_y, pc_x, pc_y, dsi0_154, \
                         dsi1_154, dph_351, dph_352, dph_353, dph_354, \
                         dph_355 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_462[k] = pb_y[k] * dsi0_154[k]
                   - f_11 * pc_y[k] * dsi1_154[k];

        t_463[k] = f_4 * pc_x[k] * dph_351[k];

        t_464[k] = f_4 * pc_x[k] * dph_352[k];

        t_465[k] = f_4 * pc_x[k] * dph_353[k];

        t_466[k] = f_4 * pc_x[k] * dph_354[k];

        t_467[k] = f_4 * pc_x[k] * dph_355[k];
    }

#pragma omp simd aligned(t_468, t_469, t_470, pb_y, pc_x, pc_y, dsi0_161, dsi0_162, dsh_120, \
                         dsh_121, dsi1_161, dsi1_162, dph_356 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_468[k] = f_4 * pc_x[k] * dph_356[k];

        t_469[k] = pb_y[k] * dsi0_161[k]
                   + f_14 * dsh_120[k]
                   - f_11 * pc_y[k] * dsi1_161[k];

        t_470[k] = pb_y[k] * dsi0_162[k]
                   + f_17 * dsh_121[k]
                   - f_11 * pc_y[k] * dsi1_162[k];
    }

#pragma omp simd aligned(t_471, t_472, t_473, pb_y, pc_y, dsi0_163, dsi0_164, dsi0_165, \
                         dsh_122, dsh_123, dsh_124, dsi1_163, dsi1_164, \
                         dsi1_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_471[k] = pb_y[k] * dsi0_163[k]
                   + f_13 * dsh_122[k]
                   - f_11 * pc_y[k] * dsi1_163[k];

        t_472[k] = pb_y[k] * dsi0_164[k]
                   + f_12 * dsh_123[k]
                   - f_11 * pc_y[k] * dsi1_164[k];

        t_473[k] = pb_y[k] * dsi0_165[k]
                   + f_0 * dsh_124[k]
                   - f_11 * pc_y[k] * dsi1_165[k];
    }

#pragma omp simd aligned(t_474, t_475, t_476, t_477, pb_y, pc_x, pc_y, dsi0_167, dsh_125, \
                         dsi1_167, dpg0_255, dpg1_255, dph_356, \
                         dph_357 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_474[k] = f_1 * dsh_125[k]
                   + f_4 * pc_y[k] * dph_356[k];

        t_475[k] = pb_y[k] * dsi0_167[k]
                   - f_11 * pc_y[k] * dsi1_167[k];

        t_476[k] = f_2 * dpg0_255[k]
                   - f_3 * dpg1_255[k]
                   + f_4 * pc_x[k] * dph_357[k];

        t_477[k] = f_4 * pc_y[k] * dph_357[k];
    }

#pragma omp simd aligned(t_478, t_479, t_480, t_481, pc_x, pc_y, dpg0_257, dpg0_258, dpg0_260, \
                         dpg1_257, dpg1_258, dpg1_260, dph_359, dph_360, \
                         dph_362 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_478[k] = f_15 * dpg0_257[k]
                   - f_16 * dpg1_257[k]
                   + f_4 * pc_x[k] * dph_359[k];

        t_479[k] = f_9 * dpg0_258[k]
                   - f_10 * dpg1_258[k]
                   + f_4 * pc_x[k] * dph_360[k];

        t_480[k] = f_4 * pc_y[k] * dph_359[k];

        t_481[k] = f_9 * dpg0_260[k]
                   - f_10 * dpg1_260[k]
                   + f_4 * pc_x[k] * dph_362[k];
    }

#pragma omp simd aligned(t_482, t_483, t_484, t_485, pc_x, pc_y, dpg0_261, dpg0_262, dpg0_264, \
                         dpg1_261, dpg1_262, dpg1_264, dph_362, dph_363, dph_364, \
                         dph_366 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_482[k] = f_7 * dpg0_261[k]
                   - f_8 * dpg1_261[k]
                   + f_4 * pc_x[k] * dph_363[k];

        t_483[k] = f_7 * dpg0_262[k]
                   - f_8 * dpg1_262[k]
                   + f_4 * pc_x[k] * dph_364[k];

        t_484[k] = f_4 * pc_y[k] * dph_362[k];

        t_485[k] = f_7 * dpg0_264[k]
                   - f_8 * dpg1_264[k]
                   + f_4 * pc_x[k] * dph_366[k];
    }

#pragma omp simd aligned(t_486, t_487, t_488, t_489, pc_x, pc_y, dpg0_265, dpg0_266, dpg0_267, \
                         dpg1_265, dpg1_266, dpg1_267, dph_366, dph_367, dph_368, \
                         dph_369 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_486[k] = f_5 * dpg0_265[k]
                   - f_6 * dpg1_265[k]
                   + f_4 * pc_x[k] * dph_367[k];

        t_487[k] = f_5 * dpg0_266[k]
                   - f_6 * dpg1_266[k]
                   + f_4 * pc_x[k] * dph_368[k];

        t_488[k] = f_5 * dpg0_267[k]
                   - f_6 * dpg1_267[k]
                   + f_4 * pc_x[k] * dph_369[k];

        t_489[k] = f_4 * pc_y[k] * dph_366[k];
    }
}

static auto
compute_prim_dpi_three_center_electron_repulsion_0_piece4(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pc,
                                                          const size_t pph, const size_t dsh,
                                                          const size_t dpg0, const size_t dpg1,
                                                          const size_t dph, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / q;
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
    const auto f_15 = 2.0 / gamma;
    const auto f_16 = 2.0 * p / (gamma * q);

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

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *pph_188 = buffer.data(pph + 188);

    const auto *dsh_125 = buffer.data(dsh + 125);

    const auto *dpg0_265 = buffer.data(dpg0 + 265);
    const auto *dpg0_266 = buffer.data(dpg0 + 266);
    const auto *dpg0_267 = buffer.data(dpg0 + 267);
    const auto *dpg0_268 = buffer.data(dpg0 + 268);
    const auto *dpg0_269 = buffer.data(dpg0 + 269);

    const auto *dpg1_265 = buffer.data(dpg1 + 265);
    const auto *dpg1_266 = buffer.data(dpg1 + 266);
    const auto *dpg1_267 = buffer.data(dpg1 + 267);
    const auto *dpg1_268 = buffer.data(dpg1 + 268);
    const auto *dpg1_269 = buffer.data(dpg1 + 269);

    const auto *dph_371 = buffer.data(dph + 371);
    const auto *dph_372 = buffer.data(dph + 372);
    const auto *dph_373 = buffer.data(dph + 373);
    const auto *dph_374 = buffer.data(dph + 374);
    const auto *dph_375 = buffer.data(dph + 375);
    const auto *dph_376 = buffer.data(dph + 376);
    const auto *dph_377 = buffer.data(dph + 377);

#pragma omp simd aligned(t_490, t_491, t_492, t_493, t_494, t_495, pc_x, dpg0_269, dpg1_269, \
                         dph_371, dph_372, dph_373, dph_374, dph_375, \
                         dph_376 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_490[k] = f_5 * dpg0_269[k]
                   - f_6 * dpg1_269[k]
                   + f_4 * pc_x[k] * dph_371[k];

        t_491[k] = f_4 * pc_x[k] * dph_372[k];

        t_492[k] = f_4 * pc_x[k] * dph_373[k];

        t_493[k] = f_4 * pc_x[k] * dph_374[k];

        t_494[k] = f_4 * pc_x[k] * dph_375[k];

        t_495[k] = f_4 * pc_x[k] * dph_376[k];
    }

#pragma omp simd aligned(t_496, t_497, t_498, t_499, pc_x, pc_y, dpg0_265, dpg0_266, dpg0_267, \
                         dpg1_265, dpg1_266, dpg1_267, dph_372, dph_373, dph_374, \
                         dph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_496[k] = f_4 * pc_x[k] * dph_377[k];

        t_497[k] = f_2 * dpg0_265[k]
                   - f_3 * dpg1_265[k]
                   + f_4 * pc_y[k] * dph_372[k];

        t_498[k] = f_15 * dpg0_266[k]
                   - f_16 * dpg1_266[k]
                   + f_4 * pc_y[k] * dph_373[k];

        t_499[k] = f_9 * dpg0_267[k]
                   - f_10 * dpg1_267[k]
                   + f_4 * pc_y[k] * dph_374[k];
    }

#pragma omp simd aligned(t_500, t_501, t_502, t_503, pc_y, pc_z, pph_188, dsh_125, dpg0_268, \
                         dpg0_269, dpg1_268, dpg1_269, dph_375, dph_376, \
                         dph_377 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_500[k] = f_7 * dpg0_268[k]
                   - f_8 * dpg1_268[k]
                   + f_4 * pc_y[k] * dph_375[k];

        t_501[k] = f_5 * dpg0_269[k]
                   - f_6 * dpg1_269[k]
                   + f_4 * pc_y[k] * dph_376[k];

        t_502[k] = f_4 * pc_y[k] * dph_377[k];

        t_503[k] = f_0 * pph_188[k]
                   + f_1 * dsh_125[k]
                   + f_2 * dpg0_269[k]
                   - f_3 * dpg1_269[k]
                   + f_4 * pc_z[k] * dph_377[k];
    }
}

auto
compute_prim_dpi_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t ppi0,
                                                   const size_t pph, const size_t ppi1,
                                                   const size_t dsi0, const size_t dsh,
                                                   const size_t dsi1, const size_t dpg0,
                                                   const size_t dpg1, const size_t dph,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_dpi_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, ppi0,
                                                              pph, ppi1, dsi0, dsh, dsi1, dpg0,
                                                              dpg1, dph, ncols, gamma, p, q);

    compute_prim_dpi_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, ppi0,
                                                              pph, ppi1, dsi0, dsh, dsi1, dpg0,
                                                              dpg1, dph, ncols, gamma, p, q);

    compute_prim_dpi_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, ppi0,
                                                              pph, ppi1, dsi0, dsh, dsi1, dpg0,
                                                              dpg1, dph, ncols, gamma, p, q);

    compute_prim_dpi_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, ppi0,
                                                              pph, ppi1, dsi0, dsh, dsi1, dpg0,
                                                              dpg1, dph, ncols, gamma, p, q);

    compute_prim_dpi_three_center_electron_repulsion_0_piece4(buffer, target, pc, pph, dsh,
                                                              dpg0, dpg1, dph, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
