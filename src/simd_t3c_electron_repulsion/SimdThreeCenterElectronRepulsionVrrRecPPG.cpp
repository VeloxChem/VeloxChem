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


#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_prim_ppg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t spg0,
                                                   const size_t spf, const size_t spg1,
                                                   const size_t psg0, const size_t psf,
                                                   const size_t psg1, const size_t ppd0,
                                                   const size_t ppd1, const size_t ppf,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / gamma;
    const auto f_2 = 1.5 * p / (gamma * q);
    const auto f_3 = p / q;
    const auto f_4 = 0.5 / gamma;
    const auto f_5 = 0.5 * p / (gamma * q);
    const auto f_6 = gamma / q;
    const auto f_7 = 1.0 / q;
    const auto f_8 = 1.0 / gamma;
    const auto f_9 = p / (gamma * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *spg0_0 = buffer.data(spg0 + 0);
    const auto *spg0_1 = buffer.data(spg0 + 1);
    const auto *spg0_2 = buffer.data(spg0 + 2);
    const auto *spg0_3 = buffer.data(spg0 + 3);
    const auto *spg0_5 = buffer.data(spg0 + 5);
    const auto *spg0_10 = buffer.data(spg0 + 10);
    const auto *spg0_14 = buffer.data(spg0 + 14);
    const auto *spg0_15 = buffer.data(spg0 + 15);
    const auto *spg0_18 = buffer.data(spg0 + 18);
    const auto *spg0_25 = buffer.data(spg0 + 25);
    const auto *spg0_26 = buffer.data(spg0 + 26);
    const auto *spg0_27 = buffer.data(spg0 + 27);
    const auto *spg0_29 = buffer.data(spg0 + 29);
    const auto *spg0_30 = buffer.data(spg0 + 30);
    const auto *spg0_35 = buffer.data(spg0 + 35);
    const auto *spg0_40 = buffer.data(spg0 + 40);
    const auto *spg0_42 = buffer.data(spg0 + 42);
    const auto *spg0_44 = buffer.data(spg0 + 44);

    const auto *spf_0 = buffer.data(spf + 0);
    const auto *spf_1 = buffer.data(spf + 1);
    const auto *spf_2 = buffer.data(spf + 2);
    const auto *spf_6 = buffer.data(spf + 6);
    const auto *spf_9 = buffer.data(spf + 9);
    const auto *spf_13 = buffer.data(spf + 13);
    const auto *spf_16 = buffer.data(spf + 16);
    const auto *spf_17 = buffer.data(spf + 17);
    const auto *spf_19 = buffer.data(spf + 19);
    const auto *spf_25 = buffer.data(spf + 25);
    const auto *spf_26 = buffer.data(spf + 26);
    const auto *spf_28 = buffer.data(spf + 28);
    const auto *spf_29 = buffer.data(spf + 29);

    const auto *spg1_0 = buffer.data(spg1 + 0);
    const auto *spg1_1 = buffer.data(spg1 + 1);
    const auto *spg1_2 = buffer.data(spg1 + 2);
    const auto *spg1_3 = buffer.data(spg1 + 3);
    const auto *spg1_5 = buffer.data(spg1 + 5);
    const auto *spg1_10 = buffer.data(spg1 + 10);
    const auto *spg1_14 = buffer.data(spg1 + 14);
    const auto *spg1_15 = buffer.data(spg1 + 15);
    const auto *spg1_18 = buffer.data(spg1 + 18);
    const auto *spg1_25 = buffer.data(spg1 + 25);
    const auto *spg1_26 = buffer.data(spg1 + 26);
    const auto *spg1_27 = buffer.data(spg1 + 27);
    const auto *spg1_29 = buffer.data(spg1 + 29);
    const auto *spg1_30 = buffer.data(spg1 + 30);
    const auto *spg1_35 = buffer.data(spg1 + 35);
    const auto *spg1_40 = buffer.data(spg1 + 40);
    const auto *spg1_42 = buffer.data(spg1 + 42);
    const auto *spg1_44 = buffer.data(spg1 + 44);

    const auto *psg0_0 = buffer.data(psg0 + 0);
    const auto *psg0_3 = buffer.data(psg0 + 3);
    const auto *psg0_5 = buffer.data(psg0 + 5);
    const auto *psg0_16 = buffer.data(psg0 + 16);
    const auto *psg0_18 = buffer.data(psg0 + 18);
    const auto *psg0_25 = buffer.data(psg0 + 25);
    const auto *psg0_27 = buffer.data(psg0 + 27);
    const auto *psg0_32 = buffer.data(psg0 + 32);
    const auto *psg0_35 = buffer.data(psg0 + 35);
    const auto *psg0_41 = buffer.data(psg0 + 41);
    const auto *psg0_42 = buffer.data(psg0 + 42);
    const auto *psg0_44 = buffer.data(psg0 + 44);

    const auto *psf_0 = buffer.data(psf + 0);
    const auto *psf_2 = buffer.data(psf + 2);
    const auto *psf_3 = buffer.data(psf + 3);
    const auto *psf_5 = buffer.data(psf + 5);
    const auto *psf_6 = buffer.data(psf + 6);
    const auto *psf_9 = buffer.data(psf + 9);
    const auto *psf_10 = buffer.data(psf + 10);
    const auto *psf_11 = buffer.data(psf + 11);
    const auto *psf_16 = buffer.data(psf + 16);
    const auto *psf_17 = buffer.data(psf + 17);
    const auto *psf_18 = buffer.data(psf + 18);
    const auto *psf_19 = buffer.data(psf + 19);
    const auto *psf_20 = buffer.data(psf + 20);
    const auto *psf_22 = buffer.data(psf + 22);
    const auto *psf_26 = buffer.data(psf + 26);
    const auto *psf_27 = buffer.data(psf + 27);
    const auto *psf_28 = buffer.data(psf + 28);
    const auto *psf_29 = buffer.data(psf + 29);

    const auto *psg1_0 = buffer.data(psg1 + 0);
    const auto *psg1_3 = buffer.data(psg1 + 3);
    const auto *psg1_5 = buffer.data(psg1 + 5);
    const auto *psg1_16 = buffer.data(psg1 + 16);
    const auto *psg1_18 = buffer.data(psg1 + 18);
    const auto *psg1_25 = buffer.data(psg1 + 25);
    const auto *psg1_27 = buffer.data(psg1 + 27);
    const auto *psg1_32 = buffer.data(psg1 + 32);
    const auto *psg1_35 = buffer.data(psg1 + 35);
    const auto *psg1_41 = buffer.data(psg1 + 41);
    const auto *psg1_42 = buffer.data(psg1 + 42);
    const auto *psg1_44 = buffer.data(psg1 + 44);

    const auto *ppd0_0 = buffer.data(ppd0 + 0);
    const auto *ppd0_3 = buffer.data(ppd0 + 3);
    const auto *ppd0_5 = buffer.data(ppd0 + 5);
    const auto *ppd0_24 = buffer.data(ppd0 + 24);
    const auto *ppd0_25 = buffer.data(ppd0 + 25);
    const auto *ppd0_27 = buffer.data(ppd0 + 27);
    const auto *ppd0_29 = buffer.data(ppd0 + 29);
    const auto *ppd0_48 = buffer.data(ppd0 + 48);
    const auto *ppd0_50 = buffer.data(ppd0 + 50);
    const auto *ppd0_51 = buffer.data(ppd0 + 51);
    const auto *ppd0_52 = buffer.data(ppd0 + 52);
    const auto *ppd0_53 = buffer.data(ppd0 + 53);

    const auto *ppd1_0 = buffer.data(ppd1 + 0);
    const auto *ppd1_3 = buffer.data(ppd1 + 3);
    const auto *ppd1_5 = buffer.data(ppd1 + 5);
    const auto *ppd1_24 = buffer.data(ppd1 + 24);
    const auto *ppd1_25 = buffer.data(ppd1 + 25);
    const auto *ppd1_27 = buffer.data(ppd1 + 27);
    const auto *ppd1_29 = buffer.data(ppd1 + 29);
    const auto *ppd1_48 = buffer.data(ppd1 + 48);
    const auto *ppd1_50 = buffer.data(ppd1 + 50);
    const auto *ppd1_51 = buffer.data(ppd1 + 51);
    const auto *ppd1_52 = buffer.data(ppd1 + 52);
    const auto *ppd1_53 = buffer.data(ppd1 + 53);

    const auto *ppf_0 = buffer.data(ppf + 0);
    const auto *ppf_1 = buffer.data(ppf + 1);
    const auto *ppf_2 = buffer.data(ppf + 2);
    const auto *ppf_3 = buffer.data(ppf + 3);
    const auto *ppf_5 = buffer.data(ppf + 5);
    const auto *ppf_6 = buffer.data(ppf + 6);
    const auto *ppf_8 = buffer.data(ppf + 8);
    const auto *ppf_9 = buffer.data(ppf + 9);
    const auto *ppf_10 = buffer.data(ppf + 10);
    const auto *ppf_12 = buffer.data(ppf + 12);
    const auto *ppf_13 = buffer.data(ppf + 13);
    const auto *ppf_15 = buffer.data(ppf + 15);
    const auto *ppf_16 = buffer.data(ppf + 16);
    const auto *ppf_19 = buffer.data(ppf + 19);
    const auto *ppf_20 = buffer.data(ppf + 20);
    const auto *ppf_22 = buffer.data(ppf + 22);
    const auto *ppf_23 = buffer.data(ppf + 23);
    const auto *ppf_25 = buffer.data(ppf + 25);
    const auto *ppf_26 = buffer.data(ppf + 26);
    const auto *ppf_29 = buffer.data(ppf + 29);
    const auto *ppf_30 = buffer.data(ppf + 30);
    const auto *ppf_31 = buffer.data(ppf + 31);
    const auto *ppf_36 = buffer.data(ppf + 36);
    const auto *ppf_37 = buffer.data(ppf + 37);
    const auto *ppf_38 = buffer.data(ppf + 38);
    const auto *ppf_39 = buffer.data(ppf + 39);
    const auto *ppf_40 = buffer.data(ppf + 40);
    const auto *ppf_41 = buffer.data(ppf + 41);
    const auto *ppf_43 = buffer.data(ppf + 43);
    const auto *ppf_45 = buffer.data(ppf + 45);
    const auto *ppf_46 = buffer.data(ppf + 46);
    const auto *ppf_47 = buffer.data(ppf + 47);
    const auto *ppf_48 = buffer.data(ppf + 48);
    const auto *ppf_49 = buffer.data(ppf + 49);
    const auto *ppf_50 = buffer.data(ppf + 50);
    const auto *ppf_51 = buffer.data(ppf + 51);
    const auto *ppf_56 = buffer.data(ppf + 56);
    const auto *ppf_57 = buffer.data(ppf + 57);
    const auto *ppf_58 = buffer.data(ppf + 58);
    const auto *ppf_59 = buffer.data(ppf + 59);
    const auto *ppf_60 = buffer.data(ppf + 60);
    const auto *ppf_62 = buffer.data(ppf + 62);
    const auto *ppf_66 = buffer.data(ppf + 66);
    const auto *ppf_67 = buffer.data(ppf + 67);
    const auto *ppf_68 = buffer.data(ppf + 68);
    const auto *ppf_69 = buffer.data(ppf + 69);
    const auto *ppf_70 = buffer.data(ppf + 70);
    const auto *ppf_72 = buffer.data(ppf + 72);
    const auto *ppf_76 = buffer.data(ppf + 76);
    const auto *ppf_77 = buffer.data(ppf + 77);
    const auto *ppf_78 = buffer.data(ppf + 78);
    const auto *ppf_79 = buffer.data(ppf + 79);
    const auto *ppf_80 = buffer.data(ppf + 80);
    const auto *ppf_82 = buffer.data(ppf + 82);
    const auto *ppf_83 = buffer.data(ppf + 83);
    const auto *ppf_85 = buffer.data(ppf + 85);
    const auto *ppf_86 = buffer.data(ppf + 86);
    const auto *ppf_87 = buffer.data(ppf + 87);
    const auto *ppf_88 = buffer.data(ppf + 88);
    const auto *ppf_89 = buffer.data(ppf + 89);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, spf_0, psf_0, ppd0_0, \
                         ppd1_0, ppf_0, ppf_1, ppf_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * spf_0[k]
                 + f_0 * psf_0[k]
                 + f_1 * ppd0_0[k]
                 - f_2 * ppd1_0[k]
                 + f_3 * pc_x[k] * ppf_0[k];

        t_1[k] = f_3 * pc_y[k] * ppf_0[k];

        t_2[k] = f_3 * pc_z[k] * ppf_0[k];

        t_3[k] = f_4 * ppd0_0[k]
                 - f_5 * ppd1_0[k]
                 + f_3 * pc_y[k] * ppf_1[k];

        t_4[k] = f_3 * pc_y[k] * ppf_2[k];

        t_5[k] = f_4 * ppd0_0[k]
                 - f_5 * ppd1_0[k]
                 + f_3 * pc_z[k] * ppf_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, spf_6, spf_9, psf_6, psf_9, \
                         ppf_3, ppf_5, ppf_6, ppf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * spf_6[k]
                 + f_0 * psf_6[k]
                 + f_3 * pc_x[k] * ppf_6[k];

        t_7[k] = f_3 * pc_z[k] * ppf_3[k];

        t_8[k] = f_3 * pc_y[k] * ppf_5[k];

        t_9[k] = f_0 * spf_9[k]
                 + f_0 * psf_9[k]
                 + f_3 * pc_x[k] * ppf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, ppd0_3, ppd0_5, ppd1_3, \
                         ppd1_5, ppf_6, ppf_8, ppf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_1 * ppd0_3[k]
                  - f_2 * ppd1_3[k]
                  + f_3 * pc_y[k] * ppf_6[k];

        t_11[k] = f_3 * pc_z[k] * ppf_6[k];

        t_12[k] = f_4 * ppd0_5[k]
                  - f_5 * ppd1_5[k]
                  + f_3 * pc_y[k] * ppf_8[k];

        t_13[k] = f_3 * pc_y[k] * ppf_9[k];

        t_14[k] = f_1 * ppd0_5[k]
                  - f_2 * ppd1_5[k]
                  + f_3 * pc_z[k] * ppf_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pa_x, pb_y, pc_x, pc_y, pc_z, spg0_18, \
                         spf_13, spg1_18, psg0_0, psf_0, psg1_0, \
                         ppf_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * psg0_0[k]
                  - f_6 * pc_y[k] * psg1_0[k];

        t_16[k] = f_0 * psf_0[k]
                  + f_3 * pc_y[k] * ppf_10[k];

        t_17[k] = f_3 * pc_z[k] * ppf_10[k];

        t_18[k] = pa_x[k] * spg0_18[k]
                  + f_7 * spf_13[k]
                  - f_6 * pc_x[k] * spg1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, spf_16, psg0_5, \
                         psf_2, psg1_5, ppf_12, ppf_13, ppf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * psf_2[k]
                  + f_3 * pc_y[k] * ppf_12[k];

        t_20[k] = pb_y[k] * psg0_5[k]
                  - f_6 * pc_y[k] * psg1_5[k];

        t_21[k] = f_0 * spf_16[k]
                  + f_3 * pc_x[k] * ppf_16[k];

        t_22[k] = f_3 * pc_z[k] * ppf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_x, pc_x, pc_y, pc_z, spg0_25, spf_19, \
                         spg1_25, psf_5, ppf_15, ppf_16, ppf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * psf_5[k]
                  + f_3 * pc_y[k] * ppf_15[k];

        t_24[k] = f_0 * spf_19[k]
                  + f_3 * pc_x[k] * ppf_19[k];

        t_25[k] = pa_x[k] * spg0_25[k]
                  - f_6 * pc_x[k] * spg1_25[k];

        t_26[k] = f_3 * pc_z[k] * ppf_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pc_x, pc_y, spg0_27, spg0_29, spg1_27, \
                         spg1_29, psf_9, ppf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = pa_x[k] * spg0_27[k]
                  - f_6 * pc_x[k] * spg1_27[k];

        t_28[k] = f_0 * psf_9[k]
                  + f_3 * pc_y[k] * ppf_19[k];

        t_29[k] = pa_x[k] * spg0_29[k]
                  - f_6 * pc_x[k] * spg1_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, psg0_0, psg0_3, \
                         psf_0, psg1_0, psg1_3, ppf_20, ppf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * psg0_0[k]
                  - f_6 * pc_z[k] * psg1_0[k];

        t_31[k] = f_3 * pc_y[k] * ppf_20[k];

        t_32[k] = f_0 * psf_0[k]
                  + f_3 * pc_z[k] * ppf_20[k];

        t_33[k] = pb_z[k] * psg0_3[k]
                  - f_6 * pc_z[k] * psg1_3[k];

        t_34[k] = f_3 * pc_y[k] * ppf_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pa_x, pc_x, pc_y, pc_z, spg0_35, spf_25, \
                         spf_26, spg1_35, psf_3, ppf_23, ppf_25, \
                         ppf_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pa_x[k] * spg0_35[k]
                  + f_7 * spf_25[k]
                  - f_6 * pc_x[k] * spg1_35[k];

        t_36[k] = f_0 * spf_26[k]
                  + f_3 * pc_x[k] * ppf_26[k];

        t_37[k] = f_0 * psf_3[k]
                  + f_3 * pc_z[k] * ppf_23[k];

        t_38[k] = f_3 * pc_y[k] * ppf_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, t_42, pa_x, pc_x, pc_z, spg0_40, spg0_42, spf_29, \
                         spg1_40, spg1_42, psf_6, ppf_26, ppf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * spf_29[k]
                  + f_3 * pc_x[k] * ppf_29[k];

        t_40[k] = pa_x[k] * spg0_40[k]
                  - f_6 * pc_x[k] * spg1_40[k];

        t_41[k] = f_0 * psf_6[k]
                  + f_3 * pc_z[k] * ppf_26[k];

        t_42[k] = pa_x[k] * spg0_42[k]
                  - f_6 * pc_x[k] * spg1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, t_46, pa_x, pa_y, pc_x, pc_y, spg0_0, spg0_1, \
                         spg0_44, spf_0, spg1_0, spg1_1, spg1_44, \
                         ppf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_3 * pc_y[k] * ppf_29[k];

        t_44[k] = pa_x[k] * spg0_44[k]
                  - f_6 * pc_x[k] * spg1_44[k];

        t_45[k] = pa_y[k] * spg0_0[k]
                  - f_6 * pc_y[k] * spg1_0[k];

        t_46[k] = pa_y[k] * spg0_1[k]
                  + f_0 * spf_0[k]
                  - f_6 * pc_y[k] * spg1_1[k];
    }

#pragma omp simd aligned(t_47, t_48, t_49, t_50, pa_y, pc_y, pc_z, spg0_3, spg0_5, spf_1, \
                         spg1_3, spg1_5, ppf_30, ppf_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_3 * pc_z[k] * ppf_30[k];

        t_48[k] = pa_y[k] * spg0_3[k]
                  + f_7 * spf_1[k]
                  - f_6 * pc_y[k] * spg1_3[k];

        t_49[k] = f_3 * pc_z[k] * ppf_31[k];

        t_50[k] = pa_y[k] * spg0_5[k]
                  - f_6 * pc_y[k] * spg1_5[k];
    }

#pragma omp simd aligned(t_51, t_52, t_53, t_54, pc_x, psf_16, psf_17, psf_18, psf_19, ppf_36, \
                         ppf_37, ppf_38, ppf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_0 * psf_16[k]
                  + f_3 * pc_x[k] * ppf_36[k];

        t_52[k] = f_0 * psf_17[k]
                  + f_3 * pc_x[k] * ppf_37[k];

        t_53[k] = f_0 * psf_18[k]
                  + f_3 * pc_x[k] * ppf_38[k];

        t_54[k] = f_0 * psf_19[k]
                  + f_3 * pc_x[k] * ppf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, t_58, pb_x, pc_x, pc_y, pc_z, spf_9, psg0_25, \
                         psg0_27, psg1_25, psg1_27, ppf_36, ppf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = pb_x[k] * psg0_25[k]
                  - f_6 * pc_x[k] * psg1_25[k];

        t_56[k] = f_3 * pc_z[k] * ppf_36[k];

        t_57[k] = pb_x[k] * psg0_27[k]
                  - f_6 * pc_x[k] * psg1_27[k];

        t_58[k] = f_0 * spf_9[k]
                  + f_3 * pc_y[k] * ppf_39[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, t_62, pa_y, pc_x, pc_y, pc_z, spg0_14, spg1_14, \
                         ppd0_24, ppd0_25, ppd1_24, ppd1_25, ppf_40, \
                         ppf_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = pa_y[k] * spg0_14[k]
                  - f_6 * pc_y[k] * spg1_14[k];

        t_60[k] = f_1 * ppd0_24[k]
                  - f_2 * ppd1_24[k]
                  + f_3 * pc_x[k] * ppf_40[k];

        t_61[k] = f_8 * ppd0_25[k]
                  - f_9 * ppd1_25[k]
                  + f_3 * pc_x[k] * ppf_41[k];

        t_62[k] = f_3 * pc_z[k] * ppf_40[k];
    }

#pragma omp simd aligned(t_63, t_64, t_65, t_66, t_67, pc_x, pc_z, ppd0_27, ppd0_29, ppd1_27, \
                         ppd1_29, ppf_41, ppf_43, ppf_45, ppf_46, \
                         ppf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_4 * ppd0_27[k]
                  - f_5 * ppd1_27[k]
                  + f_3 * pc_x[k] * ppf_43[k];

        t_64[k] = f_3 * pc_z[k] * ppf_41[k];

        t_65[k] = f_4 * ppd0_29[k]
                  - f_5 * ppd1_29[k]
                  + f_3 * pc_x[k] * ppf_45[k];

        t_66[k] = f_3 * pc_x[k] * ppf_46[k];

        t_67[k] = f_3 * pc_x[k] * ppf_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, t_71, t_72, pc_x, pc_y, pc_z, spf_16, psf_16, \
                         ppd0_27, ppd1_27, ppf_46, ppf_47, ppf_48, \
                         ppf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * pc_x[k] * ppf_48[k];

        t_69[k] = f_3 * pc_x[k] * ppf_49[k];

        t_70[k] = f_0 * spf_16[k]
                  + f_0 * psf_16[k]
                  + f_1 * ppd0_27[k]
                  - f_2 * ppd1_27[k]
                  + f_3 * pc_y[k] * ppf_46[k];

        t_71[k] = f_3 * pc_z[k] * ppf_46[k];

        t_72[k] = f_4 * ppd0_27[k]
                  - f_5 * ppd1_27[k]
                  + f_3 * pc_z[k] * ppf_47[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pa_y, pc_y, pc_z, spg0_30, spf_19, spg1_30, psf_19, \
                         ppd0_29, ppd1_29, ppf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_0 * spf_19[k]
                  + f_0 * psf_19[k]
                  + f_3 * pc_y[k] * ppf_49[k];

        t_74[k] = f_1 * ppd0_29[k]
                  - f_2 * ppd1_29[k]
                  + f_3 * pc_z[k] * ppf_49[k];

        t_75[k] = pa_y[k] * spg0_30[k]
                  - f_6 * pc_y[k] * spg1_30[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, t_79, pb_z, pc_z, psg0_16, psg0_18, psf_10, psf_11, \
                         psg1_16, psg1_18, ppf_50, ppf_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = pb_z[k] * psg0_16[k]
                  - f_6 * pc_z[k] * psg1_16[k];

        t_77[k] = f_0 * psf_10[k]
                  + f_3 * pc_z[k] * ppf_50[k];

        t_78[k] = pb_z[k] * psg0_18[k]
                  - f_6 * pc_z[k] * psg1_18[k];

        t_79[k] = f_0 * psf_11[k]
                  + f_3 * pc_z[k] * ppf_51[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, t_83, t_84, pa_y, pc_x, pc_y, spg0_35, spg1_35, \
                         ppf_56, ppf_57, ppf_58, ppf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * spg0_35[k]
                  - f_6 * pc_y[k] * spg1_35[k];

        t_81[k] = f_3 * pc_x[k] * ppf_56[k];

        t_82[k] = f_3 * pc_x[k] * ppf_57[k];

        t_83[k] = f_3 * pc_x[k] * ppf_58[k];

        t_84[k] = f_3 * pc_x[k] * ppf_59[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, pa_y, pb_z, pc_y, pc_z, spg0_42, spf_28, spg1_42, \
                         psg0_25, psf_16, psg1_25, ppf_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_z[k] * psg0_25[k]
                  - f_6 * pc_z[k] * psg1_25[k];

        t_86[k] = f_0 * psf_16[k]
                  + f_3 * pc_z[k] * ppf_56[k];

        t_87[k] = pa_y[k] * spg0_42[k]
                  + f_7 * spf_28[k]
                  - f_6 * pc_y[k] * spg1_42[k];
    }

#pragma omp simd aligned(t_88, t_89, t_90, t_91, pa_y, pa_z, pc_y, pc_z, spg0_0, spg0_44, \
                         spf_29, spg1_0, spg1_44, ppf_59, ppf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_0 * spf_29[k]
                  + f_3 * pc_y[k] * ppf_59[k];

        t_89[k] = pa_y[k] * spg0_44[k]
                  - f_6 * pc_y[k] * spg1_44[k];

        t_90[k] = pa_z[k] * spg0_0[k]
                  - f_6 * pc_z[k] * spg1_0[k];

        t_91[k] = f_3 * pc_y[k] * ppf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_z, pc_y, pc_z, spg0_2, spg0_3, spg0_5, \
                         spf_0, spf_2, spg1_2, spg1_3, spg1_5, ppf_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pa_z[k] * spg0_2[k]
                  + f_0 * spf_0[k]
                  - f_6 * pc_z[k] * spg1_2[k];

        t_93[k] = pa_z[k] * spg0_3[k]
                  - f_6 * pc_z[k] * spg1_3[k];

        t_94[k] = f_3 * pc_y[k] * ppf_62[k];

        t_95[k] = pa_z[k] * spg0_5[k]
                  + f_7 * spf_2[k]
                  - f_6 * pc_z[k] * spg1_5[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, pc_x, psf_26, psf_27, psf_28, psf_29, ppf_66, \
                         ppf_67, ppf_68, ppf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_0 * psf_26[k]
                  + f_3 * pc_x[k] * ppf_66[k];

        t_97[k] = f_0 * psf_27[k]
                  + f_3 * pc_x[k] * ppf_67[k];

        t_98[k] = f_0 * psf_28[k]
                  + f_3 * pc_x[k] * ppf_68[k];

        t_99[k] = f_0 * psf_29[k]
                  + f_3 * pc_x[k] * ppf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, t_103, pa_z, pb_x, pc_x, pc_y, pc_z, spg0_10, \
                         spg1_10, psg0_41, psg0_42, psg1_41, psg1_42, \
                         ppf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = pa_z[k] * spg0_10[k]
                   - f_6 * pc_z[k] * spg1_10[k];

        t_101[k] = pb_x[k] * psg0_41[k]
                   - f_6 * pc_x[k] * psg1_41[k];

        t_102[k] = pb_x[k] * psg0_42[k]
                   - f_6 * pc_x[k] * psg1_42[k];

        t_103[k] = f_3 * pc_y[k] * ppf_69[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_z, pb_x, pc_x, pc_y, pc_z, spg0_15, spg1_15, \
                         psg0_44, psf_20, psg1_44, ppf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pb_x[k] * psg0_44[k]
                   - f_6 * pc_x[k] * psg1_44[k];

        t_105[k] = pa_z[k] * spg0_15[k]
                   - f_6 * pc_z[k] * spg1_15[k];

        t_106[k] = f_0 * psf_20[k]
                   + f_3 * pc_y[k] * ppf_70[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, t_110, pa_z, pb_y, pc_y, pc_z, spg0_18, spg1_18, \
                         psg0_32, psg0_35, psf_22, psg1_32, psg1_35, \
                         ppf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = pb_y[k] * psg0_32[k]
                   - f_6 * pc_y[k] * psg1_32[k];

        t_108[k] = pa_z[k] * spg0_18[k]
                   - f_6 * pc_z[k] * spg1_18[k];

        t_109[k] = f_0 * psf_22[k]
                   + f_3 * pc_y[k] * ppf_72[k];

        t_110[k] = pb_y[k] * psg0_35[k]
                   - f_6 * pc_y[k] * psg1_35[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, t_115, pa_z, pc_x, pc_z, spg0_25, \
                         spg1_25, ppf_76, ppf_77, ppf_78, ppf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_3 * pc_x[k] * ppf_76[k];

        t_112[k] = f_3 * pc_x[k] * ppf_77[k];

        t_113[k] = f_3 * pc_x[k] * ppf_78[k];

        t_114[k] = f_3 * pc_x[k] * ppf_79[k];

        t_115[k] = pa_z[k] * spg0_25[k]
                   - f_6 * pc_z[k] * spg1_25[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, pa_z, pc_y, pc_z, spg0_26, spg0_27, spf_16, \
                         spf_17, spg1_26, spg1_27, psf_29, ppf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_z[k] * spg0_26[k]
                   + f_0 * spf_16[k]
                   - f_6 * pc_z[k] * spg1_26[k];

        t_117[k] = pa_z[k] * spg0_27[k]
                   + f_7 * spf_17[k]
                   - f_6 * pc_z[k] * spg1_27[k];

        t_118[k] = f_0 * psf_29[k]
                   + f_3 * pc_y[k] * ppf_79[k];
    }

#pragma omp simd aligned(t_119, t_120, t_121, t_122, pb_y, pc_x, pc_y, psg0_44, psg1_44, \
                         ppd0_48, ppd0_50, ppd1_48, ppd1_50, ppf_80, \
                         ppf_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_119[k] = pb_y[k] * psg0_44[k]
                   - f_6 * pc_y[k] * psg1_44[k];

        t_120[k] = f_1 * ppd0_48[k]
                   - f_2 * ppd1_48[k]
                   + f_3 * pc_x[k] * ppf_80[k];

        t_121[k] = f_3 * pc_y[k] * ppf_80[k];

        t_122[k] = f_8 * ppd0_50[k]
                   - f_9 * ppd1_50[k]
                   + f_3 * pc_x[k] * ppf_82[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, t_126, t_127, pc_x, pc_y, ppd0_51, ppd0_53, \
                         ppd1_51, ppd1_53, ppf_82, ppf_83, ppf_85, ppf_86, \
                         ppf_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_4 * ppd0_51[k]
                   - f_5 * ppd1_51[k]
                   + f_3 * pc_x[k] * ppf_83[k];

        t_124[k] = f_3 * pc_y[k] * ppf_82[k];

        t_125[k] = f_4 * ppd0_53[k]
                   - f_5 * ppd1_53[k]
                   + f_3 * pc_x[k] * ppf_85[k];

        t_126[k] = f_3 * pc_x[k] * ppf_86[k];

        t_127[k] = f_3 * pc_x[k] * ppf_87[k];
    }

#pragma omp simd aligned(t_128, t_129, t_130, t_131, pc_x, pc_y, ppd0_51, ppd0_52, ppd1_51, \
                         ppd1_52, ppf_86, ppf_87, ppf_88, ppf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_128[k] = f_3 * pc_x[k] * ppf_88[k];

        t_129[k] = f_3 * pc_x[k] * ppf_89[k];

        t_130[k] = f_1 * ppd0_51[k]
                   - f_2 * ppd1_51[k]
                   + f_3 * pc_y[k] * ppf_86[k];

        t_131[k] = f_8 * ppd0_52[k]
                   - f_9 * ppd1_52[k]
                   + f_3 * pc_y[k] * ppf_87[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, pc_y, pc_z, spf_29, psf_29, ppd0_53, ppd1_53, \
                         ppf_88, ppf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = f_4 * ppd0_53[k]
                   - f_5 * ppd1_53[k]
                   + f_3 * pc_y[k] * ppf_88[k];

        t_133[k] = f_3 * pc_y[k] * ppf_89[k];

        t_134[k] = f_0 * spf_29[k]
                   + f_0 * psf_29[k]
                   + f_1 * ppd0_53[k]
                   - f_2 * ppd1_53[k]
                   + f_3 * pc_z[k] * ppf_89[k];
    }
}

}  // namespace simdt3ceri
