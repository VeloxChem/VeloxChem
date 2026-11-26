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


#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_pfg_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfg0, const size_t sff,
                                                          const size_t sfg1, const size_t pdg0,
                                                          const size_t pdf, const size_t pdg1,
                                                          const size_t pfd0, const size_t pfd1,
                                                          const size_t pff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;

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

    const auto *pa_x = buffer.data(pa + 0);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg0_90 = buffer.data(sfg0 + 90);
    const auto *sfg0_93 = buffer.data(sfg0 + 93);
    const auto *sfg0_100 = buffer.data(sfg0 + 100);
    const auto *sfg0_102 = buffer.data(sfg0 + 102);
    const auto *sfg0_104 = buffer.data(sfg0 + 104);
    const auto *sfg0_110 = buffer.data(sfg0 + 110);
    const auto *sfg0_115 = buffer.data(sfg0 + 115);
    const auto *sfg0_117 = buffer.data(sfg0 + 117);
    const auto *sfg0_119 = buffer.data(sfg0 + 119);
    const auto *sfg0_123 = buffer.data(sfg0 + 123);

    const auto *sff_0 = buffer.data(sff + 0);
    const auto *sff_6 = buffer.data(sff + 6);
    const auto *sff_9 = buffer.data(sff + 9);
    const auto *sff_16 = buffer.data(sff + 16);
    const auto *sff_29 = buffer.data(sff + 29);
    const auto *sff_30 = buffer.data(sff + 30);
    const auto *sff_36 = buffer.data(sff + 36);
    const auto *sff_39 = buffer.data(sff + 39);
    const auto *sff_50 = buffer.data(sff + 50);
    const auto *sff_56 = buffer.data(sff + 56);
    const auto *sff_59 = buffer.data(sff + 59);
    const auto *sff_60 = buffer.data(sff + 60);
    const auto *sff_63 = buffer.data(sff + 63);
    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_75 = buffer.data(sff + 75);
    const auto *sff_76 = buffer.data(sff + 76);
    const auto *sff_79 = buffer.data(sff + 79);
    const auto *sff_83 = buffer.data(sff + 83);

    const auto *sfg1_90 = buffer.data(sfg1 + 90);
    const auto *sfg1_93 = buffer.data(sfg1 + 93);
    const auto *sfg1_100 = buffer.data(sfg1 + 100);
    const auto *sfg1_102 = buffer.data(sfg1 + 102);
    const auto *sfg1_104 = buffer.data(sfg1 + 104);
    const auto *sfg1_110 = buffer.data(sfg1 + 110);
    const auto *sfg1_115 = buffer.data(sfg1 + 115);
    const auto *sfg1_117 = buffer.data(sfg1 + 117);
    const auto *sfg1_119 = buffer.data(sfg1 + 119);
    const auto *sfg1_123 = buffer.data(sfg1 + 123);

    const auto *pdg0_0 = buffer.data(pdg0 + 0);
    const auto *pdg0_3 = buffer.data(pdg0 + 3);
    const auto *pdg0_5 = buffer.data(pdg0 + 5);
    const auto *pdg0_6 = buffer.data(pdg0 + 6);
    const auto *pdg0_9 = buffer.data(pdg0 + 9);
    const auto *pdg0_10 = buffer.data(pdg0 + 10);
    const auto *pdg0_14 = buffer.data(pdg0 + 14);
    const auto *pdg0_18 = buffer.data(pdg0 + 18);
    const auto *pdg0_21 = buffer.data(pdg0 + 21);
    const auto *pdg0_30 = buffer.data(pdg0 + 30);
    const auto *pdg0_35 = buffer.data(pdg0 + 35);
    const auto *pdg0_39 = buffer.data(pdg0 + 39);
    const auto *pdg0_45 = buffer.data(pdg0 + 45);
    const auto *pdg0_48 = buffer.data(pdg0 + 48);
    const auto *pdg0_75 = buffer.data(pdg0 + 75);

    const auto *pdf_0 = buffer.data(pdf + 0);
    const auto *pdf_1 = buffer.data(pdf + 1);
    const auto *pdf_2 = buffer.data(pdf + 2);
    const auto *pdf_3 = buffer.data(pdf + 3);
    const auto *pdf_5 = buffer.data(pdf + 5);
    const auto *pdf_6 = buffer.data(pdf + 6);
    const auto *pdf_8 = buffer.data(pdf + 8);
    const auto *pdf_9 = buffer.data(pdf + 9);
    const auto *pdf_10 = buffer.data(pdf + 10);
    const auto *pdf_11 = buffer.data(pdf + 11);
    const auto *pdf_12 = buffer.data(pdf + 12);
    const auto *pdf_13 = buffer.data(pdf + 13);
    const auto *pdf_15 = buffer.data(pdf + 15);
    const auto *pdf_16 = buffer.data(pdf + 16);
    const auto *pdf_18 = buffer.data(pdf + 18);
    const auto *pdf_19 = buffer.data(pdf + 19);
    const auto *pdf_20 = buffer.data(pdf + 20);
    const auto *pdf_22 = buffer.data(pdf + 22);
    const auto *pdf_23 = buffer.data(pdf + 23);
    const auto *pdf_25 = buffer.data(pdf + 25);
    const auto *pdf_26 = buffer.data(pdf + 26);
    const auto *pdf_28 = buffer.data(pdf + 28);
    const auto *pdf_29 = buffer.data(pdf + 29);
    const auto *pdf_30 = buffer.data(pdf + 30);
    const auto *pdf_32 = buffer.data(pdf + 32);
    const auto *pdf_33 = buffer.data(pdf + 33);
    const auto *pdf_35 = buffer.data(pdf + 35);
    const auto *pdf_36 = buffer.data(pdf + 36);
    const auto *pdf_39 = buffer.data(pdf + 39);
    const auto *pdf_40 = buffer.data(pdf + 40);
    const auto *pdf_42 = buffer.data(pdf + 42);
    const auto *pdf_45 = buffer.data(pdf + 45);
    const auto *pdf_49 = buffer.data(pdf + 49);
    const auto *pdf_50 = buffer.data(pdf + 50);
    const auto *pdf_52 = buffer.data(pdf + 52);
    const auto *pdf_56 = buffer.data(pdf + 56);
    const auto *pdf_59 = buffer.data(pdf + 59);

    const auto *pdg1_0 = buffer.data(pdg1 + 0);
    const auto *pdg1_3 = buffer.data(pdg1 + 3);
    const auto *pdg1_5 = buffer.data(pdg1 + 5);
    const auto *pdg1_6 = buffer.data(pdg1 + 6);
    const auto *pdg1_9 = buffer.data(pdg1 + 9);
    const auto *pdg1_10 = buffer.data(pdg1 + 10);
    const auto *pdg1_14 = buffer.data(pdg1 + 14);
    const auto *pdg1_18 = buffer.data(pdg1 + 18);
    const auto *pdg1_21 = buffer.data(pdg1 + 21);
    const auto *pdg1_30 = buffer.data(pdg1 + 30);
    const auto *pdg1_35 = buffer.data(pdg1 + 35);
    const auto *pdg1_39 = buffer.data(pdg1 + 39);
    const auto *pdg1_45 = buffer.data(pdg1 + 45);
    const auto *pdg1_48 = buffer.data(pdg1 + 48);
    const auto *pdg1_75 = buffer.data(pdg1 + 75);

    const auto *pfd0_0 = buffer.data(pfd0 + 0);
    const auto *pfd0_3 = buffer.data(pfd0 + 3);
    const auto *pfd0_5 = buffer.data(pfd0 + 5);
    const auto *pfd0_9 = buffer.data(pfd0 + 9);
    const auto *pfd0_11 = buffer.data(pfd0 + 11);
    const auto *pfd0_17 = buffer.data(pfd0 + 17);
    const auto *pfd0_18 = buffer.data(pfd0 + 18);
    const auto *pfd0_21 = buffer.data(pfd0 + 21);
    const auto *pfd0_23 = buffer.data(pfd0 + 23);
    const auto *pfd0_27 = buffer.data(pfd0 + 27);
    const auto *pfd0_29 = buffer.data(pfd0 + 29);
    const auto *pfd0_30 = buffer.data(pfd0 + 30);
    const auto *pfd0_33 = buffer.data(pfd0 + 33);
    const auto *pfd0_35 = buffer.data(pfd0 + 35);
    const auto *pfd0_36 = buffer.data(pfd0 + 36);

    const auto *pfd1_0 = buffer.data(pfd1 + 0);
    const auto *pfd1_3 = buffer.data(pfd1 + 3);
    const auto *pfd1_5 = buffer.data(pfd1 + 5);
    const auto *pfd1_9 = buffer.data(pfd1 + 9);
    const auto *pfd1_11 = buffer.data(pfd1 + 11);
    const auto *pfd1_17 = buffer.data(pfd1 + 17);
    const auto *pfd1_18 = buffer.data(pfd1 + 18);
    const auto *pfd1_21 = buffer.data(pfd1 + 21);
    const auto *pfd1_23 = buffer.data(pfd1 + 23);
    const auto *pfd1_27 = buffer.data(pfd1 + 27);
    const auto *pfd1_29 = buffer.data(pfd1 + 29);
    const auto *pfd1_30 = buffer.data(pfd1 + 30);
    const auto *pfd1_33 = buffer.data(pfd1 + 33);
    const auto *pfd1_35 = buffer.data(pfd1 + 35);
    const auto *pfd1_36 = buffer.data(pfd1 + 36);

    const auto *pff_0 = buffer.data(pff + 0);
    const auto *pff_1 = buffer.data(pff + 1);
    const auto *pff_2 = buffer.data(pff + 2);
    const auto *pff_3 = buffer.data(pff + 3);
    const auto *pff_5 = buffer.data(pff + 5);
    const auto *pff_6 = buffer.data(pff + 6);
    const auto *pff_8 = buffer.data(pff + 8);
    const auto *pff_9 = buffer.data(pff + 9);
    const auto *pff_10 = buffer.data(pff + 10);
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
    const auto *pff_43 = buffer.data(pff + 43);
    const auto *pff_45 = buffer.data(pff + 45);
    const auto *pff_46 = buffer.data(pff + 46);
    const auto *pff_48 = buffer.data(pff + 48);
    const auto *pff_49 = buffer.data(pff + 49);
    const auto *pff_50 = buffer.data(pff + 50);
    const auto *pff_51 = buffer.data(pff + 51);
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

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pc_x, pc_y, pc_z, sff_0, pdf_0, pfd0_0, \
                         pfd1_0, pff_0, pff_1, pff_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * sff_0[k]
                 + f_1 * pdf_0[k]
                 + f_2 * pfd0_0[k]
                 - f_3 * pfd1_0[k]
                 + f_4 * pc_x[k] * pff_0[k];

        t_1[k] = f_4 * pc_y[k] * pff_0[k];

        t_2[k] = f_4 * pc_z[k] * pff_0[k];

        t_3[k] = f_5 * pfd0_0[k]
                 - f_6 * pfd1_0[k]
                 + f_4 * pc_y[k] * pff_1[k];

        t_4[k] = f_4 * pc_y[k] * pff_2[k];

        t_5[k] = f_5 * pfd0_0[k]
                 - f_6 * pfd1_0[k]
                 + f_4 * pc_z[k] * pff_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pc_x, pc_y, pc_z, sff_6, sff_9, pdf_6, pdf_9, \
                         pff_3, pff_5, pff_6, pff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_0 * sff_6[k]
                 + f_1 * pdf_6[k]
                 + f_4 * pc_x[k] * pff_6[k];

        t_7[k] = f_4 * pc_z[k] * pff_3[k];

        t_8[k] = f_4 * pc_y[k] * pff_5[k];

        t_9[k] = f_0 * sff_9[k]
                 + f_1 * pdf_9[k]
                 + f_4 * pc_x[k] * pff_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, t_13, t_14, pc_y, pc_z, pfd0_3, pfd0_5, pfd1_3, \
                         pfd1_5, pff_6, pff_8, pff_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = f_2 * pfd0_3[k]
                  - f_3 * pfd1_3[k]
                  + f_4 * pc_y[k] * pff_6[k];

        t_11[k] = f_4 * pc_z[k] * pff_6[k];

        t_12[k] = f_5 * pfd0_5[k]
                  - f_6 * pfd1_5[k]
                  + f_4 * pc_y[k] * pff_8[k];

        t_13[k] = f_4 * pc_y[k] * pff_9[k];

        t_14[k] = f_2 * pfd0_5[k]
                  - f_3 * pfd1_5[k]
                  + f_4 * pc_z[k] * pff_9[k];
    }

#pragma omp simd aligned(t_15, t_16, t_17, t_18, pb_y, pc_y, pc_z, pdg0_0, pdg0_3, pdf_0, \
                         pdf_1, pdg1_0, pdg1_3, pff_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = pb_y[k] * pdg0_0[k]
                  - f_7 * pc_y[k] * pdg1_0[k];

        t_16[k] = f_0 * pdf_0[k]
                  + f_4 * pc_y[k] * pff_10[k];

        t_17[k] = f_4 * pc_z[k] * pff_10[k];

        t_18[k] = pb_y[k] * pdg0_3[k]
                  + f_8 * pdf_1[k]
                  - f_7 * pc_y[k] * pdg1_3[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pb_y, pc_x, pc_y, pc_z, sff_16, pdg0_5, \
                         pdf_2, pdf_16, pdg1_5, pff_12, pff_13, \
                         pff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = f_0 * pdf_2[k]
                  + f_4 * pc_y[k] * pff_12[k];

        t_20[k] = pb_y[k] * pdg0_5[k]
                  - f_7 * pc_y[k] * pdg1_5[k];

        t_21[k] = f_0 * sff_16[k]
                  + f_8 * pdf_16[k]
                  + f_4 * pc_x[k] * pff_16[k];

        t_22[k] = f_4 * pc_z[k] * pff_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pb_y, pc_y, pc_z, pdg0_9, pdf_5, pdf_6, \
                         pdg1_9, pfd0_9, pfd1_9, pff_15, pff_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_0 * pdf_5[k]
                  + f_4 * pc_y[k] * pff_15[k];

        t_24[k] = pb_y[k] * pdg0_9[k]
                  - f_7 * pc_y[k] * pdg1_9[k];

        t_25[k] = f_0 * pdf_6[k]
                  + f_2 * pfd0_9[k]
                  - f_3 * pfd1_9[k]
                  + f_4 * pc_y[k] * pff_16[k];

        t_26[k] = f_4 * pc_z[k] * pff_16[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pb_y, pc_y, pdg0_14, pdf_8, pdf_9, pdg1_14, \
                         pfd0_11, pfd1_11, pff_18, pff_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_0 * pdf_8[k]
                  + f_5 * pfd0_11[k]
                  - f_6 * pfd1_11[k]
                  + f_4 * pc_y[k] * pff_18[k];

        t_28[k] = f_0 * pdf_9[k]
                  + f_4 * pc_y[k] * pff_19[k];

        t_29[k] = pb_y[k] * pdg0_14[k]
                  - f_7 * pc_y[k] * pdg1_14[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, t_34, pb_z, pc_y, pc_z, pdg0_0, pdg0_3, \
                         pdf_0, pdg1_0, pdg1_3, pff_20, pff_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pb_z[k] * pdg0_0[k]
                  - f_7 * pc_z[k] * pdg1_0[k];

        t_31[k] = f_4 * pc_y[k] * pff_20[k];

        t_32[k] = f_0 * pdf_0[k]
                  + f_4 * pc_z[k] * pff_20[k];

        t_33[k] = pb_z[k] * pdg0_3[k]
                  - f_7 * pc_z[k] * pdg1_3[k];

        t_34[k] = f_4 * pc_y[k] * pff_22[k];
    }

#pragma omp simd aligned(t_35, t_36, t_37, t_38, pb_z, pc_y, pc_z, pdg0_5, pdg0_6, pdf_2, \
                         pdf_3, pdg1_5, pdg1_6, pff_23, pff_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = pb_z[k] * pdg0_5[k]
                  + f_8 * pdf_2[k]
                  - f_7 * pc_z[k] * pdg1_5[k];

        t_36[k] = pb_z[k] * pdg0_6[k]
                  - f_7 * pc_z[k] * pdg1_6[k];

        t_37[k] = f_0 * pdf_3[k]
                  + f_4 * pc_z[k] * pff_23[k];

        t_38[k] = f_4 * pc_y[k] * pff_25[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pb_z, pc_x, pc_z, sff_29, pdg0_10, pdf_6, pdf_29, \
                         pdg1_10, pff_26, pff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_0 * sff_29[k]
                  + f_8 * pdf_29[k]
                  + f_4 * pc_x[k] * pff_29[k];

        t_40[k] = pb_z[k] * pdg0_10[k]
                  - f_7 * pc_z[k] * pdg1_10[k];

        t_41[k] = f_0 * pdf_6[k]
                  + f_4 * pc_z[k] * pff_26[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, pc_y, pc_z, pdf_9, pfd0_17, pfd1_17, pff_28, \
                         pff_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = f_5 * pfd0_17[k]
                  - f_6 * pfd1_17[k]
                  + f_4 * pc_y[k] * pff_28[k];

        t_43[k] = f_4 * pc_y[k] * pff_29[k];

        t_44[k] = f_0 * pdf_9[k]
                  + f_2 * pfd0_17[k]
                  - f_3 * pfd1_17[k]
                  + f_4 * pc_z[k] * pff_29[k];
    }

#pragma omp simd aligned(t_45, t_46, t_47, t_48, pc_x, pc_y, pc_z, sff_30, pdf_10, pdf_11, \
                         pdf_30, pfd0_18, pfd1_18, pff_30, pff_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_0 * sff_30[k]
                  + f_0 * pdf_30[k]
                  + f_2 * pfd0_18[k]
                  - f_3 * pfd1_18[k]
                  + f_4 * pc_x[k] * pff_30[k];

        t_46[k] = f_8 * pdf_10[k]
                  + f_4 * pc_y[k] * pff_30[k];

        t_47[k] = f_4 * pc_z[k] * pff_30[k];

        t_48[k] = f_8 * pdf_11[k]
                  + f_5 * pfd0_18[k]
                  - f_6 * pfd1_18[k]
                  + f_4 * pc_y[k] * pff_31[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, pc_x, pc_y, pc_z, sff_36, pdf_12, pdf_36, \
                         pfd0_18, pfd1_18, pff_32, pff_33, pff_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_8 * pdf_12[k]
                  + f_4 * pc_y[k] * pff_32[k];

        t_50[k] = f_5 * pfd0_18[k]
                  - f_6 * pfd1_18[k]
                  + f_4 * pc_z[k] * pff_32[k];

        t_51[k] = f_0 * sff_36[k]
                  + f_0 * pdf_36[k]
                  + f_4 * pc_x[k] * pff_36[k];

        t_52[k] = f_4 * pc_z[k] * pff_33[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, t_56, pc_x, pc_y, pc_z, sff_39, pdf_15, pdf_16, \
                         pdf_39, pfd0_21, pfd1_21, pff_35, pff_36, \
                         pff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_8 * pdf_15[k]
                  + f_4 * pc_y[k] * pff_35[k];

        t_54[k] = f_0 * sff_39[k]
                  + f_0 * pdf_39[k]
                  + f_4 * pc_x[k] * pff_39[k];

        t_55[k] = f_8 * pdf_16[k]
                  + f_2 * pfd0_21[k]
                  - f_3 * pfd1_21[k]
                  + f_4 * pc_y[k] * pff_36[k];

        t_56[k] = f_4 * pc_z[k] * pff_36[k];
    }

#pragma omp simd aligned(t_57, t_58, t_59, t_60, pb_y, pc_y, pc_z, pdg0_30, pdf_18, pdf_19, \
                         pdg1_30, pfd0_23, pfd1_23, pff_38, pff_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_8 * pdf_18[k]
                  + f_5 * pfd0_23[k]
                  - f_6 * pfd1_23[k]
                  + f_4 * pc_y[k] * pff_38[k];

        t_58[k] = f_8 * pdf_19[k]
                  + f_4 * pc_y[k] * pff_39[k];

        t_59[k] = f_2 * pfd0_23[k]
                  - f_3 * pfd1_23[k]
                  + f_4 * pc_z[k] * pff_39[k];

        t_60[k] = pb_y[k] * pdg0_30[k]
                  - f_7 * pc_y[k] * pdg1_30[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_z, pc_y, pc_z, pdg0_18, pdf_10, pdf_20, \
                         pdf_22, pdg1_18, pff_40, pff_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_0 * pdf_20[k]
                  + f_4 * pc_y[k] * pff_40[k];

        t_62[k] = f_0 * pdf_10[k]
                  + f_4 * pc_z[k] * pff_40[k];

        t_63[k] = pb_z[k] * pdg0_18[k]
                  - f_7 * pc_z[k] * pdg1_18[k];

        t_64[k] = f_0 * pdf_22[k]
                  + f_4 * pc_y[k] * pff_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pb_y, pb_z, pc_y, pc_z, pdg0_21, pdg0_35, \
                         pdf_13, pdf_25, pdg1_21, pdg1_35, pff_43, \
                         pff_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pb_y[k] * pdg0_35[k]
                  - f_7 * pc_y[k] * pdg1_35[k];

        t_66[k] = pb_z[k] * pdg0_21[k]
                  - f_7 * pc_z[k] * pdg1_21[k];

        t_67[k] = f_0 * pdf_13[k]
                  + f_4 * pc_z[k] * pff_43[k];

        t_68[k] = f_0 * pdf_25[k]
                  + f_4 * pc_y[k] * pff_45[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, pb_y, pc_y, pc_z, pdg0_39, pdf_16, pdf_26, pdg1_39, \
                         pfd0_27, pfd1_27, pff_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * pdg0_39[k]
                  - f_7 * pc_y[k] * pdg1_39[k];

        t_70[k] = f_0 * pdf_26[k]
                  + f_2 * pfd0_27[k]
                  - f_3 * pfd1_27[k]
                  + f_4 * pc_y[k] * pff_46[k];

        t_71[k] = f_0 * pdf_16[k]
                  + f_4 * pc_z[k] * pff_46[k];
    }

#pragma omp simd aligned(t_72, t_73, t_74, pc_y, pc_z, pdf_19, pdf_28, pdf_29, pfd0_29, \
                         pfd1_29, pff_48, pff_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_72[k] = f_0 * pdf_28[k]
                  + f_5 * pfd0_29[k]
                  - f_6 * pfd1_29[k]
                  + f_4 * pc_y[k] * pff_48[k];

        t_73[k] = f_0 * pdf_29[k]
                  + f_4 * pc_y[k] * pff_49[k];

        t_74[k] = f_0 * pdf_19[k]
                  + f_2 * pfd0_29[k]
                  - f_3 * pfd1_29[k]
                  + f_4 * pc_z[k] * pff_49[k];
    }

#pragma omp simd aligned(t_75, t_76, t_77, t_78, t_79, pc_x, pc_y, pc_z, sff_50, pdf_20, \
                         pdf_50, pfd0_30, pfd1_30, pff_50, pff_51, \
                         pff_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_75[k] = f_0 * sff_50[k]
                  + f_0 * pdf_50[k]
                  + f_2 * pfd0_30[k]
                  - f_3 * pfd1_30[k]
                  + f_4 * pc_x[k] * pff_50[k];

        t_76[k] = f_4 * pc_y[k] * pff_50[k];

        t_77[k] = f_8 * pdf_20[k]
                  + f_4 * pc_z[k] * pff_50[k];

        t_78[k] = f_5 * pfd0_30[k]
                  - f_6 * pfd1_30[k]
                  + f_4 * pc_y[k] * pff_51[k];

        t_79[k] = f_4 * pc_y[k] * pff_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pc_x, pc_z, sff_56, pdf_22, pdf_23, pdf_56, \
                         pfd0_30, pfd1_30, pff_52, pff_53, pff_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_8 * pdf_22[k]
                  + f_5 * pfd0_30[k]
                  - f_6 * pfd1_30[k]
                  + f_4 * pc_z[k] * pff_52[k];

        t_81[k] = f_0 * sff_56[k]
                  + f_0 * pdf_56[k]
                  + f_4 * pc_x[k] * pff_56[k];

        t_82[k] = f_8 * pdf_23[k]
                  + f_4 * pc_z[k] * pff_53[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, t_86, pc_x, pc_y, pc_z, sff_59, pdf_26, pdf_59, \
                         pfd0_33, pfd1_33, pff_55, pff_56, pff_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_4 * pc_y[k] * pff_55[k];

        t_84[k] = f_0 * sff_59[k]
                  + f_0 * pdf_59[k]
                  + f_4 * pc_x[k] * pff_59[k];

        t_85[k] = f_2 * pfd0_33[k]
                  - f_3 * pfd1_33[k]
                  + f_4 * pc_y[k] * pff_56[k];

        t_86[k] = f_8 * pdf_26[k]
                  + f_4 * pc_z[k] * pff_56[k];
    }

#pragma omp simd aligned(t_87, t_88, t_89, t_90, pa_x, pc_x, pc_y, pc_z, sfg0_90, sff_60, \
                         sfg1_90, pdf_29, pfd0_35, pfd1_35, pff_58, \
                         pff_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_87[k] = f_5 * pfd0_35[k]
                  - f_6 * pfd1_35[k]
                  + f_4 * pc_y[k] * pff_58[k];

        t_88[k] = f_4 * pc_y[k] * pff_59[k];

        t_89[k] = f_8 * pdf_29[k]
                  + f_2 * pfd0_35[k]
                  - f_3 * pfd1_35[k]
                  + f_4 * pc_z[k] * pff_59[k];

        t_90[k] = pa_x[k] * sfg0_90[k]
                  + f_9 * sff_60[k]
                  - f_7 * pc_x[k] * sfg1_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, t_94, pa_x, pc_x, pc_y, pc_z, sfg0_93, sff_63, \
                         sfg1_93, pdf_30, pdf_32, pff_60, pff_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_1 * pdf_30[k]
                  + f_4 * pc_y[k] * pff_60[k];

        t_92[k] = f_4 * pc_z[k] * pff_60[k];

        t_93[k] = pa_x[k] * sfg0_93[k]
                  + f_8 * sff_63[k]
                  - f_7 * pc_x[k] * sfg1_93[k];

        t_94[k] = f_1 * pdf_32[k]
                  + f_4 * pc_y[k] * pff_62[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, t_98, pc_x, pc_y, pc_z, sff_66, pdf_35, pfd0_36, \
                         pfd1_36, pff_62, pff_63, pff_65, pff_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_5 * pfd0_36[k]
                  - f_6 * pfd1_36[k]
                  + f_4 * pc_z[k] * pff_62[k];

        t_96[k] = f_0 * sff_66[k]
                  + f_4 * pc_x[k] * pff_66[k];

        t_97[k] = f_4 * pc_z[k] * pff_63[k];

        t_98[k] = f_1 * pdf_35[k]
                  + f_4 * pc_y[k] * pff_65[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pc_x, pc_z, sfg0_100, sfg0_102, \
                         sff_69, sfg1_100, sfg1_102, pff_66, pff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_0 * sff_69[k]
                  + f_4 * pc_x[k] * pff_69[k];

        t_100[k] = pa_x[k] * sfg0_100[k]
                   - f_7 * pc_x[k] * sfg1_100[k];

        t_101[k] = f_4 * pc_z[k] * pff_66[k];

        t_102[k] = pa_x[k] * sfg0_102[k]
                   - f_7 * pc_x[k] * sfg1_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, pa_x, pb_z, pc_x, pc_y, pc_z, sfg0_104, \
                         sfg1_104, pdg0_45, pdf_39, pdg1_45, pff_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_1 * pdf_39[k]
                   + f_4 * pc_y[k] * pff_69[k];

        t_104[k] = pa_x[k] * sfg0_104[k]
                   - f_7 * pc_x[k] * sfg1_104[k];

        t_105[k] = pb_z[k] * pdg0_45[k]
                   - f_7 * pc_z[k] * pdg1_45[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, t_109, pb_z, pc_y, pc_z, pdg0_48, pdf_30, \
                         pdf_40, pdf_42, pdg1_48, pff_70, pff_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = f_8 * pdf_40[k]
                   + f_4 * pc_y[k] * pff_70[k];

        t_107[k] = f_0 * pdf_30[k]
                   + f_4 * pc_z[k] * pff_70[k];

        t_108[k] = pb_z[k] * pdg0_48[k]
                   - f_7 * pc_z[k] * pdg1_48[k];

        t_109[k] = f_8 * pdf_42[k]
                   + f_4 * pc_y[k] * pff_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_x, pc_x, pc_z, sfg0_110, sff_75, sff_76, \
                         sfg1_110, pdf_33, pff_73, pff_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = pa_x[k] * sfg0_110[k]
                   + f_8 * sff_75[k]
                   - f_7 * pc_x[k] * sfg1_110[k];

        t_111[k] = f_0 * sff_76[k]
                   + f_4 * pc_x[k] * pff_76[k];

        t_112[k] = f_0 * pdf_33[k]
                   + f_4 * pc_z[k] * pff_73[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pc_x, pc_y, pc_z, sfg0_115, sff_79, \
                         sfg1_115, pdf_36, pdf_45, pff_75, pff_76, \
                         pff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_8 * pdf_45[k]
                   + f_4 * pc_y[k] * pff_75[k];

        t_114[k] = f_0 * sff_79[k]
                   + f_4 * pc_x[k] * pff_79[k];

        t_115[k] = pa_x[k] * sfg0_115[k]
                   - f_7 * pc_x[k] * sfg1_115[k];

        t_116[k] = f_0 * pdf_36[k]
                   + f_4 * pc_z[k] * pff_76[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pb_y, pc_x, pc_y, sfg0_117, \
                         sfg0_119, sfg1_117, sfg1_119, pdg0_75, pdf_49, pdg1_75, \
                         pff_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_x[k] * sfg0_117[k]
                   - f_7 * pc_x[k] * sfg1_117[k];

        t_118[k] = f_8 * pdf_49[k]
                   + f_4 * pc_y[k] * pff_79[k];

        t_119[k] = pa_x[k] * sfg0_119[k]
                   - f_7 * pc_x[k] * sfg1_119[k];

        t_120[k] = pb_y[k] * pdg0_75[k]
                   - f_7 * pc_y[k] * pdg1_75[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, t_124, pa_x, pc_x, pc_y, pc_z, sfg0_123, sff_83, \
                         sfg1_123, pdf_40, pdf_50, pdf_52, pff_80, \
                         pff_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_0 * pdf_50[k]
                   + f_4 * pc_y[k] * pff_80[k];

        t_122[k] = f_8 * pdf_40[k]
                   + f_4 * pc_z[k] * pff_80[k];

        t_123[k] = pa_x[k] * sfg0_123[k]
                   + f_8 * sff_83[k]
                   - f_7 * pc_x[k] * sfg1_123[k];

        t_124[k] = f_0 * pdf_52[k]
                   + f_4 * pc_y[k] * pff_82[k];
    }
}

static auto
compute_prim_pfg_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfg0, const size_t sff,
                                                          const size_t sfg1, const size_t ppg0,
                                                          const size_t ppg1, const size_t pdg0,
                                                          const size_t pdf, const size_t pdg1,
                                                          const size_t pfd0, const size_t pfd1,
                                                          const size_t pff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg0_0 = buffer.data(sfg0 + 0);
    const auto *sfg0_1 = buffer.data(sfg0 + 1);
    const auto *sfg0_3 = buffer.data(sfg0 + 3);
    const auto *sfg0_5 = buffer.data(sfg0 + 5);
    const auto *sfg0_10 = buffer.data(sfg0 + 10);
    const auto *sfg0_14 = buffer.data(sfg0 + 14);
    const auto *sfg0_30 = buffer.data(sfg0 + 30);
    const auto *sfg0_35 = buffer.data(sfg0 + 35);
    const auto *sfg0_42 = buffer.data(sfg0 + 42);
    const auto *sfg0_44 = buffer.data(sfg0 + 44);
    const auto *sfg0_75 = buffer.data(sfg0 + 75);
    const auto *sfg0_76 = buffer.data(sfg0 + 76);
    const auto *sfg0_78 = buffer.data(sfg0 + 78);
    const auto *sfg0_80 = buffer.data(sfg0 + 80);
    const auto *sfg0_89 = buffer.data(sfg0 + 89);
    const auto *sfg0_130 = buffer.data(sfg0 + 130);
    const auto *sfg0_132 = buffer.data(sfg0 + 132);
    const auto *sfg0_134 = buffer.data(sfg0 + 134);
    const auto *sfg0_135 = buffer.data(sfg0 + 135);
    const auto *sfg0_140 = buffer.data(sfg0 + 140);
    const auto *sfg0_145 = buffer.data(sfg0 + 145);
    const auto *sfg0_147 = buffer.data(sfg0 + 147);
    const auto *sfg0_149 = buffer.data(sfg0 + 149);

    const auto *sff_0 = buffer.data(sff + 0);
    const auto *sff_1 = buffer.data(sff + 1);
    const auto *sff_6 = buffer.data(sff + 6);
    const auto *sff_9 = buffer.data(sff + 9);
    const auto *sff_19 = buffer.data(sff + 19);
    const auto *sff_28 = buffer.data(sff + 28);
    const auto *sff_29 = buffer.data(sff + 29);
    const auto *sff_50 = buffer.data(sff + 50);
    const auto *sff_51 = buffer.data(sff + 51);
    const auto *sff_59 = buffer.data(sff + 59);
    const auto *sff_86 = buffer.data(sff + 86);
    const auto *sff_89 = buffer.data(sff + 89);
    const auto *sff_90 = buffer.data(sff + 90);
    const auto *sff_95 = buffer.data(sff + 95);
    const auto *sff_96 = buffer.data(sff + 96);
    const auto *sff_99 = buffer.data(sff + 99);

    const auto *sfg1_0 = buffer.data(sfg1 + 0);
    const auto *sfg1_1 = buffer.data(sfg1 + 1);
    const auto *sfg1_3 = buffer.data(sfg1 + 3);
    const auto *sfg1_5 = buffer.data(sfg1 + 5);
    const auto *sfg1_10 = buffer.data(sfg1 + 10);
    const auto *sfg1_14 = buffer.data(sfg1 + 14);
    const auto *sfg1_30 = buffer.data(sfg1 + 30);
    const auto *sfg1_35 = buffer.data(sfg1 + 35);
    const auto *sfg1_42 = buffer.data(sfg1 + 42);
    const auto *sfg1_44 = buffer.data(sfg1 + 44);
    const auto *sfg1_75 = buffer.data(sfg1 + 75);
    const auto *sfg1_76 = buffer.data(sfg1 + 76);
    const auto *sfg1_78 = buffer.data(sfg1 + 78);
    const auto *sfg1_80 = buffer.data(sfg1 + 80);
    const auto *sfg1_89 = buffer.data(sfg1 + 89);
    const auto *sfg1_130 = buffer.data(sfg1 + 130);
    const auto *sfg1_132 = buffer.data(sfg1 + 132);
    const auto *sfg1_134 = buffer.data(sfg1 + 134);
    const auto *sfg1_135 = buffer.data(sfg1 + 135);
    const auto *sfg1_140 = buffer.data(sfg1 + 140);
    const auto *sfg1_145 = buffer.data(sfg1 + 145);
    const auto *sfg1_147 = buffer.data(sfg1 + 147);
    const auto *sfg1_149 = buffer.data(sfg1 + 149);

    const auto *ppg0_70 = buffer.data(ppg0 + 70);

    const auto *ppg1_70 = buffer.data(ppg1 + 70);

    const auto *pdg0_80 = buffer.data(pdg0 + 80);
    const auto *pdg0_91 = buffer.data(pdg0 + 91);
    const auto *pdg0_93 = buffer.data(pdg0 + 93);
    const auto *pdg0_100 = buffer.data(pdg0 + 100);
    const auto *pdg0_106 = buffer.data(pdg0 + 106);
    const auto *pdg0_108 = buffer.data(pdg0 + 108);
    const auto *pdg0_115 = buffer.data(pdg0 + 115);
    const auto *pdg0_135 = buffer.data(pdg0 + 135);
    const auto *pdg0_136 = buffer.data(pdg0 + 136);
    const auto *pdg0_138 = buffer.data(pdg0 + 138);
    const auto *pdg0_140 = buffer.data(pdg0 + 140);
    const auto *pdg0_145 = buffer.data(pdg0 + 145);
    const auto *pdg0_147 = buffer.data(pdg0 + 147);
    const auto *pdg0_148 = buffer.data(pdg0 + 148);
    const auto *pdg0_149 = buffer.data(pdg0 + 149);
    const auto *pdg0_160 = buffer.data(pdg0 + 160);
    const auto *pdg0_162 = buffer.data(pdg0 + 162);
    const auto *pdg0_163 = buffer.data(pdg0 + 163);
    const auto *pdg0_164 = buffer.data(pdg0 + 164);
    const auto *pdg0_175 = buffer.data(pdg0 + 175);
    const auto *pdg0_177 = buffer.data(pdg0 + 177);

    const auto *pdf_43 = buffer.data(pdf + 43);
    const auto *pdf_46 = buffer.data(pdf + 46);
    const auto *pdf_50 = buffer.data(pdf + 50);
    const auto *pdf_53 = buffer.data(pdf + 53);
    const auto *pdf_55 = buffer.data(pdf + 55);
    const auto *pdf_56 = buffer.data(pdf + 56);
    const auto *pdf_59 = buffer.data(pdf + 59);
    const auto *pdf_60 = buffer.data(pdf + 60);
    const auto *pdf_61 = buffer.data(pdf + 61);
    const auto *pdf_66 = buffer.data(pdf + 66);
    const auto *pdf_67 = buffer.data(pdf + 67);
    const auto *pdf_68 = buffer.data(pdf + 68);
    const auto *pdf_69 = buffer.data(pdf + 69);
    const auto *pdf_70 = buffer.data(pdf + 70);
    const auto *pdf_71 = buffer.data(pdf + 71);
    const auto *pdf_73 = buffer.data(pdf + 73);
    const auto *pdf_75 = buffer.data(pdf + 75);
    const auto *pdf_76 = buffer.data(pdf + 76);
    const auto *pdf_77 = buffer.data(pdf + 77);
    const auto *pdf_78 = buffer.data(pdf + 78);
    const auto *pdf_79 = buffer.data(pdf + 79);
    const auto *pdf_80 = buffer.data(pdf + 80);
    const auto *pdf_81 = buffer.data(pdf + 81);
    const auto *pdf_86 = buffer.data(pdf + 86);
    const auto *pdf_87 = buffer.data(pdf + 87);
    const auto *pdf_88 = buffer.data(pdf + 88);
    const auto *pdf_89 = buffer.data(pdf + 89);
    const auto *pdf_90 = buffer.data(pdf + 90);
    const auto *pdf_91 = buffer.data(pdf + 91);
    const auto *pdf_93 = buffer.data(pdf + 93);
    const auto *pdf_95 = buffer.data(pdf + 95);
    const auto *pdf_96 = buffer.data(pdf + 96);
    const auto *pdf_97 = buffer.data(pdf + 97);
    const auto *pdf_98 = buffer.data(pdf + 98);
    const auto *pdf_99 = buffer.data(pdf + 99);
    const auto *pdf_100 = buffer.data(pdf + 100);
    const auto *pdf_105 = buffer.data(pdf + 105);
    const auto *pdf_106 = buffer.data(pdf + 106);
    const auto *pdf_107 = buffer.data(pdf + 107);
    const auto *pdf_108 = buffer.data(pdf + 108);
    const auto *pdf_109 = buffer.data(pdf + 109);
    const auto *pdf_116 = buffer.data(pdf + 116);
    const auto *pdf_117 = buffer.data(pdf + 117);
    const auto *pdf_118 = buffer.data(pdf + 118);
    const auto *pdf_119 = buffer.data(pdf + 119);

    const auto *pdg1_80 = buffer.data(pdg1 + 80);
    const auto *pdg1_91 = buffer.data(pdg1 + 91);
    const auto *pdg1_93 = buffer.data(pdg1 + 93);
    const auto *pdg1_100 = buffer.data(pdg1 + 100);
    const auto *pdg1_106 = buffer.data(pdg1 + 106);
    const auto *pdg1_108 = buffer.data(pdg1 + 108);
    const auto *pdg1_115 = buffer.data(pdg1 + 115);
    const auto *pdg1_135 = buffer.data(pdg1 + 135);
    const auto *pdg1_136 = buffer.data(pdg1 + 136);
    const auto *pdg1_138 = buffer.data(pdg1 + 138);
    const auto *pdg1_140 = buffer.data(pdg1 + 140);
    const auto *pdg1_145 = buffer.data(pdg1 + 145);
    const auto *pdg1_147 = buffer.data(pdg1 + 147);
    const auto *pdg1_148 = buffer.data(pdg1 + 148);
    const auto *pdg1_149 = buffer.data(pdg1 + 149);
    const auto *pdg1_160 = buffer.data(pdg1 + 160);
    const auto *pdg1_162 = buffer.data(pdg1 + 162);
    const auto *pdg1_163 = buffer.data(pdg1 + 163);
    const auto *pdg1_164 = buffer.data(pdg1 + 164);
    const auto *pdg1_175 = buffer.data(pdg1 + 175);
    const auto *pdg1_177 = buffer.data(pdg1 + 177);

    const auto *pfd0_54 = buffer.data(pfd0 + 54);
    const auto *pfd0_63 = buffer.data(pfd0 + 63);
    const auto *pfd0_66 = buffer.data(pfd0 + 66);
    const auto *pfd0_67 = buffer.data(pfd0 + 67);
    const auto *pfd0_69 = buffer.data(pfd0 + 69);
    const auto *pfd0_71 = buffer.data(pfd0 + 71);
    const auto *pfd0_84 = buffer.data(pfd0 + 84);
    const auto *pfd0_89 = buffer.data(pfd0 + 89);
    const auto *pfd0_96 = buffer.data(pfd0 + 96);
    const auto *pfd0_97 = buffer.data(pfd0 + 97);
    const auto *pfd0_99 = buffer.data(pfd0 + 99);
    const auto *pfd0_101 = buffer.data(pfd0 + 101);

    const auto *pfd1_54 = buffer.data(pfd1 + 54);
    const auto *pfd1_63 = buffer.data(pfd1 + 63);
    const auto *pfd1_66 = buffer.data(pfd1 + 66);
    const auto *pfd1_67 = buffer.data(pfd1 + 67);
    const auto *pfd1_69 = buffer.data(pfd1 + 69);
    const auto *pfd1_71 = buffer.data(pfd1 + 71);
    const auto *pfd1_84 = buffer.data(pfd1 + 84);
    const auto *pfd1_89 = buffer.data(pfd1 + 89);
    const auto *pfd1_96 = buffer.data(pfd1 + 96);
    const auto *pfd1_97 = buffer.data(pfd1 + 97);
    const auto *pfd1_99 = buffer.data(pfd1 + 99);
    const auto *pfd1_101 = buffer.data(pfd1 + 101);

    const auto *pff_83 = buffer.data(pff + 83);
    const auto *pff_85 = buffer.data(pff + 85);
    const auto *pff_86 = buffer.data(pff + 86);
    const auto *pff_89 = buffer.data(pff + 89);
    const auto *pff_90 = buffer.data(pff + 90);
    const auto *pff_91 = buffer.data(pff + 91);
    const auto *pff_92 = buffer.data(pff + 92);
    const auto *pff_93 = buffer.data(pff + 93);
    const auto *pff_95 = buffer.data(pff + 95);
    const auto *pff_96 = buffer.data(pff + 96);
    const auto *pff_99 = buffer.data(pff + 99);
    const auto *pff_100 = buffer.data(pff + 100);
    const auto *pff_101 = buffer.data(pff + 101);
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
    const auto *pff_120 = buffer.data(pff + 120);
    const auto *pff_121 = buffer.data(pff + 121);
    const auto *pff_126 = buffer.data(pff + 126);
    const auto *pff_127 = buffer.data(pff + 127);
    const auto *pff_128 = buffer.data(pff + 128);
    const auto *pff_129 = buffer.data(pff + 129);
    const auto *pff_130 = buffer.data(pff + 130);
    const auto *pff_131 = buffer.data(pff + 131);
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

#pragma omp simd aligned(t_125, t_126, t_127, t_128, pb_y, pc_x, pc_y, pc_z, sff_86, pdg0_80, \
                         pdf_43, pdf_55, pdg1_80, pff_83, pff_85, \
                         pff_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_125[k] = pb_y[k] * pdg0_80[k]
                   - f_7 * pc_y[k] * pdg1_80[k];

        t_126[k] = f_0 * sff_86[k]
                   + f_4 * pc_x[k] * pff_86[k];

        t_127[k] = f_8 * pdf_43[k]
                   + f_4 * pc_z[k] * pff_83[k];

        t_128[k] = f_0 * pdf_55[k]
                   + f_4 * pc_y[k] * pff_85[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, pa_x, pc_x, pc_z, sfg0_130, sfg0_132, \
                         sff_89, sfg1_130, sfg1_132, pdf_46, pff_86, \
                         pff_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = f_0 * sff_89[k]
                   + f_4 * pc_x[k] * pff_89[k];

        t_130[k] = pa_x[k] * sfg0_130[k]
                   - f_7 * pc_x[k] * sfg1_130[k];

        t_131[k] = f_8 * pdf_46[k]
                   + f_4 * pc_z[k] * pff_86[k];

        t_132[k] = pa_x[k] * sfg0_132[k]
                   - f_7 * pc_x[k] * sfg1_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, t_136, pa_x, pc_x, pc_y, sfg0_134, sfg0_135, \
                         sff_90, sfg1_134, sfg1_135, pdf_59, pff_89, \
                         pff_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_0 * pdf_59[k]
                   + f_4 * pc_y[k] * pff_89[k];

        t_134[k] = pa_x[k] * sfg0_134[k]
                   - f_7 * pc_x[k] * sfg1_134[k];

        t_135[k] = pa_x[k] * sfg0_135[k]
                   + f_9 * sff_90[k]
                   - f_7 * pc_x[k] * sfg1_135[k];

        t_136[k] = f_4 * pc_y[k] * pff_90[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pc_y, pc_z, pdf_50, pfd0_54, pfd1_54, pff_90, \
                         pff_91, pff_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_1 * pdf_50[k]
                   + f_4 * pc_z[k] * pff_90[k];

        t_138[k] = f_5 * pfd0_54[k]
                   - f_6 * pfd1_54[k]
                   + f_4 * pc_y[k] * pff_91[k];

        t_139[k] = f_4 * pc_y[k] * pff_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pc_x, pc_y, pc_z, sfg0_140, sff_95, \
                         sff_96, sfg1_140, pdf_53, pff_93, pff_95, \
                         pff_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = pa_x[k] * sfg0_140[k]
                   + f_8 * sff_95[k]
                   - f_7 * pc_x[k] * sfg1_140[k];

        t_141[k] = f_0 * sff_96[k]
                   + f_4 * pc_x[k] * pff_96[k];

        t_142[k] = f_1 * pdf_53[k]
                   + f_4 * pc_z[k] * pff_93[k];

        t_143[k] = f_4 * pc_y[k] * pff_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pc_x, pc_z, sfg0_145, sfg0_147, \
                         sff_99, sfg1_145, sfg1_147, pdf_56, pff_96, \
                         pff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_0 * sff_99[k]
                   + f_4 * pc_x[k] * pff_99[k];

        t_145[k] = pa_x[k] * sfg0_145[k]
                   - f_7 * pc_x[k] * sfg1_145[k];

        t_146[k] = f_1 * pdf_56[k]
                   + f_4 * pc_z[k] * pff_96[k];

        t_147[k] = pa_x[k] * sfg0_147[k]
                   - f_7 * pc_x[k] * sfg1_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, t_151, pa_x, pa_y, pc_x, pc_y, sfg0_0, sfg0_1, \
                         sfg0_149, sff_0, sfg1_0, sfg1_1, sfg1_149, \
                         pff_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_4 * pc_y[k] * pff_99[k];

        t_149[k] = pa_x[k] * sfg0_149[k]
                   - f_7 * pc_x[k] * sfg1_149[k];

        t_150[k] = pa_y[k] * sfg0_0[k]
                   - f_7 * pc_y[k] * sfg1_0[k];

        t_151[k] = pa_y[k] * sfg0_1[k]
                   + f_0 * sff_0[k]
                   - f_7 * pc_y[k] * sfg1_1[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_y, pc_y, pc_z, sfg0_3, sfg0_5, sff_1, \
                         sfg1_3, sfg1_5, pff_100, pff_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_4 * pc_z[k] * pff_100[k];

        t_153[k] = pa_y[k] * sfg0_3[k]
                   + f_8 * sff_1[k]
                   - f_7 * pc_y[k] * sfg1_3[k];

        t_154[k] = f_4 * pc_z[k] * pff_101[k];

        t_155[k] = pa_y[k] * sfg0_5[k]
                   - f_7 * pc_y[k] * sfg1_5[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, t_159, pc_x, pdf_66, pdf_67, pdf_68, pdf_69, \
                         pff_106, pff_107, pff_108, pff_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_1 * pdf_66[k]
                   + f_4 * pc_x[k] * pff_106[k];

        t_157[k] = f_1 * pdf_67[k]
                   + f_4 * pc_x[k] * pff_107[k];

        t_158[k] = f_1 * pdf_68[k]
                   + f_4 * pc_x[k] * pff_108[k];

        t_159[k] = f_1 * pdf_69[k]
                   + f_4 * pc_x[k] * pff_109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_y, pc_y, pc_z, sfg0_10, sff_6, sff_9, \
                         sfg1_10, pfd0_63, pfd1_63, pff_106, pff_107, \
                         pff_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_y[k] * sfg0_10[k]
                   + f_9 * sff_6[k]
                   - f_7 * pc_y[k] * sfg1_10[k];

        t_161[k] = f_4 * pc_z[k] * pff_106[k];

        t_162[k] = f_5 * pfd0_63[k]
                   - f_6 * pfd1_63[k]
                   + f_4 * pc_z[k] * pff_107[k];

        t_163[k] = f_0 * sff_9[k]
                   + f_4 * pc_y[k] * pff_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pa_y, pc_x, pc_y, sfg0_14, sfg1_14, pdf_70, \
                         pdf_71, pfd0_66, pfd0_67, pfd1_66, pfd1_67, pff_110, \
                         pff_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_y[k] * sfg0_14[k]
                   - f_7 * pc_y[k] * sfg1_14[k];

        t_165[k] = f_8 * pdf_70[k]
                   + f_2 * pfd0_66[k]
                   - f_3 * pfd1_66[k]
                   + f_4 * pc_x[k] * pff_110[k];

        t_166[k] = f_8 * pdf_71[k]
                   + f_10 * pfd0_67[k]
                   - f_11 * pfd1_67[k]
                   + f_4 * pc_x[k] * pff_111[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, t_170, pc_x, pc_z, pdf_73, pdf_75, pfd0_69, \
                         pfd0_71, pfd1_69, pfd1_71, pff_110, pff_111, pff_113, \
                         pff_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_4 * pc_z[k] * pff_110[k];

        t_168[k] = f_8 * pdf_73[k]
                   + f_5 * pfd0_69[k]
                   - f_6 * pfd1_69[k]
                   + f_4 * pc_x[k] * pff_113[k];

        t_169[k] = f_4 * pc_z[k] * pff_111[k];

        t_170[k] = f_8 * pdf_75[k]
                   + f_5 * pfd0_71[k]
                   - f_6 * pfd1_71[k]
                   + f_4 * pc_x[k] * pff_115[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, t_174, pc_x, pdf_76, pdf_77, pdf_78, pdf_79, \
                         pff_116, pff_117, pff_118, pff_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_8 * pdf_76[k]
                   + f_4 * pc_x[k] * pff_116[k];

        t_172[k] = f_8 * pdf_77[k]
                   + f_4 * pc_x[k] * pff_117[k];

        t_173[k] = f_8 * pdf_78[k]
                   + f_4 * pc_x[k] * pff_118[k];

        t_174[k] = f_8 * pdf_79[k]
                   + f_4 * pc_x[k] * pff_119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pb_x, pc_x, pc_z, ppg0_70, ppg1_70, pdg0_115, \
                         pdg1_115, pfd0_69, pfd1_69, pff_116, pff_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_12 * ppg0_70[k]
                   - f_13 * ppg1_70[k]
                   + pb_x[k] * pdg0_115[k]
                   - f_7 * pc_x[k] * pdg1_115[k];

        t_176[k] = f_4 * pc_z[k] * pff_116[k];

        t_177[k] = f_5 * pfd0_69[k]
                   - f_6 * pfd1_69[k]
                   + f_4 * pc_z[k] * pff_117[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pa_y, pc_y, pc_z, sfg0_30, sff_19, sfg1_30, \
                         pdf_69, pfd0_71, pfd1_71, pff_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_0 * sff_19[k]
                   + f_0 * pdf_69[k]
                   + f_4 * pc_y[k] * pff_119[k];

        t_179[k] = f_2 * pfd0_71[k]
                   - f_3 * pfd1_71[k]
                   + f_4 * pc_z[k] * pff_119[k];

        t_180[k] = pa_y[k] * sfg0_30[k]
                   - f_7 * pc_y[k] * sfg1_30[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pb_z, pc_z, pdg0_91, pdg0_93, pdf_60, \
                         pdf_61, pdg1_91, pdg1_93, pff_120, pff_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = pb_z[k] * pdg0_91[k]
                   - f_7 * pc_z[k] * pdg1_91[k];

        t_182[k] = f_0 * pdf_60[k]
                   + f_4 * pc_z[k] * pff_120[k];

        t_183[k] = pb_z[k] * pdg0_93[k]
                   - f_7 * pc_z[k] * pdg1_93[k];

        t_184[k] = f_0 * pdf_61[k]
                   + f_4 * pc_z[k] * pff_121[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, t_188, pa_y, pc_x, pc_y, sfg0_35, sfg1_35, \
                         pdf_86, pdf_87, pdf_88, pff_126, pff_127, \
                         pff_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = pa_y[k] * sfg0_35[k]
                   - f_7 * pc_y[k] * sfg1_35[k];

        t_186[k] = f_8 * pdf_86[k]
                   + f_4 * pc_x[k] * pff_126[k];

        t_187[k] = f_8 * pdf_87[k]
                   + f_4 * pc_x[k] * pff_127[k];

        t_188[k] = f_8 * pdf_88[k]
                   + f_4 * pc_x[k] * pff_128[k];
    }

#pragma omp simd aligned(t_189, t_190, t_191, pb_z, pc_x, pc_z, pdg0_100, pdf_66, pdf_89, \
                         pdg1_100, pff_126, pff_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_189[k] = f_8 * pdf_89[k]
                   + f_4 * pc_x[k] * pff_129[k];

        t_190[k] = pb_z[k] * pdg0_100[k]
                   - f_7 * pc_z[k] * pdg1_100[k];

        t_191[k] = f_0 * pdf_66[k]
                   + f_4 * pc_z[k] * pff_126[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pc_y, sfg0_42, sfg0_44, sff_28, sff_29, \
                         sfg1_42, sfg1_44, pff_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = pa_y[k] * sfg0_42[k]
                   + f_8 * sff_28[k]
                   - f_7 * pc_y[k] * sfg1_42[k];

        t_193[k] = f_0 * sff_29[k]
                   + f_4 * pc_y[k] * pff_129[k];

        t_194[k] = pa_y[k] * sfg0_44[k]
                   - f_7 * pc_y[k] * sfg1_44[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pb_x, pc_x, pc_z, pdg0_135, pdg0_136, pdf_90, \
                         pdf_91, pdg1_135, pdg1_136, pff_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pb_x[k] * pdg0_135[k]
                   + f_9 * pdf_90[k]
                   - f_7 * pc_x[k] * pdg1_135[k];

        t_196[k] = pb_x[k] * pdg0_136[k]
                   + f_1 * pdf_91[k]
                   - f_7 * pc_x[k] * pdg1_136[k];

        t_197[k] = f_4 * pc_z[k] * pff_130[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pb_x, pc_x, pc_z, pdg0_138, pdg0_140, \
                         pdf_93, pdf_95, pdf_96, pdg1_138, pdg1_140, pff_131, \
                         pff_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = pb_x[k] * pdg0_138[k]
                   + f_8 * pdf_93[k]
                   - f_7 * pc_x[k] * pdg1_138[k];

        t_199[k] = f_4 * pc_z[k] * pff_131[k];

        t_200[k] = pb_x[k] * pdg0_140[k]
                   + f_8 * pdf_95[k]
                   - f_7 * pc_x[k] * pdg1_140[k];

        t_201[k] = f_0 * pdf_96[k]
                   + f_4 * pc_x[k] * pff_136[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pb_x, pc_x, pdg0_145, pdf_97, pdf_98, \
                         pdf_99, pdg1_145, pff_137, pff_138, pff_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_0 * pdf_97[k]
                   + f_4 * pc_x[k] * pff_137[k];

        t_203[k] = f_0 * pdf_98[k]
                   + f_4 * pc_x[k] * pff_138[k];

        t_204[k] = f_0 * pdf_99[k]
                   + f_4 * pc_x[k] * pff_139[k];

        t_205[k] = pb_x[k] * pdg0_145[k]
                   - f_7 * pc_x[k] * pdg1_145[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pb_x, pc_x, pc_z, pdg0_147, pdg0_148, \
                         pdg0_149, pdg1_147, pdg1_148, pdg1_149, \
                         pff_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_4 * pc_z[k] * pff_136[k];

        t_207[k] = pb_x[k] * pdg0_147[k]
                   - f_7 * pc_x[k] * pdg1_147[k];

        t_208[k] = pb_x[k] * pdg0_148[k]
                   - f_7 * pc_x[k] * pdg1_148[k];

        t_209[k] = pb_x[k] * pdg0_149[k]
                   - f_7 * pc_x[k] * pdg1_149[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pb_z, pc_x, pc_z, pdg0_106, pdg0_108, \
                         pdf_70, pdf_100, pdg1_106, pdg1_108, pfd0_84, pfd1_84, \
                         pff_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_0 * pdf_100[k]
                   + f_2 * pfd0_84[k]
                   - f_3 * pfd1_84[k]
                   + f_4 * pc_x[k] * pff_140[k];

        t_211[k] = pb_z[k] * pdg0_106[k]
                   - f_7 * pc_z[k] * pdg1_106[k];

        t_212[k] = f_0 * pdf_70[k]
                   + f_4 * pc_z[k] * pff_140[k];

        t_213[k] = pb_z[k] * pdg0_108[k]
                   - f_7 * pc_z[k] * pdg1_108[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pc_x, pc_z, pdf_71, pdf_105, pdf_106, \
                         pdf_107, pfd0_89, pfd1_89, pff_141, pff_145, pff_146, \
                         pff_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_0 * pdf_71[k]
                   + f_4 * pc_z[k] * pff_141[k];

        t_215[k] = f_0 * pdf_105[k]
                   + f_5 * pfd0_89[k]
                   - f_6 * pfd1_89[k]
                   + f_4 * pc_x[k] * pff_145[k];

        t_216[k] = f_0 * pdf_106[k]
                   + f_4 * pc_x[k] * pff_146[k];

        t_217[k] = f_0 * pdf_107[k]
                   + f_4 * pc_x[k] * pff_147[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pb_x, pc_x, pc_z, pdg0_160, pdf_76, \
                         pdf_108, pdf_109, pdg1_160, pff_146, pff_148, \
                         pff_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_0 * pdf_108[k]
                   + f_4 * pc_x[k] * pff_148[k];

        t_219[k] = f_0 * pdf_109[k]
                   + f_4 * pc_x[k] * pff_149[k];

        t_220[k] = pb_x[k] * pdg0_160[k]
                   - f_7 * pc_x[k] * pdg1_160[k];

        t_221[k] = f_0 * pdf_76[k]
                   + f_4 * pc_z[k] * pff_146[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, t_225, pa_y, pb_x, pc_x, pc_y, sfg0_75, sfg1_75, \
                         pdg0_162, pdg0_163, pdg0_164, pdg1_162, pdg1_163, \
                         pdg1_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = pb_x[k] * pdg0_162[k]
                   - f_7 * pc_x[k] * pdg1_162[k];

        t_223[k] = pb_x[k] * pdg0_163[k]
                   - f_7 * pc_x[k] * pdg1_163[k];

        t_224[k] = pb_x[k] * pdg0_164[k]
                   - f_7 * pc_x[k] * pdg1_164[k];

        t_225[k] = pa_y[k] * sfg0_75[k]
                   - f_7 * pc_y[k] * sfg1_75[k];
    }

#pragma omp simd aligned(t_226, t_227, t_228, pa_y, pc_y, pc_z, sfg0_76, sfg0_78, sff_50, \
                         sff_51, sfg1_76, sfg1_78, pdf_80, pff_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_226[k] = pa_y[k] * sfg0_76[k]
                   + f_0 * sff_50[k]
                   - f_7 * pc_y[k] * sfg1_76[k];

        t_227[k] = f_8 * pdf_80[k]
                   + f_4 * pc_z[k] * pff_150[k];

        t_228[k] = pa_y[k] * sfg0_78[k]
                   + f_8 * sff_51[k]
                   - f_7 * pc_y[k] * sfg1_78[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, pa_y, pc_x, pc_y, pc_z, sfg0_80, sfg1_80, \
                         pdf_81, pdf_116, pdf_117, pff_151, pff_156, \
                         pff_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_8 * pdf_81[k]
                   + f_4 * pc_z[k] * pff_151[k];

        t_230[k] = pa_y[k] * sfg0_80[k]
                   - f_7 * pc_y[k] * sfg1_80[k];

        t_231[k] = f_0 * pdf_116[k]
                   + f_4 * pc_x[k] * pff_156[k];

        t_232[k] = f_0 * pdf_117[k]
                   + f_4 * pc_x[k] * pff_157[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pb_x, pc_x, pc_z, pdg0_175, pdf_86, \
                         pdf_118, pdf_119, pdg1_175, pff_156, pff_158, \
                         pff_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_0 * pdf_118[k]
                   + f_4 * pc_x[k] * pff_158[k];

        t_234[k] = f_0 * pdf_119[k]
                   + f_4 * pc_x[k] * pff_159[k];

        t_235[k] = pb_x[k] * pdg0_175[k]
                   - f_7 * pc_x[k] * pdg1_175[k];

        t_236[k] = f_8 * pdf_86[k]
                   + f_4 * pc_z[k] * pff_156[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pb_x, pc_x, pc_y, sfg0_89, sff_59, \
                         sfg1_89, pdg0_177, pdg1_177, pff_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = pb_x[k] * pdg0_177[k]
                   - f_7 * pc_x[k] * pdg1_177[k];

        t_238[k] = f_0 * sff_59[k]
                   + f_4 * pc_y[k] * pff_159[k];

        t_239[k] = pa_y[k] * sfg0_89[k]
                   - f_7 * pc_y[k] * sfg1_89[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, t_243, t_244, pc_x, pc_z, pfd0_96, pfd0_97, \
                         pfd0_99, pfd1_96, pfd1_97, pfd1_99, pff_160, pff_161, \
                         pff_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = f_2 * pfd0_96[k]
                   - f_3 * pfd1_96[k]
                   + f_4 * pc_x[k] * pff_160[k];

        t_241[k] = f_10 * pfd0_97[k]
                   - f_11 * pfd1_97[k]
                   + f_4 * pc_x[k] * pff_161[k];

        t_242[k] = f_4 * pc_z[k] * pff_160[k];

        t_243[k] = f_5 * pfd0_99[k]
                   - f_6 * pfd1_99[k]
                   + f_4 * pc_x[k] * pff_163[k];

        t_244[k] = f_4 * pc_z[k] * pff_161[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, t_248, t_249, pc_x, pfd0_101, pfd1_101, pff_165, \
                         pff_166, pff_167, pff_168, pff_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_5 * pfd0_101[k]
                   - f_6 * pfd1_101[k]
                   + f_4 * pc_x[k] * pff_165[k];

        t_246[k] = f_4 * pc_x[k] * pff_166[k];

        t_247[k] = f_4 * pc_x[k] * pff_167[k];

        t_248[k] = f_4 * pc_x[k] * pff_168[k];

        t_249[k] = f_4 * pc_x[k] * pff_169[k];
    }
}

static auto
compute_prim_pfg_three_center_electron_repulsion_0_piece2(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfg0, const size_t sff,
                                                          const size_t sfg1, const size_t ppg0,
                                                          const size_t ppg1, const size_t pdg0,
                                                          const size_t pdf, const size_t pdg1,
                                                          const size_t pfd0, const size_t pfd1,
                                                          const size_t pff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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
    auto *t_370 = buffer.data(target + 370);
    auto *t_371 = buffer.data(target + 371);
    auto *t_372 = buffer.data(target + 372);
    auto *t_373 = buffer.data(target + 373);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *sfg0_0 = buffer.data(sfg0 + 0);
    const auto *sfg0_2 = buffer.data(sfg0 + 2);
    const auto *sfg0_3 = buffer.data(sfg0 + 3);
    const auto *sfg0_5 = buffer.data(sfg0 + 5);
    const auto *sfg0_10 = buffer.data(sfg0 + 10);
    const auto *sfg0_14 = buffer.data(sfg0 + 14);
    const auto *sfg0_15 = buffer.data(sfg0 + 15);
    const auto *sfg0_18 = buffer.data(sfg0 + 18);
    const auto *sfg0_25 = buffer.data(sfg0 + 25);
    const auto *sfg0_26 = buffer.data(sfg0 + 26);
    const auto *sfg0_27 = buffer.data(sfg0 + 27);
    const auto *sfg0_45 = buffer.data(sfg0 + 45);
    const auto *sfg0_47 = buffer.data(sfg0 + 47);
    const auto *sfg0_48 = buffer.data(sfg0 + 48);
    const auto *sfg0_50 = buffer.data(sfg0 + 50);
    const auto *sfg0_55 = buffer.data(sfg0 + 55);
    const auto *sfg0_135 = buffer.data(sfg0 + 135);
    const auto *sfg0_140 = buffer.data(sfg0 + 140);
    const auto *sfg0_145 = buffer.data(sfg0 + 145);
    const auto *sfg0_147 = buffer.data(sfg0 + 147);
    const auto *sfg0_149 = buffer.data(sfg0 + 149);

    const auto *sff_0 = buffer.data(sff + 0);
    const auto *sff_2 = buffer.data(sff + 2);
    const auto *sff_9 = buffer.data(sff + 9);
    const auto *sff_16 = buffer.data(sff + 16);
    const auto *sff_17 = buffer.data(sff + 17);
    const auto *sff_30 = buffer.data(sff + 30);
    const auto *sff_32 = buffer.data(sff + 32);
    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_79 = buffer.data(sff + 79);
    const auto *sff_86 = buffer.data(sff + 86);
    const auto *sff_89 = buffer.data(sff + 89);
    const auto *sff_96 = buffer.data(sff + 96);
    const auto *sff_98 = buffer.data(sff + 98);
    const auto *sff_99 = buffer.data(sff + 99);

    const auto *sfg1_0 = buffer.data(sfg1 + 0);
    const auto *sfg1_2 = buffer.data(sfg1 + 2);
    const auto *sfg1_3 = buffer.data(sfg1 + 3);
    const auto *sfg1_5 = buffer.data(sfg1 + 5);
    const auto *sfg1_10 = buffer.data(sfg1 + 10);
    const auto *sfg1_14 = buffer.data(sfg1 + 14);
    const auto *sfg1_15 = buffer.data(sfg1 + 15);
    const auto *sfg1_18 = buffer.data(sfg1 + 18);
    const auto *sfg1_25 = buffer.data(sfg1 + 25);
    const auto *sfg1_26 = buffer.data(sfg1 + 26);
    const auto *sfg1_27 = buffer.data(sfg1 + 27);
    const auto *sfg1_45 = buffer.data(sfg1 + 45);
    const auto *sfg1_47 = buffer.data(sfg1 + 47);
    const auto *sfg1_48 = buffer.data(sfg1 + 48);
    const auto *sfg1_50 = buffer.data(sfg1 + 50);
    const auto *sfg1_55 = buffer.data(sfg1 + 55);
    const auto *sfg1_135 = buffer.data(sfg1 + 135);
    const auto *sfg1_140 = buffer.data(sfg1 + 140);
    const auto *sfg1_145 = buffer.data(sfg1 + 145);
    const auto *sfg1_147 = buffer.data(sfg1 + 147);
    const auto *sfg1_149 = buffer.data(sfg1 + 149);

    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *pdg0_135 = buffer.data(pdg0 + 135);
    const auto *pdg0_136 = buffer.data(pdg0 + 136);
    const auto *pdg0_138 = buffer.data(pdg0 + 138);
    const auto *pdg0_145 = buffer.data(pdg0 + 145);
    const auto *pdg0_147 = buffer.data(pdg0 + 147);
    const auto *pdg0_182 = buffer.data(pdg0 + 182);
    const auto *pdg0_185 = buffer.data(pdg0 + 185);
    const auto *pdg0_194 = buffer.data(pdg0 + 194);
    const auto *pdg0_210 = buffer.data(pdg0 + 210);
    const auto *pdg0_212 = buffer.data(pdg0 + 212);
    const auto *pdg0_215 = buffer.data(pdg0 + 215);
    const auto *pdg0_224 = buffer.data(pdg0 + 224);
    const auto *pdg0_236 = buffer.data(pdg0 + 236);
    const auto *pdg0_237 = buffer.data(pdg0 + 237);
    const auto *pdg0_239 = buffer.data(pdg0 + 239);
    const auto *pdg0_250 = buffer.data(pdg0 + 250);
    const auto *pdg0_251 = buffer.data(pdg0 + 251);
    const auto *pdg0_252 = buffer.data(pdg0 + 252);

    const auto *pdf_90 = buffer.data(pdf + 90);
    const auto *pdf_91 = buffer.data(pdf + 91);
    const auto *pdf_96 = buffer.data(pdf + 96);
    const auto *pdf_97 = buffer.data(pdf + 97);
    const auto *pdf_99 = buffer.data(pdf + 99);
    const auto *pdf_100 = buffer.data(pdf + 100);
    const auto *pdf_101 = buffer.data(pdf + 101);
    const auto *pdf_106 = buffer.data(pdf + 106);
    const auto *pdf_107 = buffer.data(pdf + 107);
    const auto *pdf_109 = buffer.data(pdf + 109);
    const auto *pdf_110 = buffer.data(pdf + 110);
    const auto *pdf_111 = buffer.data(pdf + 111);
    const auto *pdf_116 = buffer.data(pdf + 116);
    const auto *pdf_119 = buffer.data(pdf + 119);
    const auto *pdf_120 = buffer.data(pdf + 120);
    const auto *pdf_122 = buffer.data(pdf + 122);
    const auto *pdf_126 = buffer.data(pdf + 126);
    const auto *pdf_127 = buffer.data(pdf + 127);
    const auto *pdf_128 = buffer.data(pdf + 128);
    const auto *pdf_129 = buffer.data(pdf + 129);
    const auto *pdf_130 = buffer.data(pdf + 130);
    const auto *pdf_132 = buffer.data(pdf + 132);
    const auto *pdf_136 = buffer.data(pdf + 136);
    const auto *pdf_137 = buffer.data(pdf + 137);
    const auto *pdf_138 = buffer.data(pdf + 138);
    const auto *pdf_139 = buffer.data(pdf + 139);
    const auto *pdf_140 = buffer.data(pdf + 140);
    const auto *pdf_142 = buffer.data(pdf + 142);
    const auto *pdf_143 = buffer.data(pdf + 143);
    const auto *pdf_145 = buffer.data(pdf + 145);
    const auto *pdf_146 = buffer.data(pdf + 146);
    const auto *pdf_147 = buffer.data(pdf + 147);
    const auto *pdf_148 = buffer.data(pdf + 148);
    const auto *pdf_149 = buffer.data(pdf + 149);
    const auto *pdf_156 = buffer.data(pdf + 156);
    const auto *pdf_157 = buffer.data(pdf + 157);
    const auto *pdf_158 = buffer.data(pdf + 158);
    const auto *pdf_159 = buffer.data(pdf + 159);
    const auto *pdf_163 = buffer.data(pdf + 163);
    const auto *pdf_166 = buffer.data(pdf + 166);
    const auto *pdf_167 = buffer.data(pdf + 167);
    const auto *pdf_168 = buffer.data(pdf + 168);
    const auto *pdf_169 = buffer.data(pdf + 169);

    const auto *pdg1_135 = buffer.data(pdg1 + 135);
    const auto *pdg1_136 = buffer.data(pdg1 + 136);
    const auto *pdg1_138 = buffer.data(pdg1 + 138);
    const auto *pdg1_145 = buffer.data(pdg1 + 145);
    const auto *pdg1_147 = buffer.data(pdg1 + 147);
    const auto *pdg1_182 = buffer.data(pdg1 + 182);
    const auto *pdg1_185 = buffer.data(pdg1 + 185);
    const auto *pdg1_194 = buffer.data(pdg1 + 194);
    const auto *pdg1_210 = buffer.data(pdg1 + 210);
    const auto *pdg1_212 = buffer.data(pdg1 + 212);
    const auto *pdg1_215 = buffer.data(pdg1 + 215);
    const auto *pdg1_224 = buffer.data(pdg1 + 224);
    const auto *pdg1_236 = buffer.data(pdg1 + 236);
    const auto *pdg1_237 = buffer.data(pdg1 + 237);
    const auto *pdg1_239 = buffer.data(pdg1 + 239);
    const auto *pdg1_250 = buffer.data(pdg1 + 250);
    const auto *pdg1_251 = buffer.data(pdg1 + 251);
    const auto *pdg1_252 = buffer.data(pdg1 + 252);

    const auto *pfd0_99 = buffer.data(pfd0 + 99);
    const auto *pfd0_101 = buffer.data(pfd0 + 101);
    const auto *pfd0_107 = buffer.data(pfd0 + 107);
    const auto *pfd0_108 = buffer.data(pfd0 + 108);
    const auto *pfd0_109 = buffer.data(pfd0 + 109);
    const auto *pfd0_111 = buffer.data(pfd0 + 111);
    const auto *pfd0_113 = buffer.data(pfd0 + 113);
    const auto *pfd0_115 = buffer.data(pfd0 + 115);
    const auto *pfd0_117 = buffer.data(pfd0 + 117);
    const auto *pfd0_124 = buffer.data(pfd0 + 124);
    const auto *pfd0_125 = buffer.data(pfd0 + 125);
    const auto *pfd0_132 = buffer.data(pfd0 + 132);
    const auto *pfd0_134 = buffer.data(pfd0 + 134);
    const auto *pfd0_135 = buffer.data(pfd0 + 135);
    const auto *pfd0_136 = buffer.data(pfd0 + 136);
    const auto *pfd0_137 = buffer.data(pfd0 + 137);
    const auto *pfd0_147 = buffer.data(pfd0 + 147);

    const auto *pfd1_99 = buffer.data(pfd1 + 99);
    const auto *pfd1_101 = buffer.data(pfd1 + 101);
    const auto *pfd1_107 = buffer.data(pfd1 + 107);
    const auto *pfd1_108 = buffer.data(pfd1 + 108);
    const auto *pfd1_109 = buffer.data(pfd1 + 109);
    const auto *pfd1_111 = buffer.data(pfd1 + 111);
    const auto *pfd1_113 = buffer.data(pfd1 + 113);
    const auto *pfd1_115 = buffer.data(pfd1 + 115);
    const auto *pfd1_117 = buffer.data(pfd1 + 117);
    const auto *pfd1_124 = buffer.data(pfd1 + 124);
    const auto *pfd1_125 = buffer.data(pfd1 + 125);
    const auto *pfd1_132 = buffer.data(pfd1 + 132);
    const auto *pfd1_134 = buffer.data(pfd1 + 134);
    const auto *pfd1_135 = buffer.data(pfd1 + 135);
    const auto *pfd1_136 = buffer.data(pfd1 + 136);
    const auto *pfd1_137 = buffer.data(pfd1 + 137);
    const auto *pfd1_147 = buffer.data(pfd1 + 147);

    const auto *pff_166 = buffer.data(pff + 166);
    const auto *pff_167 = buffer.data(pff + 167);
    const auto *pff_169 = buffer.data(pff + 169);
    const auto *pff_170 = buffer.data(pff + 170);
    const auto *pff_171 = buffer.data(pff + 171);
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
    const auto *pff_190 = buffer.data(pff + 190);
    const auto *pff_191 = buffer.data(pff + 191);
    const auto *pff_193 = buffer.data(pff + 193);
    const auto *pff_196 = buffer.data(pff + 196);
    const auto *pff_197 = buffer.data(pff + 197);
    const auto *pff_198 = buffer.data(pff + 198);
    const auto *pff_199 = buffer.data(pff + 199);
    const auto *pff_200 = buffer.data(pff + 200);
    const auto *pff_202 = buffer.data(pff + 202);
    const auto *pff_206 = buffer.data(pff + 206);
    const auto *pff_207 = buffer.data(pff + 207);
    const auto *pff_208 = buffer.data(pff + 208);
    const auto *pff_209 = buffer.data(pff + 209);
    const auto *pff_210 = buffer.data(pff + 210);
    const auto *pff_212 = buffer.data(pff + 212);
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

#pragma omp simd aligned(t_250, t_251, t_252, t_253, pc_y, pc_z, sff_66, sff_69, pdf_96, \
                         pdf_99, pfd0_99, pfd1_99, pff_166, pff_167, \
                         pff_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = f_0 * sff_66[k]
                   + f_1 * pdf_96[k]
                   + f_2 * pfd0_99[k]
                   - f_3 * pfd1_99[k]
                   + f_4 * pc_y[k] * pff_166[k];

        t_251[k] = f_4 * pc_z[k] * pff_166[k];

        t_252[k] = f_5 * pfd0_99[k]
                   - f_6 * pfd1_99[k]
                   + f_4 * pc_z[k] * pff_167[k];

        t_253[k] = f_0 * sff_69[k]
                   + f_1 * pdf_99[k]
                   + f_4 * pc_y[k] * pff_169[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, pb_z, pc_z, pdg0_135, pdg0_136, pdf_90, \
                         pdg1_135, pdg1_136, pfd0_101, pfd1_101, pff_169, \
                         pff_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_2 * pfd0_101[k]
                   - f_3 * pfd1_101[k]
                   + f_4 * pc_z[k] * pff_169[k];

        t_255[k] = pb_z[k] * pdg0_135[k]
                   - f_7 * pc_z[k] * pdg1_135[k];

        t_256[k] = pb_z[k] * pdg0_136[k]
                   - f_7 * pc_z[k] * pdg1_136[k];

        t_257[k] = f_0 * pdf_90[k]
                   + f_4 * pc_z[k] * pff_170[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, t_261, pb_z, pc_x, pc_z, pdg0_138, pdf_91, \
                         pdg1_138, pfd0_107, pfd1_107, pff_171, pff_175, \
                         pff_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = pb_z[k] * pdg0_138[k]
                   - f_7 * pc_z[k] * pdg1_138[k];

        t_259[k] = f_0 * pdf_91[k]
                   + f_4 * pc_z[k] * pff_171[k];

        t_260[k] = f_5 * pfd0_107[k]
                   - f_6 * pfd1_107[k]
                   + f_4 * pc_x[k] * pff_175[k];

        t_261[k] = f_4 * pc_x[k] * pff_176[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, t_266, pb_z, pc_x, pc_z, pdg0_145, \
                         pdf_96, pdg1_145, pff_176, pff_177, pff_178, \
                         pff_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_4 * pc_x[k] * pff_177[k];

        t_263[k] = f_4 * pc_x[k] * pff_178[k];

        t_264[k] = f_4 * pc_x[k] * pff_179[k];

        t_265[k] = pb_z[k] * pdg0_145[k]
                   - f_7 * pc_z[k] * pdg1_145[k];

        t_266[k] = f_0 * pdf_96[k]
                   + f_4 * pc_z[k] * pff_176[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pb_z, pc_y, pc_z, sff_79, pdg0_147, pdf_97, \
                         pdf_99, pdf_109, pdg1_147, pfd0_107, pfd1_107, \
                         pff_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pb_z[k] * pdg0_147[k]
                   + f_8 * pdf_97[k]
                   - f_7 * pc_z[k] * pdg1_147[k];

        t_268[k] = f_0 * sff_79[k]
                   + f_8 * pdf_109[k]
                   + f_4 * pc_y[k] * pff_179[k];

        t_269[k] = f_0 * pdf_99[k]
                   + f_2 * pfd0_107[k]
                   - f_3 * pfd1_107[k]
                   + f_4 * pc_z[k] * pff_179[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, pc_x, pc_z, pdf_100, pfd0_108, pfd0_109, \
                         pfd0_111, pfd1_108, pfd1_109, pfd1_111, pff_180, pff_181, \
                         pff_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_2 * pfd0_108[k]
                   - f_3 * pfd1_108[k]
                   + f_4 * pc_x[k] * pff_180[k];

        t_271[k] = f_10 * pfd0_109[k]
                   - f_11 * pfd1_109[k]
                   + f_4 * pc_x[k] * pff_181[k];

        t_272[k] = f_8 * pdf_100[k]
                   + f_4 * pc_z[k] * pff_180[k];

        t_273[k] = f_5 * pfd0_111[k]
                   - f_6 * pfd1_111[k]
                   + f_4 * pc_x[k] * pff_183[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, t_277, t_278, pc_x, pc_z, pdf_101, pfd0_113, \
                         pfd1_113, pff_181, pff_185, pff_186, pff_187, \
                         pff_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = f_8 * pdf_101[k]
                   + f_4 * pc_z[k] * pff_181[k];

        t_275[k] = f_5 * pfd0_113[k]
                   - f_6 * pfd1_113[k]
                   + f_4 * pc_x[k] * pff_185[k];

        t_276[k] = f_4 * pc_x[k] * pff_186[k];

        t_277[k] = f_4 * pc_x[k] * pff_187[k];

        t_278[k] = f_4 * pc_x[k] * pff_188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pc_x, pc_y, pc_z, sff_86, pdf_106, \
                         pdf_107, pdf_116, pfd0_111, pfd1_111, pff_186, pff_187, \
                         pff_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_4 * pc_x[k] * pff_189[k];

        t_280[k] = f_0 * sff_86[k]
                   + f_0 * pdf_116[k]
                   + f_2 * pfd0_111[k]
                   - f_3 * pfd1_111[k]
                   + f_4 * pc_y[k] * pff_186[k];

        t_281[k] = f_8 * pdf_106[k]
                   + f_4 * pc_z[k] * pff_186[k];

        t_282[k] = f_8 * pdf_107[k]
                   + f_5 * pfd0_111[k]
                   - f_6 * pfd1_111[k]
                   + f_4 * pc_z[k] * pff_187[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, pa_y, pc_y, pc_z, sfg0_135, sff_89, sfg1_135, \
                         pdf_109, pdf_119, pfd0_113, pfd1_113, \
                         pff_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = f_0 * sff_89[k]
                   + f_0 * pdf_119[k]
                   + f_4 * pc_y[k] * pff_189[k];

        t_284[k] = f_8 * pdf_109[k]
                   + f_2 * pfd0_113[k]
                   - f_3 * pfd1_113[k]
                   + f_4 * pc_z[k] * pff_189[k];

        t_285[k] = pa_y[k] * sfg0_135[k]
                   - f_7 * pc_y[k] * sfg1_135[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pc_x, pc_z, pdf_110, pdf_111, pfd0_115, \
                         pfd0_117, pfd1_115, pfd1_117, pff_190, pff_191, \
                         pff_193 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_10 * pfd0_115[k]
                   - f_11 * pfd1_115[k]
                   + f_4 * pc_x[k] * pff_191[k];

        t_287[k] = f_1 * pdf_110[k]
                   + f_4 * pc_z[k] * pff_190[k];

        t_288[k] = f_5 * pfd0_117[k]
                   - f_6 * pfd1_117[k]
                   + f_4 * pc_x[k] * pff_193[k];

        t_289[k] = f_1 * pdf_111[k]
                   + f_4 * pc_z[k] * pff_191[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, t_294, pa_y, pc_x, pc_y, sfg0_140, \
                         sfg1_140, pff_196, pff_197, pff_198, pff_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * sfg0_140[k]
                   - f_7 * pc_y[k] * sfg1_140[k];

        t_291[k] = f_4 * pc_x[k] * pff_196[k];

        t_292[k] = f_4 * pc_x[k] * pff_197[k];

        t_293[k] = f_4 * pc_x[k] * pff_198[k];

        t_294[k] = f_4 * pc_x[k] * pff_199[k];
    }

#pragma omp simd aligned(t_295, t_296, t_297, pa_y, pc_y, pc_z, sfg0_145, sfg0_147, sff_96, \
                         sff_98, sfg1_145, sfg1_147, pdf_116, pff_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_295[k] = pa_y[k] * sfg0_145[k]
                   + f_9 * sff_96[k]
                   - f_7 * pc_y[k] * sfg1_145[k];

        t_296[k] = f_1 * pdf_116[k]
                   + f_4 * pc_z[k] * pff_196[k];

        t_297[k] = pa_y[k] * sfg0_147[k]
                   + f_8 * sff_98[k]
                   - f_7 * pc_y[k] * sfg1_147[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pa_z, pc_y, pc_z, sfg0_0, sfg0_149, \
                         sff_99, sfg1_0, sfg1_149, pff_199, pff_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_0 * sff_99[k]
                   + f_4 * pc_y[k] * pff_199[k];

        t_299[k] = pa_y[k] * sfg0_149[k]
                   - f_7 * pc_y[k] * sfg1_149[k];

        t_300[k] = pa_z[k] * sfg0_0[k]
                   - f_7 * pc_z[k] * sfg1_0[k];

        t_301[k] = f_4 * pc_y[k] * pff_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pa_z, pc_y, pc_z, sfg0_2, sfg0_3, sfg0_5, \
                         sff_0, sff_2, sfg1_2, sfg1_3, sfg1_5, \
                         pff_202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pa_z[k] * sfg0_2[k]
                   + f_0 * sff_0[k]
                   - f_7 * pc_z[k] * sfg1_2[k];

        t_303[k] = pa_z[k] * sfg0_3[k]
                   - f_7 * pc_z[k] * sfg1_3[k];

        t_304[k] = f_4 * pc_y[k] * pff_202[k];

        t_305[k] = pa_z[k] * sfg0_5[k]
                   + f_8 * sff_2[k]
                   - f_7 * pc_z[k] * sfg1_5[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pc_x, pdf_126, pdf_127, pdf_128, pdf_129, \
                         pff_206, pff_207, pff_208, pff_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_1 * pdf_126[k]
                   + f_4 * pc_x[k] * pff_206[k];

        t_307[k] = f_1 * pdf_127[k]
                   + f_4 * pc_x[k] * pff_207[k];

        t_308[k] = f_1 * pdf_128[k]
                   + f_4 * pc_x[k] * pff_208[k];

        t_309[k] = f_1 * pdf_129[k]
                   + f_4 * pc_x[k] * pff_209[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pa_z, pc_y, pc_z, sfg0_10, sfg1_10, \
                         pfd0_124, pfd0_125, pfd1_124, pfd1_125, pff_207, pff_208, \
                         pff_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = pa_z[k] * sfg0_10[k]
                   - f_7 * pc_z[k] * sfg1_10[k];

        t_311[k] = f_10 * pfd0_124[k]
                   - f_11 * pfd1_124[k]
                   + f_4 * pc_y[k] * pff_207[k];

        t_312[k] = f_5 * pfd0_125[k]
                   - f_6 * pfd1_125[k]
                   + f_4 * pc_y[k] * pff_208[k];

        t_313[k] = f_4 * pc_y[k] * pff_209[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, pa_z, pc_y, pc_z, sfg0_14, sfg0_15, sff_9, \
                         sfg1_14, sfg1_15, pdf_120, pff_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pa_z[k] * sfg0_14[k]
                   + f_9 * sff_9[k]
                   - f_7 * pc_z[k] * sfg1_14[k];

        t_315[k] = pa_z[k] * sfg0_15[k]
                   - f_7 * pc_z[k] * sfg1_15[k];

        t_316[k] = f_0 * pdf_120[k]
                   + f_4 * pc_y[k] * pff_210[k];
    }

#pragma omp simd aligned(t_317, t_318, t_319, t_320, pa_z, pb_y, pc_y, pc_z, sfg0_18, sfg1_18, \
                         pdg0_182, pdg0_185, pdf_122, pdg1_182, pdg1_185, \
                         pff_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_317[k] = pb_y[k] * pdg0_182[k]
                   - f_7 * pc_y[k] * pdg1_182[k];

        t_318[k] = pa_z[k] * sfg0_18[k]
                   - f_7 * pc_z[k] * sfg1_18[k];

        t_319[k] = f_0 * pdf_122[k]
                   + f_4 * pc_y[k] * pff_212[k];

        t_320[k] = pb_y[k] * pdg0_185[k]
                   - f_7 * pc_y[k] * pdg1_185[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pc_x, pdf_136, pdf_137, pdf_138, pdf_139, \
                         pff_216, pff_217, pff_218, pff_219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_8 * pdf_136[k]
                   + f_4 * pc_x[k] * pff_216[k];

        t_322[k] = f_8 * pdf_137[k]
                   + f_4 * pc_x[k] * pff_217[k];

        t_323[k] = f_8 * pdf_138[k]
                   + f_4 * pc_x[k] * pff_218[k];

        t_324[k] = f_8 * pdf_139[k]
                   + f_4 * pc_x[k] * pff_219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pa_z, pc_z, sfg0_25, sfg0_26, sfg0_27, sff_16, \
                         sff_17, sfg1_25, sfg1_26, sfg1_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = pa_z[k] * sfg0_25[k]
                   - f_7 * pc_z[k] * sfg1_25[k];

        t_326[k] = pa_z[k] * sfg0_26[k]
                   + f_0 * sff_16[k]
                   - f_7 * pc_z[k] * sfg1_26[k];

        t_327[k] = pa_z[k] * sfg0_27[k]
                   + f_8 * sff_17[k]
                   - f_7 * pc_z[k] * sfg1_27[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, t_331, pb_y, pc_x, pc_y, pdg0_194, pdf_129, \
                         pdf_140, pdg1_194, pfd0_132, pfd1_132, pff_219, \
                         pff_220 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * pdf_129[k]
                   + f_4 * pc_y[k] * pff_219[k];

        t_329[k] = pb_y[k] * pdg0_194[k]
                   - f_7 * pc_y[k] * pdg1_194[k];

        t_330[k] = f_8 * pdf_140[k]
                   + f_2 * pfd0_132[k]
                   - f_3 * pfd1_132[k]
                   + f_4 * pc_x[k] * pff_220[k];

        t_331[k] = f_4 * pc_y[k] * pff_220[k];
    }

#pragma omp simd aligned(t_332, t_333, t_334, pc_x, pc_y, pdf_142, pdf_143, pfd0_134, \
                         pfd0_135, pfd1_134, pfd1_135, pff_222, \
                         pff_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_332[k] = f_8 * pdf_142[k]
                   + f_10 * pfd0_134[k]
                   - f_11 * pfd1_134[k]
                   + f_4 * pc_x[k] * pff_222[k];

        t_333[k] = f_8 * pdf_143[k]
                   + f_5 * pfd0_135[k]
                   - f_6 * pfd1_135[k]
                   + f_4 * pc_x[k] * pff_223[k];

        t_334[k] = f_4 * pc_y[k] * pff_222[k];
    }

#pragma omp simd aligned(t_335, t_336, t_337, t_338, pc_x, pdf_145, pdf_146, pdf_147, pdf_148, \
                         pfd0_137, pfd1_137, pff_225, pff_226, pff_227, \
                         pff_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_335[k] = f_8 * pdf_145[k]
                   + f_5 * pfd0_137[k]
                   - f_6 * pfd1_137[k]
                   + f_4 * pc_x[k] * pff_225[k];

        t_336[k] = f_8 * pdf_146[k]
                   + f_4 * pc_x[k] * pff_226[k];

        t_337[k] = f_8 * pdf_147[k]
                   + f_4 * pc_x[k] * pff_227[k];

        t_338[k] = f_8 * pdf_148[k]
                   + f_4 * pc_x[k] * pff_228[k];
    }

#pragma omp simd aligned(t_339, t_340, t_341, pc_x, pc_y, pdf_149, pfd0_135, pfd0_136, \
                         pfd1_135, pfd1_136, pff_226, pff_227, \
                         pff_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_339[k] = f_8 * pdf_149[k]
                   + f_4 * pc_x[k] * pff_229[k];

        t_340[k] = f_2 * pfd0_135[k]
                   - f_3 * pfd1_135[k]
                   + f_4 * pc_y[k] * pff_226[k];

        t_341[k] = f_10 * pfd0_136[k]
                   - f_11 * pfd1_136[k]
                   + f_4 * pc_y[k] * pff_227[k];
    }

#pragma omp simd aligned(t_342, t_343, t_344, pb_x, pc_x, pc_y, ppg0_134, ppg1_134, pdg0_224, \
                         pdg1_224, pfd0_137, pfd1_137, pff_228, \
                         pff_229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_342[k] = f_5 * pfd0_137[k]
                   - f_6 * pfd1_137[k]
                   + f_4 * pc_y[k] * pff_228[k];

        t_343[k] = f_4 * pc_y[k] * pff_229[k];

        t_344[k] = f_12 * ppg0_134[k]
                   - f_13 * ppg1_134[k]
                   + pb_x[k] * pdg0_224[k]
                   - f_7 * pc_x[k] * pdg1_224[k];
    }

#pragma omp simd aligned(t_345, t_346, t_347, t_348, pa_z, pc_y, pc_z, sfg0_45, sfg0_47, \
                         sfg0_48, sff_30, sfg1_45, sfg1_47, sfg1_48, pdf_130, \
                         pff_230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_345[k] = pa_z[k] * sfg0_45[k]
                   - f_7 * pc_z[k] * sfg1_45[k];

        t_346[k] = f_8 * pdf_130[k]
                   + f_4 * pc_y[k] * pff_230[k];

        t_347[k] = pa_z[k] * sfg0_47[k]
                   + f_0 * sff_30[k]
                   - f_7 * pc_z[k] * sfg1_47[k];

        t_348[k] = pa_z[k] * sfg0_48[k]
                   - f_7 * pc_z[k] * sfg1_48[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, pa_z, pc_x, pc_y, pc_z, sfg0_50, sff_32, \
                         sfg1_50, pdf_132, pdf_156, pff_232, pff_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = f_8 * pdf_132[k]
                   + f_4 * pc_y[k] * pff_232[k];

        t_350[k] = pa_z[k] * sfg0_50[k]
                   + f_8 * sff_32[k]
                   - f_7 * pc_z[k] * sfg1_50[k];

        t_351[k] = f_0 * pdf_156[k]
                   + f_4 * pc_x[k] * pff_236[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pc_x, pc_z, sfg0_55, sfg1_55, \
                         pdf_157, pdf_158, pdf_159, pff_237, pff_238, \
                         pff_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_0 * pdf_157[k]
                   + f_4 * pc_x[k] * pff_237[k];

        t_353[k] = f_0 * pdf_158[k]
                   + f_4 * pc_x[k] * pff_238[k];

        t_354[k] = f_0 * pdf_159[k]
                   + f_4 * pc_x[k] * pff_239[k];

        t_355[k] = pa_z[k] * sfg0_55[k]
                   - f_7 * pc_z[k] * sfg1_55[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, t_359, pb_x, pc_x, pc_y, pdg0_236, pdg0_237, \
                         pdg0_239, pdf_139, pdg1_236, pdg1_237, pdg1_239, \
                         pff_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = pb_x[k] * pdg0_236[k]
                   - f_7 * pc_x[k] * pdg1_236[k];

        t_357[k] = pb_x[k] * pdg0_237[k]
                   - f_7 * pc_x[k] * pdg1_237[k];

        t_358[k] = f_8 * pdf_139[k]
                   + f_4 * pc_y[k] * pff_239[k];

        t_359[k] = pb_x[k] * pdg0_239[k]
                   - f_7 * pc_x[k] * pdg1_239[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pb_y, pc_y, pdg0_210, pdg0_212, pdf_140, \
                         pdg1_210, pdg1_212, pff_240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = pb_y[k] * pdg0_210[k]
                   - f_7 * pc_y[k] * pdg1_210[k];

        t_361[k] = f_0 * pdf_140[k]
                   + f_4 * pc_y[k] * pff_240[k];

        t_362[k] = pb_y[k] * pdg0_212[k]
                   - f_7 * pc_y[k] * pdg1_212[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pb_y, pc_x, pc_y, pdg0_215, pdf_142, pdf_163, \
                         pdg1_215, pfd0_147, pfd1_147, pff_242, \
                         pff_243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = f_0 * pdf_163[k]
                   + f_5 * pfd0_147[k]
                   - f_6 * pfd1_147[k]
                   + f_4 * pc_x[k] * pff_243[k];

        t_364[k] = f_0 * pdf_142[k]
                   + f_4 * pc_y[k] * pff_242[k];

        t_365[k] = pb_y[k] * pdg0_215[k]
                   - f_7 * pc_y[k] * pdg1_215[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pc_x, pdf_166, pdf_167, pdf_168, pdf_169, \
                         pff_246, pff_247, pff_248, pff_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = f_0 * pdf_166[k]
                   + f_4 * pc_x[k] * pff_246[k];

        t_367[k] = f_0 * pdf_167[k]
                   + f_4 * pc_x[k] * pff_247[k];

        t_368[k] = f_0 * pdf_168[k]
                   + f_4 * pc_x[k] * pff_248[k];

        t_369[k] = f_0 * pdf_169[k]
                   + f_4 * pc_x[k] * pff_249[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, t_373, pb_x, pc_x, pc_y, pdg0_250, pdg0_251, \
                         pdg0_252, pdf_149, pdg1_250, pdg1_251, pdg1_252, \
                         pff_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = pb_x[k] * pdg0_250[k]
                   - f_7 * pc_x[k] * pdg1_250[k];

        t_371[k] = pb_x[k] * pdg0_251[k]
                   - f_7 * pc_x[k] * pdg1_251[k];

        t_372[k] = pb_x[k] * pdg0_252[k]
                   - f_7 * pc_x[k] * pdg1_252[k];

        t_373[k] = f_0 * pdf_149[k]
                   + f_4 * pc_y[k] * pff_249[k];
    }
}

static auto
compute_prim_pfg_three_center_electron_repulsion_0_piece3(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pb, const size_t pc,
                                                          const size_t sfg0, const size_t sff,
                                                          const size_t sfg1, const size_t ppg0,
                                                          const size_t ppg1, const size_t pdg0,
                                                          const size_t pdf, const size_t pdg1,
                                                          const size_t pfd0, const size_t pfd1,
                                                          const size_t pff, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 0.5 / q;
    const auto f_1 = 1.5 / q;
    const auto f_2 = 1.5 / gamma;
    const auto f_3 = 1.5 * p / (gamma * q);
    const auto f_4 = p / q;
    const auto f_5 = 0.5 / gamma;
    const auto f_6 = 0.5 * p / (gamma * q);
    const auto f_7 = gamma / q;
    const auto f_8 = 1.0 / q;
    const auto f_9 = 2.0 / q;
    const auto f_10 = 1.0 / gamma;
    const auto f_11 = p / (gamma * q);
    const auto f_12 = 0.5 / p;
    const auto f_13 = 0.5 * gamma / (p * q);

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

    const auto *sfg0_90 = buffer.data(sfg0 + 90);
    const auto *sfg0_93 = buffer.data(sfg0 + 93);
    const auto *sfg0_100 = buffer.data(sfg0 + 100);
    const auto *sfg0_101 = buffer.data(sfg0 + 101);
    const auto *sfg0_102 = buffer.data(sfg0 + 102);
    const auto *sfg0_104 = buffer.data(sfg0 + 104);

    const auto *sff_66 = buffer.data(sff + 66);
    const auto *sff_67 = buffer.data(sff + 67);
    const auto *sff_69 = buffer.data(sff + 69);
    const auto *sff_99 = buffer.data(sff + 99);

    const auto *sfg1_90 = buffer.data(sfg1 + 90);
    const auto *sfg1_93 = buffer.data(sfg1 + 93);
    const auto *sfg1_100 = buffer.data(sfg1 + 100);
    const auto *sfg1_101 = buffer.data(sfg1 + 101);
    const auto *sfg1_102 = buffer.data(sfg1 + 102);
    const auto *sfg1_104 = buffer.data(sfg1 + 104);

    const auto *ppg0_134 = buffer.data(ppg0 + 134);

    const auto *ppg1_134 = buffer.data(ppg1 + 134);

    const auto *pdg0_254 = buffer.data(pdg0 + 254);
    const auto *pdg0_255 = buffer.data(pdg0 + 255);
    const auto *pdg0_257 = buffer.data(pdg0 + 257);
    const auto *pdg0_258 = buffer.data(pdg0 + 258);
    const auto *pdg0_260 = buffer.data(pdg0 + 260);
    const auto *pdg0_265 = buffer.data(pdg0 + 265);
    const auto *pdg0_266 = buffer.data(pdg0 + 266);
    const auto *pdg0_267 = buffer.data(pdg0 + 267);
    const auto *pdg0_269 = buffer.data(pdg0 + 269);

    const auto *pdf_150 = buffer.data(pdf + 150);
    const auto *pdf_152 = buffer.data(pdf + 152);
    const auto *pdf_159 = buffer.data(pdf + 159);
    const auto *pdf_160 = buffer.data(pdf + 160);
    const auto *pdf_162 = buffer.data(pdf + 162);
    const auto *pdf_166 = buffer.data(pdf + 166);
    const auto *pdf_167 = buffer.data(pdf + 167);
    const auto *pdf_168 = buffer.data(pdf + 168);
    const auto *pdf_169 = buffer.data(pdf + 169);
    const auto *pdf_170 = buffer.data(pdf + 170);
    const auto *pdf_172 = buffer.data(pdf + 172);
    const auto *pdf_173 = buffer.data(pdf + 173);
    const auto *pdf_175 = buffer.data(pdf + 175);
    const auto *pdf_176 = buffer.data(pdf + 176);
    const auto *pdf_177 = buffer.data(pdf + 177);
    const auto *pdf_178 = buffer.data(pdf + 178);
    const auto *pdf_179 = buffer.data(pdf + 179);

    const auto *pdg1_254 = buffer.data(pdg1 + 254);
    const auto *pdg1_255 = buffer.data(pdg1 + 255);
    const auto *pdg1_257 = buffer.data(pdg1 + 257);
    const auto *pdg1_258 = buffer.data(pdg1 + 258);
    const auto *pdg1_260 = buffer.data(pdg1 + 260);
    const auto *pdg1_265 = buffer.data(pdg1 + 265);
    const auto *pdg1_266 = buffer.data(pdg1 + 266);
    const auto *pdg1_267 = buffer.data(pdg1 + 267);
    const auto *pdg1_269 = buffer.data(pdg1 + 269);

    const auto *pfd0_158 = buffer.data(pfd0 + 158);
    const auto *pfd0_161 = buffer.data(pfd0 + 161);
    const auto *pfd0_162 = buffer.data(pfd0 + 162);
    const auto *pfd0_164 = buffer.data(pfd0 + 164);
    const auto *pfd0_165 = buffer.data(pfd0 + 165);
    const auto *pfd0_166 = buffer.data(pfd0 + 166);
    const auto *pfd0_167 = buffer.data(pfd0 + 167);
    const auto *pfd0_171 = buffer.data(pfd0 + 171);
    const auto *pfd0_174 = buffer.data(pfd0 + 174);
    const auto *pfd0_176 = buffer.data(pfd0 + 176);
    const auto *pfd0_177 = buffer.data(pfd0 + 177);
    const auto *pfd0_178 = buffer.data(pfd0 + 178);
    const auto *pfd0_179 = buffer.data(pfd0 + 179);

    const auto *pfd1_158 = buffer.data(pfd1 + 158);
    const auto *pfd1_161 = buffer.data(pfd1 + 161);
    const auto *pfd1_162 = buffer.data(pfd1 + 162);
    const auto *pfd1_164 = buffer.data(pfd1 + 164);
    const auto *pfd1_165 = buffer.data(pfd1 + 165);
    const auto *pfd1_166 = buffer.data(pfd1 + 166);
    const auto *pfd1_167 = buffer.data(pfd1 + 167);
    const auto *pfd1_171 = buffer.data(pfd1 + 171);
    const auto *pfd1_174 = buffer.data(pfd1 + 174);
    const auto *pfd1_176 = buffer.data(pfd1 + 176);
    const auto *pfd1_177 = buffer.data(pfd1 + 177);
    const auto *pfd1_178 = buffer.data(pfd1 + 178);
    const auto *pfd1_179 = buffer.data(pfd1 + 179);

    const auto *pff_250 = buffer.data(pff + 250);
    const auto *pff_252 = buffer.data(pff + 252);
    const auto *pff_256 = buffer.data(pff + 256);
    const auto *pff_257 = buffer.data(pff + 257);
    const auto *pff_258 = buffer.data(pff + 258);
    const auto *pff_259 = buffer.data(pff + 259);
    const auto *pff_260 = buffer.data(pff + 260);
    const auto *pff_262 = buffer.data(pff + 262);
    const auto *pff_265 = buffer.data(pff + 265);
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
    const auto *pff_280 = buffer.data(pff + 280);
    const auto *pff_282 = buffer.data(pff + 282);
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

#pragma omp simd aligned(t_374, t_375, t_376, t_377, pb_x, pc_x, pc_y, pdg0_254, pdg0_255, \
                         pdg0_257, pdf_170, pdf_172, pdg1_254, pdg1_255, pdg1_257, \
                         pff_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = pb_x[k] * pdg0_254[k]
                   - f_7 * pc_x[k] * pdg1_254[k];

        t_375[k] = pb_x[k] * pdg0_255[k]
                   + f_9 * pdf_170[k]
                   - f_7 * pc_x[k] * pdg1_255[k];

        t_376[k] = f_4 * pc_y[k] * pff_250[k];

        t_377[k] = pb_x[k] * pdg0_257[k]
                   + f_1 * pdf_172[k]
                   - f_7 * pc_x[k] * pdg1_257[k];
    }

#pragma omp simd aligned(t_378, t_379, t_380, t_381, pb_x, pc_x, pc_y, pdg0_258, pdg0_260, \
                         pdf_173, pdf_175, pdf_176, pdg1_258, pdg1_260, pff_252, \
                         pff_256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_378[k] = pb_x[k] * pdg0_258[k]
                   + f_8 * pdf_173[k]
                   - f_7 * pc_x[k] * pdg1_258[k];

        t_379[k] = f_4 * pc_y[k] * pff_252[k];

        t_380[k] = pb_x[k] * pdg0_260[k]
                   + f_8 * pdf_175[k]
                   - f_7 * pc_x[k] * pdg1_260[k];

        t_381[k] = f_0 * pdf_176[k]
                   + f_4 * pc_x[k] * pff_256[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pb_x, pc_x, pdg0_265, pdf_177, pdf_178, \
                         pdf_179, pdg1_265, pff_257, pff_258, pff_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_0 * pdf_177[k]
                   + f_4 * pc_x[k] * pff_257[k];

        t_383[k] = f_0 * pdf_178[k]
                   + f_4 * pc_x[k] * pff_258[k];

        t_384[k] = f_0 * pdf_179[k]
                   + f_4 * pc_x[k] * pff_259[k];

        t_385[k] = pb_x[k] * pdg0_265[k]
                   - f_7 * pc_x[k] * pdg1_265[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, t_389, pb_x, pc_x, pc_y, pdg0_266, pdg0_267, \
                         pdg0_269, pdg1_266, pdg1_267, pdg1_269, \
                         pff_259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = pb_x[k] * pdg0_266[k]
                   - f_7 * pc_x[k] * pdg1_266[k];

        t_387[k] = pb_x[k] * pdg0_267[k]
                   - f_7 * pc_x[k] * pdg1_267[k];

        t_388[k] = f_4 * pc_y[k] * pff_259[k];

        t_389[k] = pb_x[k] * pdg0_269[k]
                   - f_7 * pc_x[k] * pdg1_269[k];
    }

#pragma omp simd aligned(t_390, t_391, t_392, pa_z, pc_x, pc_y, pc_z, sfg0_90, sfg1_90, \
                         pdf_150, pfd0_158, pfd1_158, pff_260, \
                         pff_262 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_390[k] = pa_z[k] * sfg0_90[k]
                   - f_7 * pc_z[k] * sfg1_90[k];

        t_391[k] = f_1 * pdf_150[k]
                   + f_4 * pc_y[k] * pff_260[k];

        t_392[k] = f_10 * pfd0_158[k]
                   - f_11 * pfd1_158[k]
                   + f_4 * pc_x[k] * pff_262[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_z, pc_x, pc_y, pc_z, sfg0_93, sfg1_93, \
                         pdf_152, pfd0_161, pfd1_161, pff_262, pff_265, \
                         pff_266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = pa_z[k] * sfg0_93[k]
                   - f_7 * pc_z[k] * sfg1_93[k];

        t_394[k] = f_1 * pdf_152[k]
                   + f_4 * pc_y[k] * pff_262[k];

        t_395[k] = f_5 * pfd0_161[k]
                   - f_6 * pfd1_161[k]
                   + f_4 * pc_x[k] * pff_265[k];

        t_396[k] = f_4 * pc_x[k] * pff_266[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, t_401, pa_z, pc_x, pc_z, sfg0_100, \
                         sfg0_101, sff_66, sfg1_100, sfg1_101, pff_267, pff_268, \
                         pff_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_4 * pc_x[k] * pff_267[k];

        t_398[k] = f_4 * pc_x[k] * pff_268[k];

        t_399[k] = f_4 * pc_x[k] * pff_269[k];

        t_400[k] = pa_z[k] * sfg0_100[k]
                   - f_7 * pc_z[k] * sfg1_100[k];

        t_401[k] = pa_z[k] * sfg0_101[k]
                   + f_0 * sff_66[k]
                   - f_7 * pc_z[k] * sfg1_101[k];
    }

#pragma omp simd aligned(t_402, t_403, t_404, pa_z, pc_y, pc_z, sfg0_102, sfg0_104, sff_67, \
                         sff_69, sfg1_102, sfg1_104, pdf_159, pff_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_402[k] = pa_z[k] * sfg0_102[k]
                   + f_8 * sff_67[k]
                   - f_7 * pc_z[k] * sfg1_102[k];

        t_403[k] = f_1 * pdf_159[k]
                   + f_4 * pc_y[k] * pff_269[k];

        t_404[k] = pa_z[k] * sfg0_104[k]
                   + f_9 * sff_69[k]
                   - f_7 * pc_z[k] * sfg1_104[k];
    }

#pragma omp simd aligned(t_405, t_406, t_407, t_408, pc_x, pc_y, pdf_160, pfd0_162, pfd0_164, \
                         pfd0_165, pfd1_162, pfd1_164, pfd1_165, pff_270, pff_272, \
                         pff_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_405[k] = f_2 * pfd0_162[k]
                   - f_3 * pfd1_162[k]
                   + f_4 * pc_x[k] * pff_270[k];

        t_406[k] = f_8 * pdf_160[k]
                   + f_4 * pc_y[k] * pff_270[k];

        t_407[k] = f_10 * pfd0_164[k]
                   - f_11 * pfd1_164[k]
                   + f_4 * pc_x[k] * pff_272[k];

        t_408[k] = f_5 * pfd0_165[k]
                   - f_6 * pfd1_165[k]
                   + f_4 * pc_x[k] * pff_273[k];
    }

#pragma omp simd aligned(t_409, t_410, t_411, t_412, t_413, pc_x, pc_y, pdf_162, pfd0_167, \
                         pfd1_167, pff_272, pff_275, pff_276, pff_277, \
                         pff_278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_409[k] = f_8 * pdf_162[k]
                   + f_4 * pc_y[k] * pff_272[k];

        t_410[k] = f_5 * pfd0_167[k]
                   - f_6 * pfd1_167[k]
                   + f_4 * pc_x[k] * pff_275[k];

        t_411[k] = f_4 * pc_x[k] * pff_276[k];

        t_412[k] = f_4 * pc_x[k] * pff_277[k];

        t_413[k] = f_4 * pc_x[k] * pff_278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pc_x, pc_y, pdf_166, pdf_167, pfd0_165, \
                         pfd0_166, pfd1_165, pfd1_166, pff_276, pff_277, \
                         pff_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_4 * pc_x[k] * pff_279[k];

        t_415[k] = f_8 * pdf_166[k]
                   + f_2 * pfd0_165[k]
                   - f_3 * pfd1_165[k]
                   + f_4 * pc_y[k] * pff_276[k];

        t_416[k] = f_8 * pdf_167[k]
                   + f_10 * pfd0_166[k]
                   - f_11 * pfd1_166[k]
                   + f_4 * pc_y[k] * pff_277[k];
    }

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pc_y, ppg0_134, ppg1_134, pdg0_254, \
                         pdf_168, pdf_169, pdg1_254, pfd0_167, pfd1_167, pff_278, \
                         pff_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = f_8 * pdf_168[k]
                   + f_5 * pfd0_167[k]
                   - f_6 * pfd1_167[k]
                   + f_4 * pc_y[k] * pff_278[k];

        t_418[k] = f_8 * pdf_169[k]
                   + f_4 * pc_y[k] * pff_279[k];

        t_419[k] = f_12 * ppg0_134[k]
                   - f_13 * ppg1_134[k]
                   + pb_y[k] * pdg0_254[k]
                   - f_7 * pc_y[k] * pdg1_254[k];
    }

#pragma omp simd aligned(t_420, t_421, t_422, t_423, pb_y, pc_x, pc_y, pdg0_255, pdg0_257, \
                         pdf_170, pdg1_255, pdg1_257, pfd0_171, pfd1_171, pff_280, \
                         pff_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_420[k] = pb_y[k] * pdg0_255[k]
                   - f_7 * pc_y[k] * pdg1_255[k];

        t_421[k] = f_0 * pdf_170[k]
                   + f_4 * pc_y[k] * pff_280[k];

        t_422[k] = pb_y[k] * pdg0_257[k]
                   - f_7 * pc_y[k] * pdg1_257[k];

        t_423[k] = f_5 * pfd0_171[k]
                   - f_6 * pfd1_171[k]
                   + f_4 * pc_x[k] * pff_283[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, t_427, t_428, pb_y, pc_x, pc_y, pdg0_260, \
                         pdf_172, pdg1_260, pff_282, pff_286, pff_287, \
                         pff_288 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_0 * pdf_172[k]
                   + f_4 * pc_y[k] * pff_282[k];

        t_425[k] = pb_y[k] * pdg0_260[k]
                   - f_7 * pc_y[k] * pdg1_260[k];

        t_426[k] = f_4 * pc_x[k] * pff_286[k];

        t_427[k] = f_4 * pc_x[k] * pff_287[k];

        t_428[k] = f_4 * pc_x[k] * pff_288[k];
    }

#pragma omp simd aligned(t_429, t_430, t_431, pb_y, pc_x, pc_y, pdg0_265, pdg0_266, pdf_176, \
                         pdf_177, pdg1_265, pdg1_266, pff_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_429[k] = f_4 * pc_x[k] * pff_289[k];

        t_430[k] = pb_y[k] * pdg0_265[k]
                   + f_9 * pdf_176[k]
                   - f_7 * pc_y[k] * pdg1_265[k];

        t_431[k] = pb_y[k] * pdg0_266[k]
                   + f_1 * pdf_177[k]
                   - f_7 * pc_y[k] * pdg1_266[k];
    }

#pragma omp simd aligned(t_432, t_433, t_434, pb_y, pc_y, pdg0_267, pdg0_269, pdf_178, \
                         pdf_179, pdg1_267, pdg1_269, pff_289 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_432[k] = pb_y[k] * pdg0_267[k]
                   + f_8 * pdf_178[k]
                   - f_7 * pc_y[k] * pdg1_267[k];

        t_433[k] = f_0 * pdf_179[k]
                   + f_4 * pc_y[k] * pff_289[k];

        t_434[k] = pb_y[k] * pdg0_269[k]
                   - f_7 * pc_y[k] * pdg1_269[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, t_438, t_439, pc_x, pc_y, pfd0_174, pfd0_176, \
                         pfd0_177, pfd1_174, pfd1_176, pfd1_177, pff_290, pff_292, \
                         pff_293 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = f_2 * pfd0_174[k]
                   - f_3 * pfd1_174[k]
                   + f_4 * pc_x[k] * pff_290[k];

        t_436[k] = f_4 * pc_y[k] * pff_290[k];

        t_437[k] = f_10 * pfd0_176[k]
                   - f_11 * pfd1_176[k]
                   + f_4 * pc_x[k] * pff_292[k];

        t_438[k] = f_5 * pfd0_177[k]
                   - f_6 * pfd1_177[k]
                   + f_4 * pc_x[k] * pff_293[k];

        t_439[k] = f_4 * pc_y[k] * pff_292[k];
    }

#pragma omp simd aligned(t_440, t_441, t_442, t_443, t_444, pc_x, pfd0_179, pfd1_179, pff_295, \
                         pff_296, pff_297, pff_298, pff_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_440[k] = f_5 * pfd0_179[k]
                   - f_6 * pfd1_179[k]
                   + f_4 * pc_x[k] * pff_295[k];

        t_441[k] = f_4 * pc_x[k] * pff_296[k];

        t_442[k] = f_4 * pc_x[k] * pff_297[k];

        t_443[k] = f_4 * pc_x[k] * pff_298[k];

        t_444[k] = f_4 * pc_x[k] * pff_299[k];
    }

#pragma omp simd aligned(t_445, t_446, t_447, t_448, pc_y, pfd0_177, pfd0_178, pfd0_179, \
                         pfd1_177, pfd1_178, pfd1_179, pff_296, pff_297, pff_298, \
                         pff_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_445[k] = f_2 * pfd0_177[k]
                   - f_3 * pfd1_177[k]
                   + f_4 * pc_y[k] * pff_296[k];

        t_446[k] = f_10 * pfd0_178[k]
                   - f_11 * pfd1_178[k]
                   + f_4 * pc_y[k] * pff_297[k];

        t_447[k] = f_5 * pfd0_179[k]
                   - f_6 * pfd1_179[k]
                   + f_4 * pc_y[k] * pff_298[k];

        t_448[k] = f_4 * pc_y[k] * pff_299[k];
    }

#pragma omp simd aligned(t_449, pc_z, sff_99, pdf_179, pfd0_179, pfd1_179, \
                         pff_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_449[k] = f_0 * sff_99[k]
                   + f_1 * pdf_179[k]
                   + f_2 * pfd0_179[k]
                   - f_3 * pfd1_179[k]
                   + f_4 * pc_z[k] * pff_299[k];
    }
}

auto
compute_prim_pfg_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pb,
                                                   const size_t pc, const size_t sfg0,
                                                   const size_t sff, const size_t sfg1,
                                                   const size_t ppg0, const size_t ppg1,
                                                   const size_t pdg0, const size_t pdf,
                                                   const size_t pdg1, const size_t pfd0,
                                                   const size_t pfd1, const size_t pff,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_pfg_three_center_electron_repulsion_0_piece0(buffer, target, pa, pb, pc, sfg0,
                                                              sff, sfg1, pdg0, pdf, pdg1, pfd0,
                                                              pfd1, pff, ncols, gamma, p, q);

    compute_prim_pfg_three_center_electron_repulsion_0_piece1(buffer, target, pa, pb, pc, sfg0,
                                                              sff, sfg1, ppg0, ppg1, pdg0, pdf,
                                                              pdg1, pfd0, pfd1, pff, ncols,
                                                              gamma, p, q);

    compute_prim_pfg_three_center_electron_repulsion_0_piece2(buffer, target, pa, pb, pc, sfg0,
                                                              sff, sfg1, ppg0, ppg1, pdg0, pdf,
                                                              pdg1, pfd0, pfd1, pff, ncols,
                                                              gamma, p, q);

    compute_prim_pfg_three_center_electron_repulsion_0_piece3(buffer, target, pa, pb, pc, sfg0,
                                                              sff, sfg1, ppg0, ppg1, pdg0, pdf,
                                                              pdg1, pfd0, pfd1, pff, ncols,
                                                              gamma, p, q);
}

}  // namespace simdt3ceri
