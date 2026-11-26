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


#include "SimdElectronRepulsionVrrRecDI.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_di_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t ph,
                                            const size_t pi, const size_t dg0, const size_t dg1,
                                            const size_t dh, const size_t ncols,
                                            const double alpha, const double beta,
                                            const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);
    const auto f_9 = 0.5 / p;
    const auto f_10 = 2.0 / p;
    const auto f_11 = 1.5 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_0 = buffer.data(ph + 0);
    const auto *ph_3 = buffer.data(ph + 3);
    const auto *ph_5 = buffer.data(ph + 5);
    const auto *ph_6 = buffer.data(ph + 6);
    const auto *ph_9 = buffer.data(ph + 9);
    const auto *ph_15 = buffer.data(ph + 15);
    const auto *ph_17 = buffer.data(ph + 17);
    const auto *ph_18 = buffer.data(ph + 18);
    const auto *ph_20 = buffer.data(ph + 20);
    const auto *ph_21 = buffer.data(ph + 21);
    const auto *ph_24 = buffer.data(ph + 24);
    const auto *ph_26 = buffer.data(ph + 26);
    const auto *ph_27 = buffer.data(ph + 27);
    const auto *ph_30 = buffer.data(ph + 30);
    const auto *ph_31 = buffer.data(ph + 31);
    const auto *ph_33 = buffer.data(ph + 33);
    const auto *ph_36 = buffer.data(ph + 36);
    const auto *ph_38 = buffer.data(ph + 38);
    const auto *ph_39 = buffer.data(ph + 39);
    const auto *ph_40 = buffer.data(ph + 40);
    const auto *ph_41 = buffer.data(ph + 41);
    const auto *ph_42 = buffer.data(ph + 42);
    const auto *ph_44 = buffer.data(ph + 44);
    const auto *ph_45 = buffer.data(ph + 45);
    const auto *ph_47 = buffer.data(ph + 47);
    const auto *ph_48 = buffer.data(ph + 48);
    const auto *ph_50 = buffer.data(ph + 50);
    const auto *ph_51 = buffer.data(ph + 51);
    const auto *ph_54 = buffer.data(ph + 54);
    const auto *ph_56 = buffer.data(ph + 56);
    const auto *ph_58 = buffer.data(ph + 58);
    const auto *ph_59 = buffer.data(ph + 59);
    const auto *ph_60 = buffer.data(ph + 60);
    const auto *ph_61 = buffer.data(ph + 61);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *pi_0 = buffer.data(pi + 0);
    const auto *pi_3 = buffer.data(pi + 3);
    const auto *pi_5 = buffer.data(pi + 5);
    const auto *pi_6 = buffer.data(pi + 6);
    const auto *pi_9 = buffer.data(pi + 9);
    const auto *pi_10 = buffer.data(pi + 10);
    const auto *pi_14 = buffer.data(pi + 14);
    const auto *pi_15 = buffer.data(pi + 15);
    const auto *pi_20 = buffer.data(pi + 20);
    const auto *pi_29 = buffer.data(pi + 29);
    const auto *pi_31 = buffer.data(pi + 31);
    const auto *pi_34 = buffer.data(pi + 34);
    const auto *pi_38 = buffer.data(pi + 38);
    const auto *pi_40 = buffer.data(pi + 40);
    const auto *pi_49 = buffer.data(pi + 49);
    const auto *pi_51 = buffer.data(pi + 51);
    const auto *pi_52 = buffer.data(pi + 52);
    const auto *pi_53 = buffer.data(pi + 53);
    const auto *pi_54 = buffer.data(pi + 54);
    const auto *pi_55 = buffer.data(pi + 55);
    const auto *pi_56 = buffer.data(pi + 56);
    const auto *pi_58 = buffer.data(pi + 58);
    const auto *pi_61 = buffer.data(pi + 61);
    const auto *pi_65 = buffer.data(pi + 65);
    const auto *pi_68 = buffer.data(pi + 68);
    const auto *pi_70 = buffer.data(pi + 70);
    const auto *pi_77 = buffer.data(pi + 77);
    const auto *pi_78 = buffer.data(pi + 78);
    const auto *pi_79 = buffer.data(pi + 79);
    const auto *pi_80 = buffer.data(pi + 80);
    const auto *pi_81 = buffer.data(pi + 81);
    const auto *pi_83 = buffer.data(pi + 83);

    const auto *dg0_0 = buffer.data(dg0 + 0);
    const auto *dg0_1 = buffer.data(dg0 + 1);
    const auto *dg0_2 = buffer.data(dg0 + 2);
    const auto *dg0_3 = buffer.data(dg0 + 3);
    const auto *dg0_5 = buffer.data(dg0 + 5);
    const auto *dg0_10 = buffer.data(dg0 + 10);
    const auto *dg0_12 = buffer.data(dg0 + 12);
    const auto *dg0_13 = buffer.data(dg0 + 13);
    const auto *dg0_14 = buffer.data(dg0 + 14);
    const auto *dg0_45 = buffer.data(dg0 + 45);
    const auto *dg0_48 = buffer.data(dg0 + 48);
    const auto *dg0_50 = buffer.data(dg0 + 50);
    const auto *dg0_51 = buffer.data(dg0 + 51);
    const auto *dg0_54 = buffer.data(dg0 + 54);
    const auto *dg0_55 = buffer.data(dg0 + 55);
    const auto *dg0_56 = buffer.data(dg0 + 56);
    const auto *dg0_57 = buffer.data(dg0 + 57);
    const auto *dg0_59 = buffer.data(dg0 + 59);
    const auto *dg0_75 = buffer.data(dg0 + 75);
    const auto *dg0_78 = buffer.data(dg0 + 78);
    const auto *dg0_80 = buffer.data(dg0 + 80);
    const auto *dg0_81 = buffer.data(dg0 + 81);
    const auto *dg0_84 = buffer.data(dg0 + 84);
    const auto *dg0_85 = buffer.data(dg0 + 85);

    const auto *dg1_0 = buffer.data(dg1 + 0);
    const auto *dg1_1 = buffer.data(dg1 + 1);
    const auto *dg1_2 = buffer.data(dg1 + 2);
    const auto *dg1_3 = buffer.data(dg1 + 3);
    const auto *dg1_5 = buffer.data(dg1 + 5);
    const auto *dg1_10 = buffer.data(dg1 + 10);
    const auto *dg1_12 = buffer.data(dg1 + 12);
    const auto *dg1_13 = buffer.data(dg1 + 13);
    const auto *dg1_14 = buffer.data(dg1 + 14);
    const auto *dg1_45 = buffer.data(dg1 + 45);
    const auto *dg1_48 = buffer.data(dg1 + 48);
    const auto *dg1_50 = buffer.data(dg1 + 50);
    const auto *dg1_51 = buffer.data(dg1 + 51);
    const auto *dg1_54 = buffer.data(dg1 + 54);
    const auto *dg1_55 = buffer.data(dg1 + 55);
    const auto *dg1_56 = buffer.data(dg1 + 56);
    const auto *dg1_57 = buffer.data(dg1 + 57);
    const auto *dg1_59 = buffer.data(dg1 + 59);
    const auto *dg1_75 = buffer.data(dg1 + 75);
    const auto *dg1_78 = buffer.data(dg1 + 78);
    const auto *dg1_80 = buffer.data(dg1 + 80);
    const auto *dg1_81 = buffer.data(dg1 + 81);
    const auto *dg1_84 = buffer.data(dg1 + 84);
    const auto *dg1_85 = buffer.data(dg1 + 85);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_1 = buffer.data(dh + 1);
    const auto *dh_2 = buffer.data(dh + 2);
    const auto *dh_3 = buffer.data(dh + 3);
    const auto *dh_5 = buffer.data(dh + 5);
    const auto *dh_6 = buffer.data(dh + 6);
    const auto *dh_8 = buffer.data(dh + 8);
    const auto *dh_9 = buffer.data(dh + 9);
    const auto *dh_10 = buffer.data(dh + 10);
    const auto *dh_14 = buffer.data(dh + 14);
    const auto *dh_15 = buffer.data(dh + 15);
    const auto *dh_17 = buffer.data(dh + 17);
    const auto *dh_18 = buffer.data(dh + 18);
    const auto *dh_19 = buffer.data(dh + 19);
    const auto *dh_20 = buffer.data(dh + 20);
    const auto *dh_21 = buffer.data(dh + 21);
    const auto *dh_22 = buffer.data(dh + 22);
    const auto *dh_24 = buffer.data(dh + 24);
    const auto *dh_26 = buffer.data(dh + 26);
    const auto *dh_27 = buffer.data(dh + 27);
    const auto *dh_30 = buffer.data(dh + 30);
    const auto *dh_31 = buffer.data(dh + 31);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_38 = buffer.data(dh + 38);
    const auto *dh_39 = buffer.data(dh + 39);
    const auto *dh_40 = buffer.data(dh + 40);
    const auto *dh_42 = buffer.data(dh + 42);
    const auto *dh_44 = buffer.data(dh + 44);
    const auto *dh_45 = buffer.data(dh + 45);
    const auto *dh_47 = buffer.data(dh + 47);
    const auto *dh_48 = buffer.data(dh + 48);
    const auto *dh_51 = buffer.data(dh + 51);
    const auto *dh_56 = buffer.data(dh + 56);
    const auto *dh_58 = buffer.data(dh + 58);
    const auto *dh_59 = buffer.data(dh + 59);
    const auto *dh_60 = buffer.data(dh + 60);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_63 = buffer.data(dh + 63);
    const auto *dh_64 = buffer.data(dh + 64);
    const auto *dh_66 = buffer.data(dh + 66);
    const auto *dh_68 = buffer.data(dh + 68);
    const auto *dh_69 = buffer.data(dh + 69);
    const auto *dh_72 = buffer.data(dh + 72);
    const auto *dh_73 = buffer.data(dh + 73);
    const auto *dh_75 = buffer.data(dh + 75);
    const auto *dh_77 = buffer.data(dh + 77);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_79 = buffer.data(dh + 79);
    const auto *dh_80 = buffer.data(dh + 80);
    const auto *dh_81 = buffer.data(dh + 81);
    const auto *dh_82 = buffer.data(dh + 82);
    const auto *dh_83 = buffer.data(dh + 83);
    const auto *dh_86 = buffer.data(dh + 86);
    const auto *dh_87 = buffer.data(dh + 87);
    const auto *dh_89 = buffer.data(dh + 89);
    const auto *dh_90 = buffer.data(dh + 90);
    const auto *dh_93 = buffer.data(dh + 93);
    const auto *dh_99 = buffer.data(dh + 99);
    const auto *dh_100 = buffer.data(dh + 100);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);
    const auto *dh_103 = buffer.data(dh + 103);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_105 = buffer.data(dh + 105);
    const auto *dh_107 = buffer.data(dh + 107);
    const auto *dh_108 = buffer.data(dh + 108);
    const auto *dh_110 = buffer.data(dh + 110);
    const auto *dh_111 = buffer.data(dh + 111);
    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_115 = buffer.data(dh + 115);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, ph_0, dg0_0, dg1_0, \
                         dh_0, dh_1, dh_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ph_0[k]
                 + f_1 * dg0_0[k]
                 - f_2 * dg1_0[k]
                 + pb_x[k] * dh_0[k];

        t_1[k] = pb_y[k] * dh_0[k];

        t_2[k] = pb_z[k] * dh_0[k];

        t_3[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pb_y[k] * dh_1[k];

        t_4[k] = pb_y[k] * dh_2[k];

        t_5[k] = f_3 * dg0_0[k]
                 - f_4 * dg1_0[k]
                 + pb_z[k] * dh_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, pb_y, pb_z, dg0_1, dg0_2, dg0_3, dg1_1, \
                         dg1_2, dg1_3, dh_3, dh_5, dh_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = f_5 * dg0_1[k]
                 - f_6 * dg1_1[k]
                 + pb_y[k] * dh_3[k];

        t_7[k] = pb_z[k] * dh_3[k];

        t_8[k] = pb_y[k] * dh_5[k];

        t_9[k] = f_5 * dg0_2[k]
                 - f_6 * dg1_2[k]
                 + pb_z[k] * dh_5[k];

        t_10[k] = f_7 * dg0_3[k]
                  - f_8 * dg1_3[k]
                  + pb_y[k] * dh_6[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, t_14, t_15, pb_x, pb_y, pb_z, ph_15, dg0_5, dg1_5, \
                         dh_6, dh_8, dh_9, dh_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pb_z[k] * dh_6[k];

        t_12[k] = f_3 * dg0_5[k]
                  - f_4 * dg1_5[k]
                  + pb_y[k] * dh_8[k];

        t_13[k] = pb_y[k] * dh_9[k];

        t_14[k] = f_7 * dg0_5[k]
                  - f_8 * dg1_5[k]
                  + pb_z[k] * dh_9[k];

        t_15[k] = f_0 * ph_15[k]
                  + pb_x[k] * dh_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, t_20, pb_x, pb_y, pb_z, ph_17, ph_18, ph_20, \
                         dh_10, dh_14, dh_17, dh_18, dh_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = pb_z[k] * dh_10[k];

        t_17[k] = f_0 * ph_17[k]
                  + pb_x[k] * dh_17[k];

        t_18[k] = f_0 * ph_18[k]
                  + pb_x[k] * dh_18[k];

        t_19[k] = pb_y[k] * dh_14[k];

        t_20[k] = f_0 * ph_20[k]
                  + pb_x[k] * dh_20[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pb_y, pb_z, dg0_10, dg0_12, dg0_13, dg1_10, \
                         dg1_12, dg1_13, dh_15, dh_17, dh_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_1 * dg0_10[k]
                  - f_2 * dg1_10[k]
                  + pb_y[k] * dh_15[k];

        t_22[k] = pb_z[k] * dh_15[k];

        t_23[k] = f_7 * dg0_12[k]
                  - f_8 * dg1_12[k]
                  + pb_y[k] * dh_17[k];

        t_24[k] = f_5 * dg0_13[k]
                  - f_6 * dg1_13[k]
                  + pb_y[k] * dh_18[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, t_30, pa_y, pb_y, pb_z, ph_0, pi_0, \
                         dg0_14, dg1_14, dh_19, dh_20, dh_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_3 * dg0_14[k]
                  - f_4 * dg1_14[k]
                  + pb_y[k] * dh_19[k];

        t_26[k] = pb_y[k] * dh_20[k];

        t_27[k] = f_1 * dg0_14[k]
                  - f_2 * dg1_14[k]
                  + pb_z[k] * dh_20[k];

        t_28[k] = pa_y[k] * pi_0[k];

        t_29[k] = f_9 * ph_0[k]
                  + pb_y[k] * dh_21[k];

        t_30[k] = pb_z[k] * dh_21[k];
    }

#pragma omp simd aligned(t_31, t_32, t_33, t_34, t_35, pa_x, pa_y, pb_z, ph_24, ph_27, pi_5, \
                         pi_31, pi_34, dh_22, dh_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_10 * ph_24[k]
                  + pa_x[k] * pi_31[k];

        t_32[k] = pb_z[k] * dh_22[k];

        t_33[k] = pa_y[k] * pi_5[k];

        t_34[k] = f_11 * ph_27[k]
                  + pa_x[k] * pi_34[k];

        t_35[k] = pb_z[k] * dh_24[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, t_39, pa_x, pa_y, pb_y, pb_z, ph_5, ph_31, pi_9, \
                         pi_38, dh_26, dh_27 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = f_9 * ph_5[k]
                  + pb_y[k] * dh_26[k];

        t_37[k] = pa_y[k] * pi_9[k];

        t_38[k] = f_0 * ph_31[k]
                  + pa_x[k] * pi_38[k];

        t_39[k] = pb_z[k] * dh_27[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, t_43, pa_x, pa_y, pb_x, pb_y, ph_9, ph_33, ph_36, \
                         pi_14, pi_40, dh_30, dh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = f_0 * ph_33[k]
                  + pa_x[k] * pi_40[k];

        t_41[k] = f_9 * ph_9[k]
                  + pb_y[k] * dh_30[k];

        t_42[k] = pa_y[k] * pi_14[k];

        t_43[k] = f_9 * ph_36[k]
                  + pb_x[k] * dh_36[k];
    }

#pragma omp simd aligned(t_44, t_45, t_46, t_47, t_48, pa_y, pb_x, pb_z, ph_38, ph_39, ph_40, \
                         pi_20, dh_31, dh_38, dh_39, dh_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_44[k] = pb_z[k] * dh_31[k];

        t_45[k] = f_9 * ph_38[k]
                  + pb_x[k] * dh_38[k];

        t_46[k] = f_9 * ph_39[k]
                  + pb_x[k] * dh_39[k];

        t_47[k] = f_9 * ph_40[k]
                  + pb_x[k] * dh_40[k];

        t_48[k] = pa_y[k] * pi_20[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, t_52, t_53, t_54, t_55, pa_x, pb_z, pi_49, pi_51, \
                         pi_52, pi_53, pi_54, pi_55, dh_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_x[k] * pi_49[k];

        t_50[k] = pb_z[k] * dh_36[k];

        t_51[k] = pa_x[k] * pi_51[k];

        t_52[k] = pa_x[k] * pi_52[k];

        t_53[k] = pa_x[k] * pi_53[k];

        t_54[k] = pa_x[k] * pi_54[k];

        t_55[k] = pa_x[k] * pi_55[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, t_59, t_60, pa_z, pb_y, pb_z, ph_0, pi_0, pi_3, \
                         dh_42, dh_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = pa_z[k] * pi_0[k];

        t_57[k] = pb_y[k] * dh_42[k];

        t_58[k] = f_9 * ph_0[k]
                  + pb_z[k] * dh_42[k];

        t_59[k] = pa_z[k] * pi_3[k];

        t_60[k] = pb_y[k] * dh_44[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pa_z, pb_y, pb_z, ph_3, ph_47, pi_6, \
                         pi_61, dh_45, dh_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_10 * ph_47[k]
                  + pa_x[k] * pi_61[k];

        t_62[k] = pa_z[k] * pi_6[k];

        t_63[k] = f_9 * ph_3[k]
                  + pb_z[k] * dh_45[k];

        t_64[k] = pb_y[k] * dh_47[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pa_z, pb_z, ph_6, ph_51, ph_54, pi_10, \
                         pi_65, pi_68, dh_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * ph_51[k]
                  + pa_x[k] * pi_65[k];

        t_66[k] = pa_z[k] * pi_10[k];

        t_67[k] = f_9 * ph_6[k]
                  + pb_z[k] * dh_48[k];

        t_68[k] = f_0 * ph_54[k]
                  + pa_x[k] * pi_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pa_z, pb_x, pb_y, ph_56, ph_58, pi_15, \
                         pi_70, dh_51, dh_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pb_y[k] * dh_51[k];

        t_70[k] = f_0 * ph_56[k]
                  + pa_x[k] * pi_70[k];

        t_71[k] = pa_z[k] * pi_15[k];

        t_72[k] = f_9 * ph_58[k]
                  + pb_x[k] * dh_58[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, t_77, pa_x, pb_x, pb_y, ph_59, ph_60, ph_62, \
                         pi_77, dh_56, dh_59, dh_60, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_9 * ph_59[k]
                  + pb_x[k] * dh_59[k];

        t_74[k] = f_9 * ph_60[k]
                  + pb_x[k] * dh_60[k];

        t_75[k] = pb_y[k] * dh_56[k];

        t_76[k] = f_9 * ph_62[k]
                  + pb_x[k] * dh_62[k];

        t_77[k] = pa_x[k] * pi_77[k];
    }

#pragma omp simd aligned(t_78, t_79, t_80, t_81, t_82, t_83, pa_x, pb_y, pi_78, pi_79, pi_80, \
                         pi_81, pi_83, dh_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_78[k] = pa_x[k] * pi_78[k];

        t_79[k] = pa_x[k] * pi_79[k];

        t_80[k] = pa_x[k] * pi_80[k];

        t_81[k] = pa_x[k] * pi_81[k];

        t_82[k] = pb_y[k] * dh_62[k];

        t_83[k] = pa_x[k] * pi_83[k];
    }

#pragma omp simd aligned(t_84, t_85, t_86, t_87, t_88, pb_x, pb_y, pb_z, ph_21, dg0_45, \
                         dg0_48, dg1_45, dg1_48, dh_63, dh_64, dh_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_1 * dg0_45[k]
                  - f_2 * dg1_45[k]
                  + pb_x[k] * dh_63[k];

        t_85[k] = f_0 * ph_21[k]
                  + pb_y[k] * dh_63[k];

        t_86[k] = pb_z[k] * dh_63[k];

        t_87[k] = f_7 * dg0_48[k]
                  - f_8 * dg1_48[k]
                  + pb_x[k] * dh_66[k];

        t_88[k] = pb_z[k] * dh_64[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, t_92, pb_x, pb_y, pb_z, ph_26, dg0_50, dg0_51, \
                         dg1_50, dg1_51, dh_66, dh_68, dh_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_7 * dg0_50[k]
                  - f_8 * dg1_50[k]
                  + pb_x[k] * dh_68[k];

        t_90[k] = f_5 * dg0_51[k]
                  - f_6 * dg1_51[k]
                  + pb_x[k] * dh_69[k];

        t_91[k] = pb_z[k] * dh_66[k];

        t_92[k] = f_0 * ph_26[k]
                  + pb_y[k] * dh_68[k];
    }

#pragma omp simd aligned(t_93, t_94, t_95, t_96, pb_x, pb_z, dg0_54, dg0_55, dg0_57, dg1_54, \
                         dg1_55, dg1_57, dh_69, dh_72, dh_73, dh_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_93[k] = f_5 * dg0_54[k]
                  - f_6 * dg1_54[k]
                  + pb_x[k] * dh_72[k];

        t_94[k] = f_3 * dg0_55[k]
                  - f_4 * dg1_55[k]
                  + pb_x[k] * dh_73[k];

        t_95[k] = pb_z[k] * dh_69[k];

        t_96[k] = f_3 * dg0_57[k]
                  - f_4 * dg1_57[k]
                  + pb_x[k] * dh_75[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, t_100, t_101, pb_x, pb_y, ph_30, dg0_59, dg1_59, \
                         dh_72, dh_77, dh_78, dh_79, dh_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_0 * ph_30[k]
                  + pb_y[k] * dh_72[k];

        t_98[k] = f_3 * dg0_59[k]
                  - f_4 * dg1_59[k]
                  + pb_x[k] * dh_77[k];

        t_99[k] = pb_x[k] * dh_78[k];

        t_100[k] = pb_x[k] * dh_79[k];

        t_101[k] = pb_x[k] * dh_80[k];
    }

#pragma omp simd aligned(t_102, t_103, t_104, t_105, t_106, pb_x, pb_y, pb_z, ph_36, dg0_55, \
                         dg1_55, dh_78, dh_81, dh_82, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_102[k] = pb_x[k] * dh_81[k];

        t_103[k] = pb_x[k] * dh_82[k];

        t_104[k] = pb_x[k] * dh_83[k];

        t_105[k] = f_0 * ph_36[k]
                   + f_1 * dg0_55[k]
                   - f_2 * dg1_55[k]
                   + pb_y[k] * dh_78[k];

        t_106[k] = pb_z[k] * dh_78[k];
    }

#pragma omp simd aligned(t_107, t_108, t_109, pb_z, dg0_55, dg0_56, dg0_57, dg1_55, dg1_56, \
                         dg1_57, dh_79, dh_80, dh_81 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_3 * dg0_55[k]
                   - f_4 * dg1_55[k]
                   + pb_z[k] * dh_79[k];

        t_108[k] = f_5 * dg0_56[k]
                   - f_6 * dg1_56[k]
                   + pb_z[k] * dh_80[k];

        t_109[k] = f_7 * dg0_57[k]
                   - f_8 * dg1_57[k]
                   + pb_z[k] * dh_81[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, t_113, t_114, pa_y, pa_z, pb_y, pb_z, ph_41, \
                         pi_29, pi_56, pi_58, dg0_59, dg1_59, dh_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_0 * ph_41[k]
                   + pb_y[k] * dh_83[k];

        t_111[k] = f_1 * dg0_59[k]
                   - f_2 * dg1_59[k]
                   + pb_z[k] * dh_83[k];

        t_112[k] = pa_y[k] * pi_56[k];

        t_113[k] = pa_z[k] * pi_29[k];

        t_114[k] = pa_y[k] * pi_58[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pa_z, pb_y, pb_z, ph_24, \
                         ph_44, pi_31, pi_34, pi_61, dh_86, dh_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * pi_31[k];

        t_116[k] = f_9 * ph_44[k]
                   + pb_y[k] * dh_86[k];

        t_117[k] = pa_y[k] * pi_61[k];

        t_118[k] = pa_z[k] * pi_34[k];

        t_119[k] = f_9 * ph_24[k]
                   + pb_z[k] * dh_87[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_y, pa_z, pb_y, pb_z, ph_27, ph_47, \
                         pi_38, pi_65, dh_89, dh_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_9 * ph_47[k]
                   + pb_y[k] * dh_89[k];

        t_121[k] = pa_y[k] * pi_65[k];

        t_122[k] = pa_z[k] * pi_38[k];

        t_123[k] = f_9 * ph_27[k]
                   + pb_z[k] * dh_90[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, t_127, t_128, pa_y, pb_x, pb_y, ph_50, ph_51, \
                         pi_68, pi_70, dh_93, dh_99, dh_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_0 * ph_50[k]
                   + pa_y[k] * pi_68[k];

        t_125[k] = f_9 * ph_51[k]
                   + pb_y[k] * dh_93[k];

        t_126[k] = pa_y[k] * pi_70[k];

        t_127[k] = pb_x[k] * dh_99[k];

        t_128[k] = pb_x[k] * dh_100[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, t_132, t_133, t_134, pa_z, pb_x, pb_z, ph_36, \
                         pi_49, dh_99, dh_101, dh_102, dh_103, dh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = pb_x[k] * dh_101[k];

        t_130[k] = pb_x[k] * dh_102[k];

        t_131[k] = pb_x[k] * dh_103[k];

        t_132[k] = pb_x[k] * dh_104[k];

        t_133[k] = pa_z[k] * pi_49[k];

        t_134[k] = f_9 * ph_36[k]
                   + pb_z[k] * dh_99[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, t_138, t_139, pa_y, pb_y, ph_59, ph_60, ph_61, \
                         ph_62, pi_79, pi_80, pi_81, pi_83, dh_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_10 * ph_59[k]
                   + pa_y[k] * pi_79[k];

        t_136[k] = f_11 * ph_60[k]
                   + pa_y[k] * pi_80[k];

        t_137[k] = f_0 * ph_61[k]
                   + pa_y[k] * pi_81[k];

        t_138[k] = f_9 * ph_62[k]
                   + pb_y[k] * dh_104[k];

        t_139[k] = pa_y[k] * pi_83[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, t_144, pb_x, pb_y, pb_z, ph_42, dg0_75, \
                         dg0_78, dg1_75, dg1_78, dh_105, dh_107, \
                         dh_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_1 * dg0_75[k]
                   - f_2 * dg1_75[k]
                   + pb_x[k] * dh_105[k];

        t_141[k] = pb_y[k] * dh_105[k];

        t_142[k] = f_0 * ph_42[k]
                   + pb_z[k] * dh_105[k];

        t_143[k] = f_7 * dg0_78[k]
                   - f_8 * dg1_78[k]
                   + pb_x[k] * dh_108[k];

        t_144[k] = pb_y[k] * dh_107[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, t_148, pb_x, pb_y, pb_z, ph_45, dg0_80, dg0_81, \
                         dg1_80, dg1_81, dh_108, dh_110, dh_111 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_7 * dg0_80[k]
                   - f_8 * dg1_80[k]
                   + pb_x[k] * dh_110[k];

        t_146[k] = f_5 * dg0_81[k]
                   - f_6 * dg1_81[k]
                   + pb_x[k] * dh_111[k];

        t_147[k] = f_0 * ph_45[k]
                   + pb_z[k] * dh_108[k];

        t_148[k] = pb_y[k] * dh_110[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pb_x, pb_z, ph_48, dg0_84, dg0_85, dg1_84, \
                         dg1_85, dh_111, dh_114, dh_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = f_5 * dg0_84[k]
                   - f_6 * dg1_84[k]
                   + pb_x[k] * dh_114[k];

        t_150[k] = f_3 * dg0_85[k]
                   - f_4 * dg1_85[k]
                   + pb_x[k] * dh_115[k];

        t_151[k] = f_0 * ph_48[k]
                   + pb_z[k] * dh_111[k];
    }
}

static auto
compute_prim_di_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pb, const size_t ph, const size_t dg0,
                                            const size_t dg1, const size_t dh,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 1.0 / p;
    const auto f_1 = 2.5 / beta;
    const auto f_2 = 2.5 * alpha / (beta * p);
    const auto f_3 = 0.5 / beta;
    const auto f_4 = 0.5 * alpha / (beta * p);
    const auto f_5 = 1.0 / beta;
    const auto f_6 = alpha / (beta * p);
    const auto f_7 = 1.5 / beta;
    const auto f_8 = 1.5 * alpha / (beta * p);

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ph_57 = buffer.data(ph + 57);
    const auto *ph_62 = buffer.data(ph + 62);

    const auto *dg0_85 = buffer.data(dg0 + 85);
    const auto *dg0_87 = buffer.data(dg0 + 87);
    const auto *dg0_88 = buffer.data(dg0 + 88);
    const auto *dg0_89 = buffer.data(dg0 + 89);

    const auto *dg1_85 = buffer.data(dg1 + 85);
    const auto *dg1_87 = buffer.data(dg1 + 87);
    const auto *dg1_88 = buffer.data(dg1 + 88);
    const auto *dg1_89 = buffer.data(dg1 + 89);

    const auto *dh_114 = buffer.data(dh + 114);
    const auto *dh_117 = buffer.data(dh + 117);
    const auto *dh_119 = buffer.data(dh + 119);
    const auto *dh_120 = buffer.data(dh + 120);
    const auto *dh_121 = buffer.data(dh + 121);
    const auto *dh_122 = buffer.data(dh + 122);
    const auto *dh_123 = buffer.data(dh + 123);
    const auto *dh_124 = buffer.data(dh + 124);
    const auto *dh_125 = buffer.data(dh + 125);

#pragma omp simd aligned(t_152, t_153, t_154, t_155, t_156, pb_x, pb_y, dg0_87, dg0_89, \
                         dg1_87, dg1_89, dh_114, dh_117, dh_119, dh_120, \
                         dh_121 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_3 * dg0_87[k]
                   - f_4 * dg1_87[k]
                   + pb_x[k] * dh_117[k];

        t_153[k] = pb_y[k] * dh_114[k];

        t_154[k] = f_3 * dg0_89[k]
                   - f_4 * dg1_89[k]
                   + pb_x[k] * dh_119[k];

        t_155[k] = pb_x[k] * dh_120[k];

        t_156[k] = pb_x[k] * dh_121[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pb_x, pb_y, dg0_85, dg1_85, \
                         dh_120, dh_122, dh_123, dh_124, dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pb_x[k] * dh_122[k];

        t_158[k] = pb_x[k] * dh_123[k];

        t_159[k] = pb_x[k] * dh_124[k];

        t_160[k] = pb_x[k] * dh_125[k];

        t_161[k] = f_1 * dg0_85[k]
                   - f_2 * dg1_85[k]
                   + pb_y[k] * dh_120[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, pb_y, pb_z, ph_57, dg0_87, dg0_88, dg1_87, \
                         dg1_88, dh_120, dh_122, dh_123 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_0 * ph_57[k]
                   + pb_z[k] * dh_120[k];

        t_163[k] = f_7 * dg0_87[k]
                   - f_8 * dg1_87[k]
                   + pb_y[k] * dh_122[k];

        t_164[k] = f_5 * dg0_88[k]
                   - f_6 * dg1_88[k]
                   + pb_y[k] * dh_123[k];
    }

#pragma omp simd aligned(t_165, t_166, t_167, pb_y, pb_z, ph_62, dg0_89, dg1_89, dh_124, \
                         dh_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_165[k] = f_3 * dg0_89[k]
                   - f_4 * dg1_89[k]
                   + pb_y[k] * dh_124[k];

        t_166[k] = pb_y[k] * dh_125[k];

        t_167[k] = f_0 * ph_62[k]
                   + f_1 * dg0_89[k]
                   - f_2 * dg1_89[k]
                   + pb_z[k] * dh_125[k];
    }
}

auto
compute_prim_di_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t ph, const size_t pi,
                                     const size_t dg0, const size_t dg1, const size_t dh,
                                     const size_t ncols, const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_di_electron_repulsion_0_piece0(buffer, target, pa, pb, ph, pi, dg0, dg1, dh,
                                                ncols, alpha, beta, p);

    compute_prim_di_electron_repulsion_0_piece1(buffer, target, pb, ph, dg0, dg1, dh, ncols,
                                                alpha, beta, p);
}

}  // namespace simdt2ceri
