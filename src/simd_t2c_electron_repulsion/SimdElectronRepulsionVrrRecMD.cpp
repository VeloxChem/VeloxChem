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


#include "SimdElectronRepulsionVrrRecMD.hpp"

#include "SimdAlign.hpp"

namespace simdt2ceri {  // simdt2ceri namespace

static auto
compute_prim_md_electron_repulsion_0_piece0(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kd0,
                                            const size_t kd1, const size_t lp, const size_t ld,
                                            const size_t ms0, const size_t ms1, const size_t mp,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 4.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.5 / p;
    const auto f_8 = 3.0 / alpha;
    const auto f_9 = 3.0 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 3.0 / p;
    const auto f_14 = 2.5 / alpha;
    const auto f_15 = 2.5 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 2.5 / p;
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.5 / p;
    const auto f_22 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kd0_0 = buffer.data(kd0 + 0);
    const auto *kd0_6 = buffer.data(kd0 + 6);
    const auto *kd0_12 = buffer.data(kd0 + 12);
    const auto *kd0_18 = buffer.data(kd0 + 18);
    const auto *kd0_21 = buffer.data(kd0 + 21);
    const auto *kd0_30 = buffer.data(kd0 + 30);
    const auto *kd0_35 = buffer.data(kd0 + 35);
    const auto *kd0_36 = buffer.data(kd0 + 36);
    const auto *kd0_39 = buffer.data(kd0 + 39);
    const auto *kd0_48 = buffer.data(kd0 + 48);
    const auto *kd0_54 = buffer.data(kd0 + 54);
    const auto *kd0_59 = buffer.data(kd0 + 59);
    const auto *kd0_60 = buffer.data(kd0 + 60);
    const auto *kd0_63 = buffer.data(kd0 + 63);
    const auto *kd0_72 = buffer.data(kd0 + 72);
    const auto *kd0_75 = buffer.data(kd0 + 75);
    const auto *kd0_77 = buffer.data(kd0 + 77);
    const auto *kd0_89 = buffer.data(kd0 + 89);
    const auto *kd0_93 = buffer.data(kd0 + 93);
    const auto *kd0_105 = buffer.data(kd0 + 105);
    const auto *kd0_107 = buffer.data(kd0 + 107);
    const auto *kd0_111 = buffer.data(kd0 + 111);
    const auto *kd0_113 = buffer.data(kd0 + 113);
    const auto *kd0_125 = buffer.data(kd0 + 125);
    const auto *kd0_129 = buffer.data(kd0 + 129);

    const auto *kd1_0 = buffer.data(kd1 + 0);
    const auto *kd1_6 = buffer.data(kd1 + 6);
    const auto *kd1_12 = buffer.data(kd1 + 12);
    const auto *kd1_18 = buffer.data(kd1 + 18);
    const auto *kd1_21 = buffer.data(kd1 + 21);
    const auto *kd1_30 = buffer.data(kd1 + 30);
    const auto *kd1_35 = buffer.data(kd1 + 35);
    const auto *kd1_36 = buffer.data(kd1 + 36);
    const auto *kd1_39 = buffer.data(kd1 + 39);
    const auto *kd1_48 = buffer.data(kd1 + 48);
    const auto *kd1_54 = buffer.data(kd1 + 54);
    const auto *kd1_59 = buffer.data(kd1 + 59);
    const auto *kd1_60 = buffer.data(kd1 + 60);
    const auto *kd1_63 = buffer.data(kd1 + 63);
    const auto *kd1_72 = buffer.data(kd1 + 72);
    const auto *kd1_75 = buffer.data(kd1 + 75);
    const auto *kd1_77 = buffer.data(kd1 + 77);
    const auto *kd1_89 = buffer.data(kd1 + 89);
    const auto *kd1_93 = buffer.data(kd1 + 93);
    const auto *kd1_105 = buffer.data(kd1 + 105);
    const auto *kd1_107 = buffer.data(kd1 + 107);
    const auto *kd1_111 = buffer.data(kd1 + 111);
    const auto *kd1_113 = buffer.data(kd1 + 113);
    const auto *kd1_125 = buffer.data(kd1 + 125);
    const auto *kd1_129 = buffer.data(kd1 + 129);

    const auto *lp_0 = buffer.data(lp + 0);
    const auto *lp_1 = buffer.data(lp + 1);
    const auto *lp_2 = buffer.data(lp + 2);
    const auto *lp_4 = buffer.data(lp + 4);
    const auto *lp_8 = buffer.data(lp + 8);
    const auto *lp_10 = buffer.data(lp + 10);
    const auto *lp_11 = buffer.data(lp + 11);
    const auto *lp_14 = buffer.data(lp + 14);
    const auto *lp_16 = buffer.data(lp + 16);
    const auto *lp_17 = buffer.data(lp + 17);
    const auto *lp_19 = buffer.data(lp + 19);
    const auto *lp_20 = buffer.data(lp + 20);
    const auto *lp_23 = buffer.data(lp + 23);
    const auto *lp_25 = buffer.data(lp + 25);
    const auto *lp_26 = buffer.data(lp + 26);
    const auto *lp_28 = buffer.data(lp + 28);
    const auto *lp_29 = buffer.data(lp + 29);
    const auto *lp_31 = buffer.data(lp + 31);
    const auto *lp_32 = buffer.data(lp + 32);
    const auto *lp_35 = buffer.data(lp + 35);
    const auto *lp_37 = buffer.data(lp + 37);
    const auto *lp_38 = buffer.data(lp + 38);
    const auto *lp_40 = buffer.data(lp + 40);
    const auto *lp_41 = buffer.data(lp + 41);
    const auto *lp_43 = buffer.data(lp + 43);
    const auto *lp_44 = buffer.data(lp + 44);
    const auto *lp_46 = buffer.data(lp + 46);
    const auto *lp_47 = buffer.data(lp + 47);
    const auto *lp_50 = buffer.data(lp + 50);
    const auto *lp_52 = buffer.data(lp + 52);
    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_55 = buffer.data(lp + 55);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_58 = buffer.data(lp + 58);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_64 = buffer.data(lp + 64);
    const auto *lp_68 = buffer.data(lp + 68);

    const auto *ld_0 = buffer.data(ld + 0);
    const auto *ld_3 = buffer.data(ld + 3);
    const auto *ld_5 = buffer.data(ld + 5);
    const auto *ld_6 = buffer.data(ld + 6);
    const auto *ld_7 = buffer.data(ld + 7);
    const auto *ld_9 = buffer.data(ld + 9);
    const auto *ld_12 = buffer.data(ld + 12);
    const auto *ld_14 = buffer.data(ld + 14);
    const auto *ld_17 = buffer.data(ld + 17);
    const auto *ld_18 = buffer.data(ld + 18);
    const auto *ld_19 = buffer.data(ld + 19);
    const auto *ld_21 = buffer.data(ld + 21);
    const auto *ld_23 = buffer.data(ld + 23);
    const auto *ld_30 = buffer.data(ld + 30);
    const auto *ld_32 = buffer.data(ld + 32);
    const auto *ld_33 = buffer.data(ld + 33);
    const auto *ld_35 = buffer.data(ld + 35);
    const auto *ld_36 = buffer.data(ld + 36);
    const auto *ld_37 = buffer.data(ld + 37);
    const auto *ld_39 = buffer.data(ld + 39);
    const auto *ld_41 = buffer.data(ld + 41);
    const auto *ld_48 = buffer.data(ld + 48);
    const auto *ld_54 = buffer.data(ld + 54);
    const auto *ld_56 = buffer.data(ld + 56);
    const auto *ld_57 = buffer.data(ld + 57);
    const auto *ld_59 = buffer.data(ld + 59);
    const auto *ld_60 = buffer.data(ld + 60);
    const auto *ld_61 = buffer.data(ld + 61);
    const auto *ld_63 = buffer.data(ld + 63);
    const auto *ld_65 = buffer.data(ld + 65);
    const auto *ld_72 = buffer.data(ld + 72);
    const auto *ld_75 = buffer.data(ld + 75);
    const auto *ld_77 = buffer.data(ld + 77);
    const auto *ld_78 = buffer.data(ld + 78);
    const auto *ld_84 = buffer.data(ld + 84);
    const auto *ld_86 = buffer.data(ld + 86);
    const auto *ld_87 = buffer.data(ld + 87);
    const auto *ld_89 = buffer.data(ld + 89);
    const auto *ld_90 = buffer.data(ld + 90);
    const auto *ld_91 = buffer.data(ld + 91);
    const auto *ld_93 = buffer.data(ld + 93);
    const auto *ld_95 = buffer.data(ld + 95);
    const auto *ld_102 = buffer.data(ld + 102);
    const auto *ld_105 = buffer.data(ld + 105);
    const auto *ld_107 = buffer.data(ld + 107);
    const auto *ld_111 = buffer.data(ld + 111);
    const auto *ld_113 = buffer.data(ld + 113);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_129 = buffer.data(ld + 129);

    const auto *ms0_0 = buffer.data(ms0 + 0);
    const auto *ms0_3 = buffer.data(ms0 + 3);
    const auto *ms0_5 = buffer.data(ms0 + 5);
    const auto *ms0_6 = buffer.data(ms0 + 6);
    const auto *ms0_9 = buffer.data(ms0 + 9);
    const auto *ms0_10 = buffer.data(ms0 + 10);
    const auto *ms0_14 = buffer.data(ms0 + 14);
    const auto *ms0_15 = buffer.data(ms0 + 15);
    const auto *ms0_20 = buffer.data(ms0 + 20);
    const auto *ms0_21 = buffer.data(ms0 + 21);

    const auto *ms1_0 = buffer.data(ms1 + 0);
    const auto *ms1_3 = buffer.data(ms1 + 3);
    const auto *ms1_5 = buffer.data(ms1 + 5);
    const auto *ms1_6 = buffer.data(ms1 + 6);
    const auto *ms1_9 = buffer.data(ms1 + 9);
    const auto *ms1_10 = buffer.data(ms1 + 10);
    const auto *ms1_14 = buffer.data(ms1 + 14);
    const auto *ms1_15 = buffer.data(ms1 + 15);
    const auto *ms1_20 = buffer.data(ms1 + 20);
    const auto *ms1_21 = buffer.data(ms1 + 21);

    const auto *mp_0 = buffer.data(mp + 0);
    const auto *mp_1 = buffer.data(mp + 1);
    const auto *mp_2 = buffer.data(mp + 2);
    const auto *mp_3 = buffer.data(mp + 3);
    const auto *mp_4 = buffer.data(mp + 4);
    const auto *mp_6 = buffer.data(mp + 6);
    const auto *mp_8 = buffer.data(mp + 8);
    const auto *mp_9 = buffer.data(mp + 9);
    const auto *mp_10 = buffer.data(mp + 10);
    const auto *mp_11 = buffer.data(mp + 11);
    const auto *mp_14 = buffer.data(mp + 14);
    const auto *mp_15 = buffer.data(mp + 15);
    const auto *mp_16 = buffer.data(mp + 16);
    const auto *mp_17 = buffer.data(mp + 17);
    const auto *mp_18 = buffer.data(mp + 18);
    const auto *mp_19 = buffer.data(mp + 19);
    const auto *mp_20 = buffer.data(mp + 20);
    const auto *mp_23 = buffer.data(mp + 23);
    const auto *mp_25 = buffer.data(mp + 25);
    const auto *mp_26 = buffer.data(mp + 26);
    const auto *mp_27 = buffer.data(mp + 27);
    const auto *mp_28 = buffer.data(mp + 28);
    const auto *mp_29 = buffer.data(mp + 29);
    const auto *mp_30 = buffer.data(mp + 30);
    const auto *mp_31 = buffer.data(mp + 31);
    const auto *mp_32 = buffer.data(mp + 32);
    const auto *mp_35 = buffer.data(mp + 35);
    const auto *mp_37 = buffer.data(mp + 37);
    const auto *mp_38 = buffer.data(mp + 38);
    const auto *mp_40 = buffer.data(mp + 40);
    const auto *mp_41 = buffer.data(mp + 41);
    const auto *mp_42 = buffer.data(mp + 42);
    const auto *mp_43 = buffer.data(mp + 43);
    const auto *mp_44 = buffer.data(mp + 44);
    const auto *mp_45 = buffer.data(mp + 45);
    const auto *mp_46 = buffer.data(mp + 46);
    const auto *mp_47 = buffer.data(mp + 47);
    const auto *mp_50 = buffer.data(mp + 50);
    const auto *mp_52 = buffer.data(mp + 52);
    const auto *mp_53 = buffer.data(mp + 53);
    const auto *mp_55 = buffer.data(mp + 55);
    const auto *mp_56 = buffer.data(mp + 56);
    const auto *mp_58 = buffer.data(mp + 58);
    const auto *mp_59 = buffer.data(mp + 59);
    const auto *mp_60 = buffer.data(mp + 60);
    const auto *mp_61 = buffer.data(mp + 61);
    const auto *mp_62 = buffer.data(mp + 62);
    const auto *mp_63 = buffer.data(mp + 63);
    const auto *mp_64 = buffer.data(mp + 64);
    const auto *mp_65 = buffer.data(mp + 65);
    const auto *mp_68 = buffer.data(mp + 68);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, t_4, t_5, pb_x, pb_y, pb_z, lp_0, ms0_0, ms1_0, \
                         mp_0, mp_1, mp_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * lp_0[k]
                 + f_1 * ms0_0[k]
                 - f_2 * ms1_0[k]
                 + pb_x[k] * mp_0[k];

        t_1[k] = pb_y[k] * mp_0[k];

        t_2[k] = pb_z[k] * mp_0[k];

        t_3[k] = f_1 * ms0_0[k]
                 - f_2 * ms1_0[k]
                 + pb_y[k] * mp_1[k];

        t_4[k] = pb_y[k] * mp_2[k];

        t_5[k] = f_1 * ms0_0[k]
                 - f_2 * ms1_0[k]
                 + pb_z[k] * mp_2[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, t_10, t_11, pa_y, pb_x, pb_z, lp_1, lp_4, ld_0, \
                         ld_3, ld_5, mp_3, mp_4 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = pa_y[k] * ld_0[k];

        t_7[k] = f_3 * lp_4[k]
                 + pb_x[k] * mp_4[k];

        t_8[k] = pb_z[k] * mp_3[k];

        t_9[k] = f_4 * lp_1[k]
                 + pa_y[k] * ld_3[k];

        t_10[k] = pb_z[k] * mp_4[k];

        t_11[k] = pa_y[k] * ld_5[k];
    }

#pragma omp simd aligned(t_12, t_13, t_14, t_15, t_16, t_17, pa_z, pb_x, pb_y, lp_2, lp_8, \
                         ld_0, ld_3, ld_5, mp_6, mp_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_12[k] = pa_z[k] * ld_0[k];

        t_13[k] = pb_y[k] * mp_6[k];

        t_14[k] = f_3 * lp_8[k]
                  + pb_x[k] * mp_8[k];

        t_15[k] = pa_z[k] * ld_3[k];

        t_16[k] = pb_y[k] * mp_8[k];

        t_17[k] = f_4 * lp_2[k]
                  + pa_z[k] * ld_5[k];
    }

#pragma omp simd aligned(t_18, t_19, t_20, pa_y, pb_x, pb_z, kd0_0, kd1_0, lp_10, ld_6, mp_9, \
                         mp_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_18[k] = f_5 * kd0_0[k]
                  - f_6 * kd1_0[k]
                  + pa_y[k] * ld_6[k];

        t_19[k] = f_7 * lp_10[k]
                  + pb_x[k] * mp_10[k];

        t_20[k] = pb_z[k] * mp_9[k];
    }

#pragma omp simd aligned(t_21, t_22, t_23, t_24, pa_x, pa_y, pb_z, kd0_21, kd1_21, ld_12, \
                         ld_21, ms0_3, ms1_3, mp_10, mp_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_8 * kd0_21[k]
                  - f_9 * kd1_21[k]
                  + pa_x[k] * ld_21[k];

        t_22[k] = pb_z[k] * mp_10[k];

        t_23[k] = f_1 * ms0_3[k]
                  - f_2 * ms1_3[k]
                  + pb_z[k] * mp_11[k];

        t_24[k] = pa_y[k] * ld_12[k];
    }

#pragma omp simd aligned(t_25, t_26, t_27, t_28, t_29, pa_y, pa_z, pb_y, lp_8, ld_7, ld_9, \
                         ld_14, ld_17, mp_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = pa_z[k] * ld_7[k];

        t_26[k] = pa_y[k] * ld_14[k];

        t_27[k] = pa_z[k] * ld_9[k];

        t_28[k] = f_10 * lp_8[k]
                  + pb_y[k] * mp_14[k];

        t_29[k] = pa_y[k] * ld_17[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, t_33, pa_z, pb_x, pb_y, kd0_0, kd1_0, lp_17, ld_12, \
                         ms0_5, ms1_5, mp_15, mp_16, mp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = f_5 * kd0_0[k]
                  - f_6 * kd1_0[k]
                  + pa_z[k] * ld_12[k];

        t_31[k] = pb_y[k] * mp_15[k];

        t_32[k] = f_7 * lp_17[k]
                  + pb_x[k] * mp_17[k];

        t_33[k] = f_1 * ms0_5[k]
                  - f_2 * ms1_5[k]
                  + pb_y[k] * mp_16[k];
    }

#pragma omp simd aligned(t_34, t_35, t_36, pa_x, pa_y, pb_y, kd0_6, kd0_35, kd1_6, kd1_35, \
                         ld_18, ld_35, mp_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_34[k] = pb_y[k] * mp_17[k];

        t_35[k] = f_8 * kd0_35[k]
                  - f_9 * kd1_35[k]
                  + pa_x[k] * ld_35[k];

        t_36[k] = f_11 * kd0_6[k]
                  - f_12 * kd1_6[k]
                  + pa_y[k] * ld_18[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, t_40, pa_x, pb_x, pb_z, kd0_39, kd1_39, lp_19, \
                         ld_39, mp_18, mp_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_13 * lp_19[k]
                  + pb_x[k] * mp_19[k];

        t_38[k] = pb_z[k] * mp_18[k];

        t_39[k] = f_14 * kd0_39[k]
                  - f_15 * kd1_39[k]
                  + pa_x[k] * ld_39[k];

        t_40[k] = pb_z[k] * mp_19[k];
    }

#pragma omp simd aligned(t_41, t_42, t_43, t_44, t_45, pa_z, pb_x, pb_z, lp_23, ld_18, ld_19, \
                         ld_21, ms0_6, ms1_6, mp_20, mp_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_1 * ms0_6[k]
                  - f_2 * ms1_6[k]
                  + pb_z[k] * mp_20[k];

        t_42[k] = pa_z[k] * ld_18[k];

        t_43[k] = pa_z[k] * ld_19[k];

        t_44[k] = f_13 * lp_23[k]
                  + pb_x[k] * mp_23[k];

        t_45[k] = pa_z[k] * ld_21[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_y, pa_z, pb_x, pb_y, lp_11, lp_14, lp_25, \
                         ld_23, ld_30, mp_23, mp_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_4 * lp_14[k]
                  + pb_y[k] * mp_23[k];

        t_47[k] = f_4 * lp_11[k]
                  + pa_z[k] * ld_23[k];

        t_48[k] = pa_y[k] * ld_30[k];

        t_49[k] = f_13 * lp_25[k]
                  + pb_x[k] * mp_25[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, t_53, pa_y, pb_y, lp_16, lp_17, ld_32, ld_33, \
                         ld_35, mp_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = pa_y[k] * ld_32[k];

        t_51[k] = f_4 * lp_16[k]
                  + pa_y[k] * ld_33[k];

        t_52[k] = f_10 * lp_17[k]
                  + pb_y[k] * mp_26[k];

        t_53[k] = pa_y[k] * ld_35[k];
    }

#pragma omp simd aligned(t_54, t_55, t_56, t_57, pa_z, pb_x, pb_y, kd0_12, kd1_12, lp_29, \
                         ld_30, ms0_9, ms1_9, mp_27, mp_28, mp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_54[k] = f_11 * kd0_12[k]
                  - f_12 * kd1_12[k]
                  + pa_z[k] * ld_30[k];

        t_55[k] = pb_y[k] * mp_27[k];

        t_56[k] = f_13 * lp_29[k]
                  + pb_x[k] * mp_29[k];

        t_57[k] = f_1 * ms0_9[k]
                  - f_2 * ms1_9[k]
                  + pb_y[k] * mp_28[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, kd0_18, kd0_59, kd1_18, kd1_59, \
                         ld_36, ld_59, mp_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = pb_y[k] * mp_29[k];

        t_59[k] = f_14 * kd0_59[k]
                  - f_15 * kd1_59[k]
                  + pa_x[k] * ld_59[k];

        t_60[k] = f_16 * kd0_18[k]
                  - f_17 * kd1_18[k]
                  + pa_y[k] * ld_36[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_x, pb_x, pb_z, kd0_63, kd1_63, lp_31, \
                         ld_63, mp_30, mp_31 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_18 * lp_31[k]
                  + pb_x[k] * mp_31[k];

        t_62[k] = pb_z[k] * mp_30[k];

        t_63[k] = f_19 * kd0_63[k]
                  - f_20 * kd1_63[k]
                  + pa_x[k] * ld_63[k];

        t_64[k] = pb_z[k] * mp_31[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, t_69, pa_z, pb_x, pb_z, lp_35, ld_36, ld_37, \
                         ld_39, ms0_10, ms1_10, mp_32, mp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_1 * ms0_10[k]
                  - f_2 * ms1_10[k]
                  + pb_z[k] * mp_32[k];

        t_66[k] = pa_z[k] * ld_36[k];

        t_67[k] = pa_z[k] * ld_37[k];

        t_68[k] = f_18 * lp_35[k]
                  + pb_x[k] * mp_35[k];

        t_69[k] = pa_z[k] * ld_39[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pa_y, pa_z, pb_y, kd0_30, kd1_30, lp_20, lp_23, \
                         ld_41, ld_48, mp_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_21 * lp_23[k]
                  + pb_y[k] * mp_35[k];

        t_71[k] = f_4 * lp_20[k]
                  + pa_z[k] * ld_41[k];

        t_72[k] = f_5 * kd0_30[k]
                  - f_6 * kd1_30[k]
                  + pa_y[k] * ld_48[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pb_x, pb_y, kd0_75, kd1_75, lp_26, \
                         lp_37, lp_38, ld_75, mp_37, mp_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_18 * lp_37[k]
                  + pb_x[k] * mp_37[k];

        t_74[k] = f_18 * lp_38[k]
                  + pb_x[k] * mp_38[k];

        t_75[k] = f_19 * kd0_75[k]
                  - f_20 * kd1_75[k]
                  + pa_x[k] * ld_75[k];

        t_76[k] = f_4 * lp_26[k]
                  + pb_y[k] * mp_38[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, t_80, pa_x, pa_y, pb_x, kd0_77, kd1_77, lp_40, \
                         ld_54, ld_56, ld_77, mp_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_19 * kd0_77[k]
                  - f_20 * kd1_77[k]
                  + pa_x[k] * ld_77[k];

        t_78[k] = pa_y[k] * ld_54[k];

        t_79[k] = f_18 * lp_40[k]
                  + pb_x[k] * mp_40[k];

        t_80[k] = pa_y[k] * ld_56[k];
    }

#pragma omp simd aligned(t_81, t_82, t_83, t_84, pa_y, pa_z, pb_y, kd0_30, kd1_30, lp_28, \
                         lp_29, ld_54, ld_57, ld_59, mp_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_81[k] = f_4 * lp_28[k]
                  + pa_y[k] * ld_57[k];

        t_82[k] = f_10 * lp_29[k]
                  + pb_y[k] * mp_41[k];

        t_83[k] = pa_y[k] * ld_59[k];

        t_84[k] = f_16 * kd0_30[k]
                  - f_17 * kd1_30[k]
                  + pa_z[k] * ld_54[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pb_x, pb_y, lp_44, ms0_14, ms1_14, mp_42, \
                         mp_43, mp_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pb_y[k] * mp_42[k];

        t_86[k] = f_18 * lp_44[k]
                  + pb_x[k] * mp_44[k];

        t_87[k] = f_1 * ms0_14[k]
                  - f_2 * ms1_14[k]
                  + pb_y[k] * mp_43[k];

        t_88[k] = pb_y[k] * mp_44[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pa_y, pb_x, kd0_36, kd0_89, kd1_36, kd1_89, \
                         lp_46, ld_60, ld_89, mp_46 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = f_19 * kd0_89[k]
                  - f_20 * kd1_89[k]
                  + pa_x[k] * ld_89[k];

        t_90[k] = f_19 * kd0_36[k]
                  - f_20 * kd1_36[k]
                  + pa_y[k] * ld_60[k];

        t_91[k] = f_22 * lp_46[k]
                  + pb_x[k] * mp_46[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, kd0_93, kd1_93, ld_93, ms0_15, \
                         ms1_15, mp_45, mp_46, mp_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = pb_z[k] * mp_45[k];

        t_93[k] = f_16 * kd0_93[k]
                  - f_17 * kd1_93[k]
                  + pa_x[k] * ld_93[k];

        t_94[k] = pb_z[k] * mp_46[k];

        t_95[k] = f_1 * ms0_15[k]
                  - f_2 * ms1_15[k]
                  + pb_z[k] * mp_47[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, t_99, t_100, pa_z, pb_x, pb_y, lp_35, lp_50, ld_60, \
                         ld_61, ld_63, mp_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = pa_z[k] * ld_60[k];

        t_97[k] = pa_z[k] * ld_61[k];

        t_98[k] = f_22 * lp_50[k]
                  + pb_x[k] * mp_50[k];

        t_99[k] = pa_z[k] * ld_63[k];

        t_100[k] = f_22 * lp_35[k]
                   + pb_y[k] * mp_50[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pa_y, pa_z, pb_x, kd0_48, kd1_48, lp_32, \
                         lp_52, lp_53, ld_65, ld_72, mp_52, mp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_4 * lp_32[k]
                   + pa_z[k] * ld_65[k];

        t_102[k] = f_11 * kd0_48[k]
                   - f_12 * kd1_48[k]
                   + pa_y[k] * ld_72[k];

        t_103[k] = f_22 * lp_52[k]
                   + pb_x[k] * mp_52[k];

        t_104[k] = f_22 * lp_53[k]
                   + pb_x[k] * mp_53[k];
    }

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_y, kd0_105, kd0_107, kd1_105, kd1_107, \
                         lp_38, ld_105, ld_107, mp_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_16 * kd0_105[k]
                   - f_17 * kd1_105[k]
                   + pa_x[k] * ld_105[k];

        t_106[k] = f_21 * lp_38[k]
                   + pb_y[k] * mp_53[k];

        t_107[k] = f_16 * kd0_107[k]
                   - f_17 * kd1_107[k]
                   + pa_x[k] * ld_107[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_y, pb_x, kd0_54, kd1_54, lp_55, lp_56, ld_78, \
                         mp_55, mp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_5 * kd0_54[k]
                   - f_6 * kd1_54[k]
                   + pa_y[k] * ld_78[k];

        t_109[k] = f_22 * lp_55[k]
                   + pb_x[k] * mp_55[k];

        t_110[k] = f_22 * lp_56[k]
                   + pb_x[k] * mp_56[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, t_114, pa_x, pa_y, pb_y, kd0_111, kd0_113, \
                         kd1_111, kd1_113, lp_41, ld_84, ld_111, ld_113, \
                         mp_56 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = f_16 * kd0_111[k]
                   - f_17 * kd1_111[k]
                   + pa_x[k] * ld_111[k];

        t_112[k] = f_4 * lp_41[k]
                   + pb_y[k] * mp_56[k];

        t_113[k] = f_16 * kd0_113[k]
                   - f_17 * kd1_113[k]
                   + pa_x[k] * ld_113[k];

        t_114[k] = pa_y[k] * ld_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, t_118, t_119, pa_y, pb_x, pb_y, lp_43, lp_44, \
                         lp_58, ld_86, ld_87, ld_89, mp_58, mp_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_22 * lp_58[k]
                   + pb_x[k] * mp_58[k];

        t_116[k] = pa_y[k] * ld_86[k];

        t_117[k] = f_4 * lp_43[k]
                   + pa_y[k] * ld_87[k];

        t_118[k] = f_10 * lp_44[k]
                   + pb_y[k] * mp_59[k];

        t_119[k] = pa_y[k] * ld_89[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, t_123, pa_z, pb_x, pb_y, kd0_54, kd1_54, lp_62, \
                         ld_84, ms0_20, ms1_20, mp_60, mp_61, mp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_19 * kd0_54[k]
                   - f_20 * kd1_54[k]
                   + pa_z[k] * ld_84[k];

        t_121[k] = pb_y[k] * mp_60[k];

        t_122[k] = f_22 * lp_62[k]
                   + pb_x[k] * mp_62[k];

        t_123[k] = f_1 * ms0_20[k]
                   - f_2 * ms1_20[k]
                   + pb_y[k] * mp_61[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pa_y, pb_y, kd0_60, kd0_125, kd1_60, \
                         kd1_125, ld_90, ld_125, mp_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = pb_y[k] * mp_62[k];

        t_125[k] = f_16 * kd0_125[k]
                   - f_17 * kd1_125[k]
                   + pa_x[k] * ld_125[k];

        t_126[k] = f_14 * kd0_60[k]
                   - f_15 * kd1_60[k]
                   + pa_y[k] * ld_90[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_x, pb_z, kd0_129, kd1_129, \
                         lp_64, ld_129, mp_63, mp_64 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_21 * lp_64[k]
                   + pb_x[k] * mp_64[k];

        t_128[k] = pb_z[k] * mp_63[k];

        t_129[k] = f_11 * kd0_129[k]
                   - f_12 * kd1_129[k]
                   + pa_x[k] * ld_129[k];

        t_130[k] = pb_z[k] * mp_64[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, t_134, t_135, pa_z, pb_x, pb_z, lp_68, ld_90, \
                         ld_91, ld_93, ms0_21, ms1_21, mp_65, mp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_1 * ms0_21[k]
                   - f_2 * ms1_21[k]
                   + pb_z[k] * mp_65[k];

        t_132[k] = pa_z[k] * ld_90[k];

        t_133[k] = pa_z[k] * ld_91[k];

        t_134[k] = f_21 * lp_68[k]
                   + pb_x[k] * mp_68[k];

        t_135[k] = pa_z[k] * ld_93[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_y, pa_z, pb_y, kd0_72, kd1_72, lp_47, lp_50, \
                         ld_95, ld_102, mp_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_18 * lp_50[k]
                   + pb_y[k] * mp_68[k];

        t_137[k] = f_4 * lp_47[k]
                   + pa_z[k] * ld_95[k];

        t_138[k] = f_16 * kd0_72[k]
                   - f_17 * kd1_72[k]
                   + pa_y[k] * ld_102[k];
    }
}

static auto
compute_prim_md_electron_repulsion_0_piece1(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kd0,
                                            const size_t kd1, const size_t lp, const size_t ld,
                                            const size_t ms0, const size_t ms1, const size_t mp,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_3 = 4.0 / p;
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_8 = 3.0 / alpha;
    const auto f_9 = 3.0 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 3.0 / p;
    const auto f_14 = 2.5 / alpha;
    const auto f_15 = 2.5 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 2.5 / p;
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.5 / p;
    const auto f_22 = 2.0 / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kd0_78 = buffer.data(kd0 + 78);
    const auto *kd0_84 = buffer.data(kd0 + 84);
    const auto *kd0_90 = buffer.data(kd0 + 90);
    const auto *kd0_102 = buffer.data(kd0 + 102);
    const auto *kd0_108 = buffer.data(kd0 + 108);
    const auto *kd0_114 = buffer.data(kd0 + 114);
    const auto *kd0_120 = buffer.data(kd0 + 120);
    const auto *kd0_141 = buffer.data(kd0 + 141);
    const auto *kd0_143 = buffer.data(kd0 + 143);
    const auto *kd0_147 = buffer.data(kd0 + 147);
    const auto *kd0_149 = buffer.data(kd0 + 149);
    const auto *kd0_153 = buffer.data(kd0 + 153);
    const auto *kd0_155 = buffer.data(kd0 + 155);
    const auto *kd0_167 = buffer.data(kd0 + 167);
    const auto *kd0_171 = buffer.data(kd0 + 171);
    const auto *kd0_183 = buffer.data(kd0 + 183);
    const auto *kd0_185 = buffer.data(kd0 + 185);
    const auto *kd0_189 = buffer.data(kd0 + 189);
    const auto *kd0_191 = buffer.data(kd0 + 191);
    const auto *kd0_195 = buffer.data(kd0 + 195);
    const auto *kd0_197 = buffer.data(kd0 + 197);
    const auto *kd0_201 = buffer.data(kd0 + 201);
    const auto *kd0_203 = buffer.data(kd0 + 203);
    const auto *kd0_215 = buffer.data(kd0 + 215);

    const auto *kd1_78 = buffer.data(kd1 + 78);
    const auto *kd1_84 = buffer.data(kd1 + 84);
    const auto *kd1_90 = buffer.data(kd1 + 90);
    const auto *kd1_102 = buffer.data(kd1 + 102);
    const auto *kd1_108 = buffer.data(kd1 + 108);
    const auto *kd1_114 = buffer.data(kd1 + 114);
    const auto *kd1_120 = buffer.data(kd1 + 120);
    const auto *kd1_141 = buffer.data(kd1 + 141);
    const auto *kd1_143 = buffer.data(kd1 + 143);
    const auto *kd1_147 = buffer.data(kd1 + 147);
    const auto *kd1_149 = buffer.data(kd1 + 149);
    const auto *kd1_153 = buffer.data(kd1 + 153);
    const auto *kd1_155 = buffer.data(kd1 + 155);
    const auto *kd1_167 = buffer.data(kd1 + 167);
    const auto *kd1_171 = buffer.data(kd1 + 171);
    const auto *kd1_183 = buffer.data(kd1 + 183);
    const auto *kd1_185 = buffer.data(kd1 + 185);
    const auto *kd1_189 = buffer.data(kd1 + 189);
    const auto *kd1_191 = buffer.data(kd1 + 191);
    const auto *kd1_195 = buffer.data(kd1 + 195);
    const auto *kd1_197 = buffer.data(kd1 + 197);
    const auto *kd1_201 = buffer.data(kd1 + 201);
    const auto *kd1_203 = buffer.data(kd1 + 203);
    const auto *kd1_215 = buffer.data(kd1 + 215);

    const auto *lp_53 = buffer.data(lp + 53);
    const auto *lp_56 = buffer.data(lp + 56);
    const auto *lp_59 = buffer.data(lp + 59);
    const auto *lp_61 = buffer.data(lp + 61);
    const auto *lp_62 = buffer.data(lp + 62);
    const auto *lp_65 = buffer.data(lp + 65);
    const auto *lp_68 = buffer.data(lp + 68);
    const auto *lp_70 = buffer.data(lp + 70);
    const auto *lp_71 = buffer.data(lp + 71);
    const auto *lp_73 = buffer.data(lp + 73);
    const auto *lp_74 = buffer.data(lp + 74);
    const auto *lp_76 = buffer.data(lp + 76);
    const auto *lp_77 = buffer.data(lp + 77);
    const auto *lp_79 = buffer.data(lp + 79);
    const auto *lp_80 = buffer.data(lp + 80);
    const auto *lp_82 = buffer.data(lp + 82);
    const auto *lp_83 = buffer.data(lp + 83);
    const auto *lp_85 = buffer.data(lp + 85);
    const auto *lp_89 = buffer.data(lp + 89);
    const auto *lp_91 = buffer.data(lp + 91);
    const auto *lp_92 = buffer.data(lp + 92);
    const auto *lp_94 = buffer.data(lp + 94);
    const auto *lp_95 = buffer.data(lp + 95);
    const auto *lp_97 = buffer.data(lp + 97);
    const auto *lp_98 = buffer.data(lp + 98);
    const auto *lp_100 = buffer.data(lp + 100);
    const auto *lp_101 = buffer.data(lp + 101);
    const auto *lp_103 = buffer.data(lp + 103);
    const auto *lp_107 = buffer.data(lp + 107);
    const auto *lp_108 = buffer.data(lp + 108);
    const auto *lp_109 = buffer.data(lp + 109);
    const auto *lp_110 = buffer.data(lp + 110);
    const auto *lp_113 = buffer.data(lp + 113);
    const auto *lp_114 = buffer.data(lp + 114);
    const auto *lp_115 = buffer.data(lp + 115);
    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_117 = buffer.data(lp + 117);
    const auto *lp_118 = buffer.data(lp + 118);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_120 = buffer.data(lp + 120);
    const auto *lp_121 = buffer.data(lp + 121);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_123 = buffer.data(lp + 123);
    const auto *lp_124 = buffer.data(lp + 124);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_126 = buffer.data(lp + 126);
    const auto *lp_127 = buffer.data(lp + 127);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_130 = buffer.data(lp + 130);
    const auto *lp_132 = buffer.data(lp + 132);
    const auto *lp_134 = buffer.data(lp + 134);

    const auto *ld_108 = buffer.data(ld + 108);
    const auto *ld_114 = buffer.data(ld + 114);
    const auto *ld_120 = buffer.data(ld + 120);
    const auto *ld_122 = buffer.data(ld + 122);
    const auto *ld_123 = buffer.data(ld + 123);
    const auto *ld_125 = buffer.data(ld + 125);
    const auto *ld_126 = buffer.data(ld + 126);
    const auto *ld_127 = buffer.data(ld + 127);
    const auto *ld_129 = buffer.data(ld + 129);
    const auto *ld_131 = buffer.data(ld + 131);
    const auto *ld_138 = buffer.data(ld + 138);
    const auto *ld_141 = buffer.data(ld + 141);
    const auto *ld_143 = buffer.data(ld + 143);
    const auto *ld_144 = buffer.data(ld + 144);
    const auto *ld_147 = buffer.data(ld + 147);
    const auto *ld_149 = buffer.data(ld + 149);
    const auto *ld_150 = buffer.data(ld + 150);
    const auto *ld_153 = buffer.data(ld + 153);
    const auto *ld_155 = buffer.data(ld + 155);
    const auto *ld_156 = buffer.data(ld + 156);
    const auto *ld_162 = buffer.data(ld + 162);
    const auto *ld_164 = buffer.data(ld + 164);
    const auto *ld_165 = buffer.data(ld + 165);
    const auto *ld_167 = buffer.data(ld + 167);
    const auto *ld_168 = buffer.data(ld + 168);
    const auto *ld_169 = buffer.data(ld + 169);
    const auto *ld_171 = buffer.data(ld + 171);
    const auto *ld_183 = buffer.data(ld + 183);
    const auto *ld_185 = buffer.data(ld + 185);
    const auto *ld_189 = buffer.data(ld + 189);
    const auto *ld_191 = buffer.data(ld + 191);
    const auto *ld_195 = buffer.data(ld + 195);
    const auto *ld_197 = buffer.data(ld + 197);
    const auto *ld_201 = buffer.data(ld + 201);
    const auto *ld_203 = buffer.data(ld + 203);
    const auto *ld_210 = buffer.data(ld + 210);
    const auto *ld_212 = buffer.data(ld + 212);
    const auto *ld_215 = buffer.data(ld + 215);
    const auto *ld_216 = buffer.data(ld + 216);
    const auto *ld_219 = buffer.data(ld + 219);
    const auto *ld_221 = buffer.data(ld + 221);
    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_226 = buffer.data(ld + 226);
    const auto *ld_227 = buffer.data(ld + 227);
    const auto *ld_228 = buffer.data(ld + 228);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_232 = buffer.data(ld + 232);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_234 = buffer.data(ld + 234);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_238 = buffer.data(ld + 238);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_240 = buffer.data(ld + 240);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_244 = buffer.data(ld + 244);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_246 = buffer.data(ld + 246);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_250 = buffer.data(ld + 250);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_252 = buffer.data(ld + 252);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_256 = buffer.data(ld + 256);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_261 = buffer.data(ld + 261);
    const auto *ld_262 = buffer.data(ld + 262);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_269 = buffer.data(ld + 269);

    const auto *ms0_27 = buffer.data(ms0 + 27);
    const auto *ms0_28 = buffer.data(ms0 + 28);
    const auto *ms0_35 = buffer.data(ms0 + 35);
    const auto *ms0_45 = buffer.data(ms0 + 45);

    const auto *ms1_27 = buffer.data(ms1 + 27);
    const auto *ms1_28 = buffer.data(ms1 + 28);
    const auto *ms1_35 = buffer.data(ms1 + 35);
    const auto *ms1_45 = buffer.data(ms1 + 45);

    const auto *mp_70 = buffer.data(mp + 70);
    const auto *mp_71 = buffer.data(mp + 71);
    const auto *mp_73 = buffer.data(mp + 73);
    const auto *mp_74 = buffer.data(mp + 74);
    const auto *mp_76 = buffer.data(mp + 76);
    const auto *mp_77 = buffer.data(mp + 77);
    const auto *mp_79 = buffer.data(mp + 79);
    const auto *mp_80 = buffer.data(mp + 80);
    const auto *mp_81 = buffer.data(mp + 81);
    const auto *mp_82 = buffer.data(mp + 82);
    const auto *mp_83 = buffer.data(mp + 83);
    const auto *mp_84 = buffer.data(mp + 84);
    const auto *mp_85 = buffer.data(mp + 85);
    const auto *mp_86 = buffer.data(mp + 86);
    const auto *mp_89 = buffer.data(mp + 89);
    const auto *mp_91 = buffer.data(mp + 91);
    const auto *mp_92 = buffer.data(mp + 92);
    const auto *mp_94 = buffer.data(mp + 94);
    const auto *mp_95 = buffer.data(mp + 95);
    const auto *mp_97 = buffer.data(mp + 97);
    const auto *mp_98 = buffer.data(mp + 98);
    const auto *mp_100 = buffer.data(mp + 100);
    const auto *mp_101 = buffer.data(mp + 101);
    const auto *mp_103 = buffer.data(mp + 103);
    const auto *mp_104 = buffer.data(mp + 104);
    const auto *mp_105 = buffer.data(mp + 105);
    const auto *mp_106 = buffer.data(mp + 106);
    const auto *mp_107 = buffer.data(mp + 107);
    const auto *mp_108 = buffer.data(mp + 108);
    const auto *mp_109 = buffer.data(mp + 109);
    const auto *mp_113 = buffer.data(mp + 113);
    const auto *mp_115 = buffer.data(mp + 115);
    const auto *mp_116 = buffer.data(mp + 116);
    const auto *mp_118 = buffer.data(mp + 118);
    const auto *mp_119 = buffer.data(mp + 119);
    const auto *mp_121 = buffer.data(mp + 121);
    const auto *mp_122 = buffer.data(mp + 122);
    const auto *mp_124 = buffer.data(mp + 124);
    const auto *mp_125 = buffer.data(mp + 125);
    const auto *mp_127 = buffer.data(mp + 127);
    const auto *mp_128 = buffer.data(mp + 128);
    const auto *mp_130 = buffer.data(mp + 130);
    const auto *mp_132 = buffer.data(mp + 132);
    const auto *mp_134 = buffer.data(mp + 134);
    const auto *mp_135 = buffer.data(mp + 135);
    const auto *mp_136 = buffer.data(mp + 136);
    const auto *mp_137 = buffer.data(mp + 137);
    const auto *mp_139 = buffer.data(mp + 139);
    const auto *mp_140 = buffer.data(mp + 140);

#pragma omp simd aligned(t_139, t_140, t_141, t_142, pa_x, pb_x, pb_y, kd0_141, kd1_141, \
                         lp_53, lp_70, lp_71, ld_141, mp_70, mp_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_21 * lp_70[k]
                   + pb_x[k] * mp_70[k];

        t_140[k] = f_21 * lp_71[k]
                   + pb_x[k] * mp_71[k];

        t_141[k] = f_11 * kd0_141[k]
                   - f_12 * kd1_141[k]
                   + pa_x[k] * ld_141[k];

        t_142[k] = f_22 * lp_53[k]
                   + pb_y[k] * mp_71[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pa_x, pa_y, pb_x, kd0_78, kd0_143, kd1_78, \
                         kd1_143, lp_73, ld_108, ld_143, mp_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_11 * kd0_143[k]
                   - f_12 * kd1_143[k]
                   + pa_x[k] * ld_143[k];

        t_144[k] = f_11 * kd0_78[k]
                   - f_12 * kd1_78[k]
                   + pa_y[k] * ld_108[k];

        t_145[k] = f_21 * lp_73[k]
                   + pb_x[k] * mp_73[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_x, pb_y, kd0_147, kd0_149, \
                         kd1_147, kd1_149, lp_56, lp_74, ld_147, ld_149, \
                         mp_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = f_21 * lp_74[k]
                   + pb_x[k] * mp_74[k];

        t_147[k] = f_11 * kd0_147[k]
                   - f_12 * kd1_147[k]
                   + pa_x[k] * ld_147[k];

        t_148[k] = f_21 * lp_56[k]
                   + pb_y[k] * mp_74[k];

        t_149[k] = f_11 * kd0_149[k]
                   - f_12 * kd1_149[k]
                   + pa_x[k] * ld_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pa_y, pb_x, kd0_84, kd1_84, lp_76, lp_77, \
                         ld_114, mp_76, mp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = f_5 * kd0_84[k]
                   - f_6 * kd1_84[k]
                   + pa_y[k] * ld_114[k];

        t_151[k] = f_21 * lp_76[k]
                   + pb_x[k] * mp_76[k];

        t_152[k] = f_21 * lp_77[k]
                   + pb_x[k] * mp_77[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pa_x, pa_y, pb_y, kd0_153, kd0_155, \
                         kd1_153, kd1_155, lp_59, ld_120, ld_153, ld_155, \
                         mp_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_11 * kd0_153[k]
                   - f_12 * kd1_153[k]
                   + pa_x[k] * ld_153[k];

        t_154[k] = f_4 * lp_59[k]
                   + pb_y[k] * mp_77[k];

        t_155[k] = f_11 * kd0_155[k]
                   - f_12 * kd1_155[k]
                   + pa_x[k] * ld_155[k];

        t_156[k] = pa_y[k] * ld_120[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, t_160, t_161, pa_y, pb_x, pb_y, lp_61, lp_62, \
                         lp_79, ld_122, ld_123, ld_125, mp_79, mp_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_21 * lp_79[k]
                   + pb_x[k] * mp_79[k];

        t_158[k] = pa_y[k] * ld_122[k];

        t_159[k] = f_4 * lp_61[k]
                   + pa_y[k] * ld_123[k];

        t_160[k] = f_10 * lp_62[k]
                   + pb_y[k] * mp_80[k];

        t_161[k] = pa_y[k] * ld_125[k];
    }

#pragma omp simd aligned(t_162, t_163, t_164, t_165, pa_z, pb_x, pb_y, kd0_84, kd1_84, lp_83, \
                         ld_120, ms0_27, ms1_27, mp_81, mp_82, mp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_162[k] = f_14 * kd0_84[k]
                   - f_15 * kd1_84[k]
                   + pa_z[k] * ld_120[k];

        t_163[k] = pb_y[k] * mp_81[k];

        t_164[k] = f_21 * lp_83[k]
                   + pb_x[k] * mp_83[k];

        t_165[k] = f_1 * ms0_27[k]
                   - f_2 * ms1_27[k]
                   + pb_y[k] * mp_82[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_x, pa_y, pb_y, kd0_90, kd0_167, kd1_90, \
                         kd1_167, ld_126, ld_167, mp_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pb_y[k] * mp_83[k];

        t_167[k] = f_11 * kd0_167[k]
                   - f_12 * kd1_167[k]
                   + pa_x[k] * ld_167[k];

        t_168[k] = f_8 * kd0_90[k]
                   - f_9 * kd1_90[k]
                   + pa_y[k] * ld_126[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, t_172, pa_x, pb_x, pb_z, kd0_171, kd1_171, \
                         lp_85, ld_171, mp_84, mp_85 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_4 * lp_85[k]
                   + pb_x[k] * mp_85[k];

        t_170[k] = pb_z[k] * mp_84[k];

        t_171[k] = f_5 * kd0_171[k]
                   - f_6 * kd1_171[k]
                   + pa_x[k] * ld_171[k];

        t_172[k] = pb_z[k] * mp_85[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, t_177, pa_z, pb_x, pb_z, lp_89, ld_126, \
                         ld_127, ld_129, ms0_28, ms1_28, mp_86, mp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_1 * ms0_28[k]
                   - f_2 * ms1_28[k]
                   + pb_z[k] * mp_86[k];

        t_174[k] = pa_z[k] * ld_126[k];

        t_175[k] = pa_z[k] * ld_127[k];

        t_176[k] = f_4 * lp_89[k]
                   + pb_x[k] * mp_89[k];

        t_177[k] = pa_z[k] * ld_129[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pa_y, pa_z, pb_y, kd0_102, kd1_102, lp_65, \
                         lp_68, ld_131, ld_138, mp_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_13 * lp_68[k]
                   + pb_y[k] * mp_89[k];

        t_179[k] = f_4 * lp_65[k]
                   + pa_z[k] * ld_131[k];

        t_180[k] = f_19 * kd0_102[k]
                   - f_20 * kd1_102[k]
                   + pa_y[k] * ld_138[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, t_184, pa_x, pb_x, pb_y, kd0_183, kd1_183, \
                         lp_71, lp_91, lp_92, ld_183, mp_91, mp_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_4 * lp_91[k]
                   + pb_x[k] * mp_91[k];

        t_182[k] = f_4 * lp_92[k]
                   + pb_x[k] * mp_92[k];

        t_183[k] = f_5 * kd0_183[k]
                   - f_6 * kd1_183[k]
                   + pa_x[k] * ld_183[k];

        t_184[k] = f_18 * lp_71[k]
                   + pb_y[k] * mp_92[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pa_x, pa_y, pb_x, kd0_108, kd0_185, kd1_108, \
                         kd1_185, lp_94, ld_144, ld_185, mp_94 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = f_5 * kd0_185[k]
                   - f_6 * kd1_185[k]
                   + pa_x[k] * ld_185[k];

        t_186[k] = f_16 * kd0_108[k]
                   - f_17 * kd1_108[k]
                   + pa_y[k] * ld_144[k];

        t_187[k] = f_4 * lp_94[k]
                   + pb_x[k] * mp_94[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pb_x, pb_y, kd0_189, kd0_191, \
                         kd1_189, kd1_191, lp_74, lp_95, ld_189, ld_191, \
                         mp_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_4 * lp_95[k]
                   + pb_x[k] * mp_95[k];

        t_189[k] = f_5 * kd0_189[k]
                   - f_6 * kd1_189[k]
                   + pa_x[k] * ld_189[k];

        t_190[k] = f_22 * lp_74[k]
                   + pb_y[k] * mp_95[k];

        t_191[k] = f_5 * kd0_191[k]
                   - f_6 * kd1_191[k]
                   + pa_x[k] * ld_191[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_y, pb_x, kd0_114, kd1_114, lp_97, lp_98, \
                         ld_150, mp_97, mp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_11 * kd0_114[k]
                   - f_12 * kd1_114[k]
                   + pa_y[k] * ld_150[k];

        t_193[k] = f_4 * lp_97[k]
                   + pb_x[k] * mp_97[k];

        t_194[k] = f_4 * lp_98[k]
                   + pb_x[k] * mp_98[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_x, pb_y, kd0_195, kd0_197, kd1_195, kd1_197, \
                         lp_77, ld_195, ld_197, mp_98 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_5 * kd0_195[k]
                   - f_6 * kd1_195[k]
                   + pa_x[k] * ld_195[k];

        t_196[k] = f_21 * lp_77[k]
                   + pb_y[k] * mp_98[k];

        t_197[k] = f_5 * kd0_197[k]
                   - f_6 * kd1_197[k]
                   + pa_x[k] * ld_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_y, pb_x, kd0_120, kd1_120, lp_100, lp_101, \
                         ld_156, mp_100, mp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_5 * kd0_120[k]
                   - f_6 * kd1_120[k]
                   + pa_y[k] * ld_156[k];

        t_199[k] = f_4 * lp_100[k]
                   + pb_x[k] * mp_100[k];

        t_200[k] = f_4 * lp_101[k]
                   + pb_x[k] * mp_101[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, t_204, pa_x, pa_y, pb_y, kd0_201, kd0_203, \
                         kd1_201, kd1_203, lp_80, ld_162, ld_201, ld_203, \
                         mp_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_5 * kd0_201[k]
                   - f_6 * kd1_201[k]
                   + pa_x[k] * ld_201[k];

        t_202[k] = f_4 * lp_80[k]
                   + pb_y[k] * mp_101[k];

        t_203[k] = f_5 * kd0_203[k]
                   - f_6 * kd1_203[k]
                   + pa_x[k] * ld_203[k];

        t_204[k] = pa_y[k] * ld_162[k];
    }

#pragma omp simd aligned(t_205, t_206, t_207, t_208, t_209, pa_y, pb_x, pb_y, lp_82, lp_83, \
                         lp_103, ld_164, ld_165, ld_167, mp_103, \
                         mp_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = f_4 * lp_103[k]
                   + pb_x[k] * mp_103[k];

        t_206[k] = pa_y[k] * ld_164[k];

        t_207[k] = f_4 * lp_82[k]
                   + pa_y[k] * ld_165[k];

        t_208[k] = f_10 * lp_83[k]
                   + pb_y[k] * mp_104[k];

        t_209[k] = pa_y[k] * ld_167[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, t_213, pa_z, pb_x, pb_y, kd0_120, kd1_120, \
                         lp_107, ld_162, ms0_35, ms1_35, mp_105, mp_106, \
                         mp_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = f_8 * kd0_120[k]
                   - f_9 * kd1_120[k]
                   + pa_z[k] * ld_162[k];

        t_211[k] = pb_y[k] * mp_105[k];

        t_212[k] = f_4 * lp_107[k]
                   + pb_x[k] * mp_107[k];

        t_213[k] = f_1 * ms0_35[k]
                   - f_2 * ms1_35[k]
                   + pb_y[k] * mp_106[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, t_217, pa_x, pb_x, pb_y, kd0_215, kd1_215, \
                         lp_108, lp_109, ld_215, ld_216, mp_107, \
                         mp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = pb_y[k] * mp_107[k];

        t_215[k] = f_5 * kd0_215[k]
                   - f_6 * kd1_215[k]
                   + pa_x[k] * ld_215[k];

        t_216[k] = f_4 * lp_108[k]
                   + pa_x[k] * ld_216[k];

        t_217[k] = f_10 * lp_109[k]
                   + pb_x[k] * mp_109[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, t_222, t_223, pa_x, pa_z, pb_z, ld_168, \
                         ld_169, ld_219, ld_221, mp_108, mp_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = pb_z[k] * mp_108[k];

        t_219[k] = pa_x[k] * ld_219[k];

        t_220[k] = pb_z[k] * mp_109[k];

        t_221[k] = pa_x[k] * ld_221[k];

        t_222[k] = pa_z[k] * ld_168[k];

        t_223[k] = pa_z[k] * ld_169[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, t_228, pa_x, pb_x, lp_113, lp_114, \
                         ld_225, ld_226, ld_227, ld_228, mp_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_10 * lp_113[k]
                   + pb_x[k] * mp_113[k];

        t_225[k] = pa_x[k] * ld_225[k];

        t_226[k] = pa_x[k] * ld_226[k];

        t_227[k] = pa_x[k] * ld_227[k];

        t_228[k] = f_4 * lp_114[k]
                   + pa_x[k] * ld_228[k];
    }

#pragma omp simd aligned(t_229, t_230, t_231, t_232, t_233, pa_x, pb_x, lp_115, lp_116, \
                         ld_231, ld_232, ld_233, mp_115, mp_116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_229[k] = f_10 * lp_115[k]
                   + pb_x[k] * mp_115[k];

        t_230[k] = f_10 * lp_116[k]
                   + pb_x[k] * mp_116[k];

        t_231[k] = pa_x[k] * ld_231[k];

        t_232[k] = pa_x[k] * ld_232[k];

        t_233[k] = pa_x[k] * ld_233[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, t_238, pa_x, pb_x, lp_117, lp_118, \
                         lp_119, ld_234, ld_237, ld_238, mp_118, \
                         mp_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_4 * lp_117[k]
                   + pa_x[k] * ld_234[k];

        t_235[k] = f_10 * lp_118[k]
                   + pb_x[k] * mp_118[k];

        t_236[k] = f_10 * lp_119[k]
                   + pb_x[k] * mp_119[k];

        t_237[k] = pa_x[k] * ld_237[k];

        t_238[k] = pa_x[k] * ld_238[k];
    }

#pragma omp simd aligned(t_239, t_240, t_241, t_242, t_243, pa_x, pb_x, lp_120, lp_121, \
                         lp_122, ld_239, ld_240, ld_243, mp_121, \
                         mp_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_239[k] = pa_x[k] * ld_239[k];

        t_240[k] = f_4 * lp_120[k]
                   + pa_x[k] * ld_240[k];

        t_241[k] = f_10 * lp_121[k]
                   + pb_x[k] * mp_121[k];

        t_242[k] = f_10 * lp_122[k]
                   + pb_x[k] * mp_122[k];

        t_243[k] = pa_x[k] * ld_243[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, t_248, pa_x, pb_x, lp_123, lp_124, \
                         lp_125, ld_244, ld_245, ld_246, mp_124, \
                         mp_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = pa_x[k] * ld_244[k];

        t_245[k] = pa_x[k] * ld_245[k];

        t_246[k] = f_4 * lp_123[k]
                   + pa_x[k] * ld_246[k];

        t_247[k] = f_10 * lp_124[k]
                   + pb_x[k] * mp_124[k];

        t_248[k] = f_10 * lp_125[k]
                   + pb_x[k] * mp_125[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, t_252, t_253, pa_x, pb_x, lp_126, lp_127, \
                         ld_249, ld_250, ld_251, ld_252, mp_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = pa_x[k] * ld_249[k];

        t_250[k] = pa_x[k] * ld_250[k];

        t_251[k] = pa_x[k] * ld_251[k];

        t_252[k] = f_4 * lp_126[k]
                   + pa_x[k] * ld_252[k];

        t_253[k] = f_10 * lp_127[k]
                   + pb_x[k] * mp_127[k];
    }

#pragma omp simd aligned(t_254, t_255, t_256, t_257, t_258, pa_x, pa_y, pb_x, lp_128, ld_210, \
                         ld_255, ld_256, ld_257, mp_128 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = f_10 * lp_128[k]
                   + pb_x[k] * mp_128[k];

        t_255[k] = pa_x[k] * ld_255[k];

        t_256[k] = pa_x[k] * ld_256[k];

        t_257[k] = pa_x[k] * ld_257[k];

        t_258[k] = pa_y[k] * ld_210[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, t_262, t_263, pa_x, pa_y, pb_x, lp_130, ld_212, \
                         ld_261, ld_262, ld_263, mp_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_10 * lp_130[k]
                   + pb_x[k] * mp_130[k];

        t_260[k] = pa_y[k] * ld_212[k];

        t_261[k] = pa_x[k] * ld_261[k];

        t_262[k] = pa_x[k] * ld_262[k];

        t_263[k] = pa_x[k] * ld_263[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, t_267, t_268, t_269, pa_x, pb_x, pb_y, lp_132, \
                         lp_134, ld_264, ld_267, ld_269, mp_132, \
                         mp_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_4 * lp_132[k]
                   + pa_x[k] * ld_264[k];

        t_265[k] = pb_y[k] * mp_132[k];

        t_266[k] = f_10 * lp_134[k]
                   + pb_x[k] * mp_134[k];

        t_267[k] = pa_x[k] * ld_267[k];

        t_268[k] = pb_y[k] * mp_134[k];

        t_269[k] = pa_x[k] * ld_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, t_273, t_274, t_275, pb_x, pb_y, pb_z, lp_109, \
                         ms0_45, ms1_45, mp_135, mp_136, mp_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_1 * ms0_45[k]
                   - f_2 * ms1_45[k]
                   + pb_x[k] * mp_135[k];

        t_271[k] = pb_x[k] * mp_136[k];

        t_272[k] = pb_x[k] * mp_137[k];

        t_273[k] = f_0 * lp_109[k]
                   + f_1 * ms0_45[k]
                   - f_2 * ms1_45[k]
                   + pb_y[k] * mp_136[k];

        t_274[k] = pb_z[k] * mp_136[k];

        t_275[k] = f_1 * ms0_45[k]
                   - f_2 * ms1_45[k]
                   + pb_z[k] * mp_137[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, t_280, t_281, pa_z, pb_x, pb_y, lp_110, \
                         lp_113, ld_216, ld_219, ld_221, mp_139, \
                         mp_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = pa_z[k] * ld_216[k];

        t_277[k] = pb_x[k] * mp_139[k];

        t_278[k] = pb_x[k] * mp_140[k];

        t_279[k] = pa_z[k] * ld_219[k];

        t_280[k] = f_3 * lp_113[k]
                   + pb_y[k] * mp_140[k];

        t_281[k] = f_4 * lp_110[k]
                   + pa_z[k] * ld_221[k];
    }
}

static auto
compute_prim_md_electron_repulsion_0_piece2(CSimdMatrix &buffer, const size_t target,
                                            const size_t pa, const size_t pb, const size_t kd0,
                                            const size_t kd1, const size_t lp, const size_t ld,
                                            const size_t ms0, const size_t ms1, const size_t mp,
                                            const size_t ncols, const double alpha,
                                            const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 4.5 / p;
    const auto f_1 = 0.5 / beta;
    const auto f_2 = 0.5 * alpha / (beta * p);
    const auto f_4 = 1.0 / p;
    const auto f_5 = 0.5 / alpha;
    const auto f_6 = 0.5 * beta / (alpha * p);
    const auto f_7 = 3.5 / p;
    const auto f_8 = 3.0 / alpha;
    const auto f_9 = 3.0 * beta / (alpha * p);
    const auto f_10 = 0.5 / p;
    const auto f_11 = 1.0 / alpha;
    const auto f_12 = beta / (alpha * p);
    const auto f_13 = 3.0 / p;
    const auto f_14 = 2.5 / alpha;
    const auto f_15 = 2.5 * beta / (alpha * p);
    const auto f_16 = 1.5 / alpha;
    const auto f_17 = 1.5 * beta / (alpha * p);
    const auto f_18 = 2.5 / p;
    const auto f_19 = 2.0 / alpha;
    const auto f_20 = 2.0 * beta / (alpha * p);
    const auto f_21 = 1.5 / p;
    const auto f_22 = 2.0 / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *kd0_171 = buffer.data(kd0 + 171);
    const auto *kd0_177 = buffer.data(kd0 + 177);
    const auto *kd0_183 = buffer.data(kd0 + 183);
    const auto *kd0_185 = buffer.data(kd0 + 185);
    const auto *kd0_189 = buffer.data(kd0 + 189);
    const auto *kd0_191 = buffer.data(kd0 + 191);
    const auto *kd0_195 = buffer.data(kd0 + 195);
    const auto *kd0_197 = buffer.data(kd0 + 197);
    const auto *kd0_201 = buffer.data(kd0 + 201);
    const auto *kd0_203 = buffer.data(kd0 + 203);
    const auto *kd0_209 = buffer.data(kd0 + 209);
    const auto *kd0_215 = buffer.data(kd0 + 215);

    const auto *kd1_171 = buffer.data(kd1 + 171);
    const auto *kd1_177 = buffer.data(kd1 + 177);
    const auto *kd1_183 = buffer.data(kd1 + 183);
    const auto *kd1_185 = buffer.data(kd1 + 185);
    const auto *kd1_189 = buffer.data(kd1 + 189);
    const auto *kd1_191 = buffer.data(kd1 + 191);
    const auto *kd1_195 = buffer.data(kd1 + 195);
    const auto *kd1_197 = buffer.data(kd1 + 197);
    const auto *kd1_201 = buffer.data(kd1 + 201);
    const auto *kd1_203 = buffer.data(kd1 + 203);
    const auto *kd1_209 = buffer.data(kd1 + 209);
    const auto *kd1_215 = buffer.data(kd1 + 215);

    const auto *lp_116 = buffer.data(lp + 116);
    const auto *lp_119 = buffer.data(lp + 119);
    const auto *lp_122 = buffer.data(lp + 122);
    const auto *lp_125 = buffer.data(lp + 125);
    const auto *lp_128 = buffer.data(lp + 128);
    const auto *lp_131 = buffer.data(lp + 131);
    const auto *lp_133 = buffer.data(lp + 133);
    const auto *lp_134 = buffer.data(lp + 134);

    const auto *ld_225 = buffer.data(ld + 225);
    const auto *ld_231 = buffer.data(ld + 231);
    const auto *ld_233 = buffer.data(ld + 233);
    const auto *ld_237 = buffer.data(ld + 237);
    const auto *ld_239 = buffer.data(ld + 239);
    const auto *ld_243 = buffer.data(ld + 243);
    const auto *ld_245 = buffer.data(ld + 245);
    const auto *ld_249 = buffer.data(ld + 249);
    const auto *ld_251 = buffer.data(ld + 251);
    const auto *ld_255 = buffer.data(ld + 255);
    const auto *ld_257 = buffer.data(ld + 257);
    const auto *ld_263 = buffer.data(ld + 263);
    const auto *ld_264 = buffer.data(ld + 264);
    const auto *ld_267 = buffer.data(ld + 267);
    const auto *ld_269 = buffer.data(ld + 269);

    const auto *ms0_47 = buffer.data(ms0 + 47);
    const auto *ms0_48 = buffer.data(ms0 + 48);
    const auto *ms0_49 = buffer.data(ms0 + 49);
    const auto *ms0_50 = buffer.data(ms0 + 50);
    const auto *ms0_51 = buffer.data(ms0 + 51);
    const auto *ms0_52 = buffer.data(ms0 + 52);
    const auto *ms0_54 = buffer.data(ms0 + 54);

    const auto *ms1_47 = buffer.data(ms1 + 47);
    const auto *ms1_48 = buffer.data(ms1 + 48);
    const auto *ms1_49 = buffer.data(ms1 + 49);
    const auto *ms1_50 = buffer.data(ms1 + 50);
    const auto *ms1_51 = buffer.data(ms1 + 51);
    const auto *ms1_52 = buffer.data(ms1 + 52);
    const auto *ms1_54 = buffer.data(ms1 + 54);

    const auto *mp_141 = buffer.data(mp + 141);
    const auto *mp_142 = buffer.data(mp + 142);
    const auto *mp_143 = buffer.data(mp + 143);
    const auto *mp_144 = buffer.data(mp + 144);
    const auto *mp_145 = buffer.data(mp + 145);
    const auto *mp_146 = buffer.data(mp + 146);
    const auto *mp_147 = buffer.data(mp + 147);
    const auto *mp_148 = buffer.data(mp + 148);
    const auto *mp_149 = buffer.data(mp + 149);
    const auto *mp_150 = buffer.data(mp + 150);
    const auto *mp_151 = buffer.data(mp + 151);
    const auto *mp_152 = buffer.data(mp + 152);
    const auto *mp_153 = buffer.data(mp + 153);
    const auto *mp_154 = buffer.data(mp + 154);
    const auto *mp_155 = buffer.data(mp + 155);
    const auto *mp_156 = buffer.data(mp + 156);
    const auto *mp_157 = buffer.data(mp + 157);
    const auto *mp_158 = buffer.data(mp + 158);
    const auto *mp_160 = buffer.data(mp + 160);
    const auto *mp_161 = buffer.data(mp + 161);
    const auto *mp_162 = buffer.data(mp + 162);
    const auto *mp_163 = buffer.data(mp + 163);
    const auto *mp_164 = buffer.data(mp + 164);

#pragma omp simd aligned(t_282, t_283, t_284, t_285, pa_z, pb_x, kd0_171, kd1_171, ld_225, \
                         ms0_47, ms1_47, mp_141, mp_142, mp_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = f_1 * ms0_47[k]
                   - f_2 * ms1_47[k]
                   + pb_x[k] * mp_141[k];

        t_283[k] = pb_x[k] * mp_142[k];

        t_284[k] = pb_x[k] * mp_143[k];

        t_285[k] = f_5 * kd0_171[k]
                   - f_6 * kd1_171[k]
                   + pa_z[k] * ld_225[k];
    }

#pragma omp simd aligned(t_286, t_287, t_288, t_289, pa_y, pb_x, pb_y, kd0_185, kd1_185, \
                         lp_116, ld_233, ms0_48, ms1_48, mp_143, mp_144, \
                         mp_145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_286[k] = f_7 * lp_116[k]
                   + pb_y[k] * mp_143[k];

        t_287[k] = f_8 * kd0_185[k]
                   - f_9 * kd1_185[k]
                   + pa_y[k] * ld_233[k];

        t_288[k] = f_1 * ms0_48[k]
                   - f_2 * ms1_48[k]
                   + pb_x[k] * mp_144[k];

        t_289[k] = pb_x[k] * mp_145[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, t_293, pa_y, pa_z, pb_x, pb_y, kd0_177, kd0_191, \
                         kd1_177, kd1_191, lp_119, ld_231, ld_239, \
                         mp_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pb_x[k] * mp_146[k];

        t_291[k] = f_11 * kd0_177[k]
                   - f_12 * kd1_177[k]
                   + pa_z[k] * ld_231[k];

        t_292[k] = f_13 * lp_119[k]
                   + pb_y[k] * mp_146[k];

        t_293[k] = f_14 * kd0_191[k]
                   - f_15 * kd1_191[k]
                   + pa_y[k] * ld_239[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, t_297, pa_z, pb_x, kd0_183, kd1_183, ld_237, \
                         ms0_49, ms1_49, mp_147, mp_148, mp_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = f_1 * ms0_49[k]
                   - f_2 * ms1_49[k]
                   + pb_x[k] * mp_147[k];

        t_295[k] = pb_x[k] * mp_148[k];

        t_296[k] = pb_x[k] * mp_149[k];

        t_297[k] = f_16 * kd0_183[k]
                   - f_17 * kd1_183[k]
                   + pa_z[k] * ld_237[k];
    }

#pragma omp simd aligned(t_298, t_299, t_300, t_301, pa_y, pb_x, pb_y, kd0_197, kd1_197, \
                         lp_122, ld_245, ms0_50, ms1_50, mp_149, mp_150, \
                         mp_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_298[k] = f_18 * lp_122[k]
                   + pb_y[k] * mp_149[k];

        t_299[k] = f_19 * kd0_197[k]
                   - f_20 * kd1_197[k]
                   + pa_y[k] * ld_245[k];

        t_300[k] = f_1 * ms0_50[k]
                   - f_2 * ms1_50[k]
                   + pb_x[k] * mp_150[k];

        t_301[k] = pb_x[k] * mp_151[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, t_305, pa_y, pa_z, pb_x, pb_y, kd0_189, kd0_203, \
                         kd1_189, kd1_203, lp_125, ld_243, ld_251, \
                         mp_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = pb_x[k] * mp_152[k];

        t_303[k] = f_19 * kd0_189[k]
                   - f_20 * kd1_189[k]
                   + pa_z[k] * ld_243[k];

        t_304[k] = f_22 * lp_125[k]
                   + pb_y[k] * mp_152[k];

        t_305[k] = f_16 * kd0_203[k]
                   - f_17 * kd1_203[k]
                   + pa_y[k] * ld_251[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, t_309, pa_z, pb_x, kd0_195, kd1_195, ld_249, \
                         ms0_51, ms1_51, mp_153, mp_154, mp_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_1 * ms0_51[k]
                   - f_2 * ms1_51[k]
                   + pb_x[k] * mp_153[k];

        t_307[k] = pb_x[k] * mp_154[k];

        t_308[k] = pb_x[k] * mp_155[k];

        t_309[k] = f_14 * kd0_195[k]
                   - f_15 * kd1_195[k]
                   + pa_z[k] * ld_249[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, t_313, pa_y, pb_x, pb_y, kd0_209, kd1_209, \
                         lp_128, ld_257, ms0_52, ms1_52, mp_155, mp_156, \
                         mp_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = f_21 * lp_128[k]
                   + pb_y[k] * mp_155[k];

        t_311[k] = f_11 * kd0_209[k]
                   - f_12 * kd1_209[k]
                   + pa_y[k] * ld_257[k];

        t_312[k] = f_1 * ms0_52[k]
                   - f_2 * ms1_52[k]
                   + pb_x[k] * mp_156[k];

        t_313[k] = pb_x[k] * mp_157[k];
    }

#pragma omp simd aligned(t_314, t_315, t_316, t_317, pa_y, pa_z, pb_x, pb_y, kd0_201, kd0_215, \
                         kd1_201, kd1_215, lp_131, ld_255, ld_263, \
                         mp_158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_314[k] = pb_x[k] * mp_158[k];

        t_315[k] = f_8 * kd0_201[k]
                   - f_9 * kd1_201[k]
                   + pa_z[k] * ld_255[k];

        t_316[k] = f_4 * lp_131[k]
                   + pb_y[k] * mp_158[k];

        t_317[k] = f_5 * kd0_215[k]
                   - f_6 * kd1_215[k]
                   + pa_y[k] * ld_263[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, t_321, t_322, t_323, pa_y, pb_x, pb_y, lp_133, \
                         lp_134, ld_264, ld_267, ld_269, mp_160, \
                         mp_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = pa_y[k] * ld_264[k];

        t_319[k] = pb_x[k] * mp_160[k];

        t_320[k] = pb_x[k] * mp_161[k];

        t_321[k] = f_4 * lp_133[k]
                   + pa_y[k] * ld_267[k];

        t_322[k] = f_10 * lp_134[k]
                   + pb_y[k] * mp_161[k];

        t_323[k] = pa_y[k] * ld_269[k];
    }

#pragma omp simd aligned(t_324, t_325, t_326, t_327, t_328, t_329, pb_x, pb_y, pb_z, lp_134, \
                         ms0_54, ms1_54, mp_162, mp_163, mp_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_324[k] = f_1 * ms0_54[k]
                   - f_2 * ms1_54[k]
                   + pb_x[k] * mp_162[k];

        t_325[k] = pb_x[k] * mp_163[k];

        t_326[k] = pb_x[k] * mp_164[k];

        t_327[k] = f_1 * ms0_54[k]
                   - f_2 * ms1_54[k]
                   + pb_y[k] * mp_163[k];

        t_328[k] = pb_y[k] * mp_164[k];

        t_329[k] = f_0 * lp_134[k]
                   + f_1 * ms0_54[k]
                   - f_2 * ms1_54[k]
                   + pb_z[k] * mp_164[k];
    }
}

auto
compute_prim_md_electron_repulsion_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                     const size_t pb, const size_t kd0, const size_t kd1,
                                     const size_t lp, const size_t ld, const size_t ms0,
                                     const size_t ms1, const size_t mp, const size_t ncols,
                                     const double alpha, const double beta,
                                     const double p) -> void
{
    compute_prim_md_electron_repulsion_0_piece0(buffer, target, pa, pb, kd0, kd1, lp, ld, ms0,
                                                ms1, mp, ncols, alpha, beta, p);

    compute_prim_md_electron_repulsion_0_piece1(buffer, target, pa, pb, kd0, kd1, lp, ld, ms0,
                                                ms1, mp, ncols, alpha, beta, p);

    compute_prim_md_electron_repulsion_0_piece2(buffer, target, pa, pb, kd0, kd1, lp, ld, ms0,
                                                ms1, mp, ncols, alpha, beta, p);
}

}  // namespace simdt2ceri
