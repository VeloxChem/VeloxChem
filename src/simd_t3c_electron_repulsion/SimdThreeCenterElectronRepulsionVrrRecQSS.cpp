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


#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"

#include "SimdAlign.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

static auto
compute_prim_qss_three_center_electron_repulsion_0_piece0(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nss0,
                                                          const size_t nss1, const size_t oss0,
                                                          const size_t oss1, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / p;
    const auto f_1 = 5.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_3 = 4.5 / p;
    const auto f_4 = 4.5 * gamma / (p * q);
    const auto f_5 = 4.0 / p;
    const auto f_6 = 4.0 * gamma / (p * q);
    const auto f_7 = 3.5 / p;
    const auto f_8 = 3.5 * gamma / (p * q);
    const auto f_9 = 3.0 / p;
    const auto f_10 = 3.0 * gamma / (p * q);
    const auto f_11 = 2.5 / p;
    const auto f_12 = 2.5 * gamma / (p * q);
    const auto f_13 = 2.0 / p;
    const auto f_14 = 2.0 * gamma / (p * q);
    const auto f_15 = 1.5 / p;
    const auto f_16 = 1.5 * gamma / (p * q);
    const auto f_17 = 1.0 / p;
    const auto f_18 = gamma / (p * q);
    const auto f_19 = 0.5 / p;
    const auto f_20 = 0.5 * gamma / (p * q);

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_x = buffer.data(pc + 0);
    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nss0_0 = buffer.data(nss0 + 0);
    const auto *nss0_3 = buffer.data(nss0 + 3);
    const auto *nss0_5 = buffer.data(nss0 + 5);
    const auto *nss0_6 = buffer.data(nss0 + 6);
    const auto *nss0_9 = buffer.data(nss0 + 9);
    const auto *nss0_10 = buffer.data(nss0 + 10);
    const auto *nss0_12 = buffer.data(nss0 + 12);
    const auto *nss0_14 = buffer.data(nss0 + 14);
    const auto *nss0_15 = buffer.data(nss0 + 15);
    const auto *nss0_17 = buffer.data(nss0 + 17);
    const auto *nss0_18 = buffer.data(nss0 + 18);
    const auto *nss0_20 = buffer.data(nss0 + 20);
    const auto *nss0_21 = buffer.data(nss0 + 21);
    const auto *nss0_23 = buffer.data(nss0 + 23);
    const auto *nss0_24 = buffer.data(nss0 + 24);
    const auto *nss0_25 = buffer.data(nss0 + 25);
    const auto *nss0_27 = buffer.data(nss0 + 27);
    const auto *nss0_28 = buffer.data(nss0 + 28);
    const auto *nss0_30 = buffer.data(nss0 + 30);
    const auto *nss0_31 = buffer.data(nss0 + 31);
    const auto *nss0_32 = buffer.data(nss0 + 32);
    const auto *nss0_33 = buffer.data(nss0 + 33);
    const auto *nss0_35 = buffer.data(nss0 + 35);
    const auto *nss0_36 = buffer.data(nss0 + 36);
    const auto *nss0_38 = buffer.data(nss0 + 38);
    const auto *nss0_39 = buffer.data(nss0 + 39);
    const auto *nss0_40 = buffer.data(nss0 + 40);
    const auto *nss0_41 = buffer.data(nss0 + 41);
    const auto *nss0_42 = buffer.data(nss0 + 42);
    const auto *nss0_44 = buffer.data(nss0 + 44);
    const auto *nss0_45 = buffer.data(nss0 + 45);
    const auto *nss0_47 = buffer.data(nss0 + 47);
    const auto *nss0_48 = buffer.data(nss0 + 48);
    const auto *nss0_49 = buffer.data(nss0 + 49);
    const auto *nss0_50 = buffer.data(nss0 + 50);
    const auto *nss0_51 = buffer.data(nss0 + 51);
    const auto *nss0_52 = buffer.data(nss0 + 52);
    const auto *nss0_54 = buffer.data(nss0 + 54);
    const auto *nss0_55 = buffer.data(nss0 + 55);
    const auto *nss0_57 = buffer.data(nss0 + 57);
    const auto *nss0_58 = buffer.data(nss0 + 58);
    const auto *nss0_59 = buffer.data(nss0 + 59);
    const auto *nss0_60 = buffer.data(nss0 + 60);
    const auto *nss0_61 = buffer.data(nss0 + 61);
    const auto *nss0_62 = buffer.data(nss0 + 62);
    const auto *nss0_63 = buffer.data(nss0 + 63);
    const auto *nss0_64 = buffer.data(nss0 + 64);
    const auto *nss0_65 = buffer.data(nss0 + 65);

    const auto *nss1_0 = buffer.data(nss1 + 0);
    const auto *nss1_3 = buffer.data(nss1 + 3);
    const auto *nss1_5 = buffer.data(nss1 + 5);
    const auto *nss1_6 = buffer.data(nss1 + 6);
    const auto *nss1_9 = buffer.data(nss1 + 9);
    const auto *nss1_10 = buffer.data(nss1 + 10);
    const auto *nss1_12 = buffer.data(nss1 + 12);
    const auto *nss1_14 = buffer.data(nss1 + 14);
    const auto *nss1_15 = buffer.data(nss1 + 15);
    const auto *nss1_17 = buffer.data(nss1 + 17);
    const auto *nss1_18 = buffer.data(nss1 + 18);
    const auto *nss1_20 = buffer.data(nss1 + 20);
    const auto *nss1_21 = buffer.data(nss1 + 21);
    const auto *nss1_23 = buffer.data(nss1 + 23);
    const auto *nss1_24 = buffer.data(nss1 + 24);
    const auto *nss1_25 = buffer.data(nss1 + 25);
    const auto *nss1_27 = buffer.data(nss1 + 27);
    const auto *nss1_28 = buffer.data(nss1 + 28);
    const auto *nss1_30 = buffer.data(nss1 + 30);
    const auto *nss1_31 = buffer.data(nss1 + 31);
    const auto *nss1_32 = buffer.data(nss1 + 32);
    const auto *nss1_33 = buffer.data(nss1 + 33);
    const auto *nss1_35 = buffer.data(nss1 + 35);
    const auto *nss1_36 = buffer.data(nss1 + 36);
    const auto *nss1_38 = buffer.data(nss1 + 38);
    const auto *nss1_39 = buffer.data(nss1 + 39);
    const auto *nss1_40 = buffer.data(nss1 + 40);
    const auto *nss1_41 = buffer.data(nss1 + 41);
    const auto *nss1_42 = buffer.data(nss1 + 42);
    const auto *nss1_44 = buffer.data(nss1 + 44);
    const auto *nss1_45 = buffer.data(nss1 + 45);
    const auto *nss1_47 = buffer.data(nss1 + 47);
    const auto *nss1_48 = buffer.data(nss1 + 48);
    const auto *nss1_49 = buffer.data(nss1 + 49);
    const auto *nss1_50 = buffer.data(nss1 + 50);
    const auto *nss1_51 = buffer.data(nss1 + 51);
    const auto *nss1_52 = buffer.data(nss1 + 52);
    const auto *nss1_54 = buffer.data(nss1 + 54);
    const auto *nss1_55 = buffer.data(nss1 + 55);
    const auto *nss1_57 = buffer.data(nss1 + 57);
    const auto *nss1_58 = buffer.data(nss1 + 58);
    const auto *nss1_59 = buffer.data(nss1 + 59);
    const auto *nss1_60 = buffer.data(nss1 + 60);
    const auto *nss1_61 = buffer.data(nss1 + 61);
    const auto *nss1_62 = buffer.data(nss1 + 62);
    const auto *nss1_63 = buffer.data(nss1 + 63);
    const auto *nss1_64 = buffer.data(nss1 + 64);
    const auto *nss1_65 = buffer.data(nss1 + 65);

    const auto *oss0_0 = buffer.data(oss0 + 0);
    const auto *oss0_2 = buffer.data(oss0 + 2);
    const auto *oss0_3 = buffer.data(oss0 + 3);
    const auto *oss0_5 = buffer.data(oss0 + 5);
    const auto *oss0_6 = buffer.data(oss0 + 6);
    const auto *oss0_9 = buffer.data(oss0 + 9);
    const auto *oss0_10 = buffer.data(oss0 + 10);
    const auto *oss0_12 = buffer.data(oss0 + 12);
    const auto *oss0_14 = buffer.data(oss0 + 14);
    const auto *oss0_15 = buffer.data(oss0 + 15);
    const auto *oss0_17 = buffer.data(oss0 + 17);
    const auto *oss0_18 = buffer.data(oss0 + 18);
    const auto *oss0_20 = buffer.data(oss0 + 20);
    const auto *oss0_21 = buffer.data(oss0 + 21);
    const auto *oss0_23 = buffer.data(oss0 + 23);
    const auto *oss0_24 = buffer.data(oss0 + 24);
    const auto *oss0_25 = buffer.data(oss0 + 25);
    const auto *oss0_27 = buffer.data(oss0 + 27);
    const auto *oss0_28 = buffer.data(oss0 + 28);
    const auto *oss0_30 = buffer.data(oss0 + 30);
    const auto *oss0_31 = buffer.data(oss0 + 31);
    const auto *oss0_32 = buffer.data(oss0 + 32);
    const auto *oss0_33 = buffer.data(oss0 + 33);
    const auto *oss0_35 = buffer.data(oss0 + 35);
    const auto *oss0_36 = buffer.data(oss0 + 36);
    const auto *oss0_38 = buffer.data(oss0 + 38);
    const auto *oss0_39 = buffer.data(oss0 + 39);
    const auto *oss0_40 = buffer.data(oss0 + 40);
    const auto *oss0_41 = buffer.data(oss0 + 41);
    const auto *oss0_42 = buffer.data(oss0 + 42);
    const auto *oss0_44 = buffer.data(oss0 + 44);
    const auto *oss0_45 = buffer.data(oss0 + 45);
    const auto *oss0_47 = buffer.data(oss0 + 47);
    const auto *oss0_48 = buffer.data(oss0 + 48);
    const auto *oss0_49 = buffer.data(oss0 + 49);
    const auto *oss0_50 = buffer.data(oss0 + 50);
    const auto *oss0_51 = buffer.data(oss0 + 51);
    const auto *oss0_52 = buffer.data(oss0 + 52);
    const auto *oss0_54 = buffer.data(oss0 + 54);
    const auto *oss0_55 = buffer.data(oss0 + 55);
    const auto *oss0_57 = buffer.data(oss0 + 57);
    const auto *oss0_58 = buffer.data(oss0 + 58);
    const auto *oss0_59 = buffer.data(oss0 + 59);
    const auto *oss0_60 = buffer.data(oss0 + 60);
    const auto *oss0_61 = buffer.data(oss0 + 61);
    const auto *oss0_62 = buffer.data(oss0 + 62);
    const auto *oss0_63 = buffer.data(oss0 + 63);
    const auto *oss0_65 = buffer.data(oss0 + 65);
    const auto *oss0_66 = buffer.data(oss0 + 66);
    const auto *oss0_67 = buffer.data(oss0 + 67);
    const auto *oss0_68 = buffer.data(oss0 + 68);
    const auto *oss0_69 = buffer.data(oss0 + 69);
    const auto *oss0_70 = buffer.data(oss0 + 70);
    const auto *oss0_71 = buffer.data(oss0 + 71);
    const auto *oss0_72 = buffer.data(oss0 + 72);
    const auto *oss0_73 = buffer.data(oss0 + 73);
    const auto *oss0_74 = buffer.data(oss0 + 74);
    const auto *oss0_75 = buffer.data(oss0 + 75);
    const auto *oss0_76 = buffer.data(oss0 + 76);
    const auto *oss0_77 = buffer.data(oss0 + 77);

    const auto *oss1_0 = buffer.data(oss1 + 0);
    const auto *oss1_2 = buffer.data(oss1 + 2);
    const auto *oss1_3 = buffer.data(oss1 + 3);
    const auto *oss1_5 = buffer.data(oss1 + 5);
    const auto *oss1_6 = buffer.data(oss1 + 6);
    const auto *oss1_9 = buffer.data(oss1 + 9);
    const auto *oss1_10 = buffer.data(oss1 + 10);
    const auto *oss1_12 = buffer.data(oss1 + 12);
    const auto *oss1_14 = buffer.data(oss1 + 14);
    const auto *oss1_15 = buffer.data(oss1 + 15);
    const auto *oss1_17 = buffer.data(oss1 + 17);
    const auto *oss1_18 = buffer.data(oss1 + 18);
    const auto *oss1_20 = buffer.data(oss1 + 20);
    const auto *oss1_21 = buffer.data(oss1 + 21);
    const auto *oss1_23 = buffer.data(oss1 + 23);
    const auto *oss1_24 = buffer.data(oss1 + 24);
    const auto *oss1_25 = buffer.data(oss1 + 25);
    const auto *oss1_27 = buffer.data(oss1 + 27);
    const auto *oss1_28 = buffer.data(oss1 + 28);
    const auto *oss1_30 = buffer.data(oss1 + 30);
    const auto *oss1_31 = buffer.data(oss1 + 31);
    const auto *oss1_32 = buffer.data(oss1 + 32);
    const auto *oss1_33 = buffer.data(oss1 + 33);
    const auto *oss1_35 = buffer.data(oss1 + 35);
    const auto *oss1_36 = buffer.data(oss1 + 36);
    const auto *oss1_38 = buffer.data(oss1 + 38);
    const auto *oss1_39 = buffer.data(oss1 + 39);
    const auto *oss1_40 = buffer.data(oss1 + 40);
    const auto *oss1_41 = buffer.data(oss1 + 41);
    const auto *oss1_42 = buffer.data(oss1 + 42);
    const auto *oss1_44 = buffer.data(oss1 + 44);
    const auto *oss1_45 = buffer.data(oss1 + 45);
    const auto *oss1_47 = buffer.data(oss1 + 47);
    const auto *oss1_48 = buffer.data(oss1 + 48);
    const auto *oss1_49 = buffer.data(oss1 + 49);
    const auto *oss1_50 = buffer.data(oss1 + 50);
    const auto *oss1_51 = buffer.data(oss1 + 51);
    const auto *oss1_52 = buffer.data(oss1 + 52);
    const auto *oss1_54 = buffer.data(oss1 + 54);
    const auto *oss1_55 = buffer.data(oss1 + 55);
    const auto *oss1_57 = buffer.data(oss1 + 57);
    const auto *oss1_58 = buffer.data(oss1 + 58);
    const auto *oss1_59 = buffer.data(oss1 + 59);
    const auto *oss1_60 = buffer.data(oss1 + 60);
    const auto *oss1_61 = buffer.data(oss1 + 61);
    const auto *oss1_62 = buffer.data(oss1 + 62);
    const auto *oss1_63 = buffer.data(oss1 + 63);
    const auto *oss1_65 = buffer.data(oss1 + 65);
    const auto *oss1_66 = buffer.data(oss1 + 66);
    const auto *oss1_67 = buffer.data(oss1 + 67);
    const auto *oss1_68 = buffer.data(oss1 + 68);
    const auto *oss1_69 = buffer.data(oss1 + 69);
    const auto *oss1_70 = buffer.data(oss1 + 70);
    const auto *oss1_71 = buffer.data(oss1 + 71);
    const auto *oss1_72 = buffer.data(oss1 + 72);
    const auto *oss1_73 = buffer.data(oss1 + 73);
    const auto *oss1_74 = buffer.data(oss1 + 74);
    const auto *oss1_75 = buffer.data(oss1 + 75);
    const auto *oss1_76 = buffer.data(oss1 + 76);
    const auto *oss1_77 = buffer.data(oss1 + 77);

#pragma omp simd aligned(t_0, t_1, t_2, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, nss0_0, nss1_0, \
                         oss0_0, oss1_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * nss0_0[k]
                 - f_1 * nss1_0[k]
                 + pa_x[k] * oss0_0[k]
                 - f_2 * pc_x[k] * oss1_0[k];

        t_1[k] = pa_y[k] * oss0_0[k]
                 - f_2 * pc_y[k] * oss1_0[k];

        t_2[k] = pa_z[k] * oss0_0[k]
                 - f_2 * pc_z[k] * oss1_0[k];
    }

#pragma omp simd aligned(t_3, t_4, pa_x, pa_y, pc_x, pc_y, nss0_3, nss1_3, oss0_2, oss0_3, \
                         oss1_2, oss1_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_3 * nss0_3[k]
                 - f_4 * nss1_3[k]
                 + pa_x[k] * oss0_3[k]
                 - f_2 * pc_x[k] * oss1_3[k];

        t_4[k] = pa_y[k] * oss0_2[k]
                 - f_2 * pc_y[k] * oss1_2[k];
    }

#pragma omp simd aligned(t_5, t_6, pa_x, pc_x, nss0_5, nss0_6, nss1_5, nss1_6, oss0_5, oss0_6, \
                         oss1_5, oss1_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_5[k] = f_3 * nss0_5[k]
                 - f_4 * nss1_5[k]
                 + pa_x[k] * oss0_5[k]
                 - f_2 * pc_x[k] * oss1_5[k];

        t_6[k] = f_5 * nss0_6[k]
                 - f_6 * nss1_6[k]
                 + pa_x[k] * oss0_6[k]
                 - f_2 * pc_x[k] * oss1_6[k];
    }

#pragma omp simd aligned(t_7, t_8, pa_y, pa_z, pc_y, pc_z, oss0_3, oss0_5, oss1_3, \
                         oss1_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = pa_z[k] * oss0_3[k]
                 - f_2 * pc_z[k] * oss1_3[k];

        t_8[k] = pa_y[k] * oss0_5[k]
                 - f_2 * pc_y[k] * oss1_5[k];
    }

#pragma omp simd aligned(t_9, t_10, pa_x, pc_x, nss0_9, nss0_10, nss1_9, nss1_10, oss0_9, \
                         oss0_10, oss1_9, oss1_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_9[k] = f_5 * nss0_9[k]
                 - f_6 * nss1_9[k]
                 + pa_x[k] * oss0_9[k]
                 - f_2 * pc_x[k] * oss1_9[k];

        t_10[k] = f_7 * nss0_10[k]
                  - f_8 * nss1_10[k]
                  + pa_x[k] * oss0_10[k]
                  - f_2 * pc_x[k] * oss1_10[k];
    }

#pragma omp simd aligned(t_11, t_12, pa_x, pa_z, pc_x, pc_z, nss0_12, nss1_12, oss0_6, \
                         oss0_12, oss1_6, oss1_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = pa_z[k] * oss0_6[k]
                  - f_2 * pc_z[k] * oss1_6[k];

        t_12[k] = f_7 * nss0_12[k]
                  - f_8 * nss1_12[k]
                  + pa_x[k] * oss0_12[k]
                  - f_2 * pc_x[k] * oss1_12[k];
    }

#pragma omp simd aligned(t_13, t_14, pa_x, pa_y, pc_x, pc_y, nss0_14, nss1_14, oss0_9, \
                         oss0_14, oss1_9, oss1_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = pa_y[k] * oss0_9[k]
                  - f_2 * pc_y[k] * oss1_9[k];

        t_14[k] = f_7 * nss0_14[k]
                  - f_8 * nss1_14[k]
                  + pa_x[k] * oss0_14[k]
                  - f_2 * pc_x[k] * oss1_14[k];
    }

#pragma omp simd aligned(t_15, t_16, pa_x, pa_z, pc_x, pc_z, nss0_15, nss1_15, oss0_10, \
                         oss0_15, oss1_10, oss1_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_15[k] = f_9 * nss0_15[k]
                  - f_10 * nss1_15[k]
                  + pa_x[k] * oss0_15[k]
                  - f_2 * pc_x[k] * oss1_15[k];

        t_16[k] = pa_z[k] * oss0_10[k]
                  - f_2 * pc_z[k] * oss1_10[k];
    }

#pragma omp simd aligned(t_17, t_18, pa_x, pc_x, nss0_17, nss0_18, nss1_17, nss1_18, oss0_17, \
                         oss0_18, oss1_17, oss1_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = f_9 * nss0_17[k]
                  - f_10 * nss1_17[k]
                  + pa_x[k] * oss0_17[k]
                  - f_2 * pc_x[k] * oss1_17[k];

        t_18[k] = f_9 * nss0_18[k]
                  - f_10 * nss1_18[k]
                  + pa_x[k] * oss0_18[k]
                  - f_2 * pc_x[k] * oss1_18[k];
    }

#pragma omp simd aligned(t_19, t_20, pa_x, pa_y, pc_x, pc_y, nss0_20, nss1_20, oss0_14, \
                         oss0_20, oss1_14, oss1_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * oss0_14[k]
                  - f_2 * pc_y[k] * oss1_14[k];

        t_20[k] = f_9 * nss0_20[k]
                  - f_10 * nss1_20[k]
                  + pa_x[k] * oss0_20[k]
                  - f_2 * pc_x[k] * oss1_20[k];
    }

#pragma omp simd aligned(t_21, t_22, pa_x, pa_z, pc_x, pc_z, nss0_21, nss1_21, oss0_15, \
                         oss0_21, oss1_15, oss1_21 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_21[k] = f_11 * nss0_21[k]
                  - f_12 * nss1_21[k]
                  + pa_x[k] * oss0_21[k]
                  - f_2 * pc_x[k] * oss1_21[k];

        t_22[k] = pa_z[k] * oss0_15[k]
                  - f_2 * pc_z[k] * oss1_15[k];
    }

#pragma omp simd aligned(t_23, t_24, pa_x, pc_x, nss0_23, nss0_24, nss1_23, nss1_24, oss0_23, \
                         oss0_24, oss1_23, oss1_24 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_11 * nss0_23[k]
                  - f_12 * nss1_23[k]
                  + pa_x[k] * oss0_23[k]
                  - f_2 * pc_x[k] * oss1_23[k];

        t_24[k] = f_11 * nss0_24[k]
                  - f_12 * nss1_24[k]
                  + pa_x[k] * oss0_24[k]
                  - f_2 * pc_x[k] * oss1_24[k];
    }

#pragma omp simd aligned(t_25, t_26, pa_x, pa_y, pc_x, pc_y, nss0_25, nss1_25, oss0_20, \
                         oss0_25, oss1_20, oss1_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_25[k] = f_11 * nss0_25[k]
                  - f_12 * nss1_25[k]
                  + pa_x[k] * oss0_25[k]
                  - f_2 * pc_x[k] * oss1_25[k];

        t_26[k] = pa_y[k] * oss0_20[k]
                  - f_2 * pc_y[k] * oss1_20[k];
    }

#pragma omp simd aligned(t_27, t_28, pa_x, pc_x, nss0_27, nss0_28, nss1_27, nss1_28, oss0_27, \
                         oss0_28, oss1_27, oss1_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_11 * nss0_27[k]
                  - f_12 * nss1_27[k]
                  + pa_x[k] * oss0_27[k]
                  - f_2 * pc_x[k] * oss1_27[k];

        t_28[k] = f_13 * nss0_28[k]
                  - f_14 * nss1_28[k]
                  + pa_x[k] * oss0_28[k]
                  - f_2 * pc_x[k] * oss1_28[k];
    }

#pragma omp simd aligned(t_29, t_30, pa_x, pa_z, pc_x, pc_z, nss0_30, nss1_30, oss0_21, \
                         oss0_30, oss1_21, oss1_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_z[k] * oss0_21[k]
                  - f_2 * pc_z[k] * oss1_21[k];

        t_30[k] = f_13 * nss0_30[k]
                  - f_14 * nss1_30[k]
                  + pa_x[k] * oss0_30[k]
                  - f_2 * pc_x[k] * oss1_30[k];
    }

#pragma omp simd aligned(t_31, t_32, pa_x, pc_x, nss0_31, nss0_32, nss1_31, nss1_32, oss0_31, \
                         oss0_32, oss1_31, oss1_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_31[k] = f_13 * nss0_31[k]
                  - f_14 * nss1_31[k]
                  + pa_x[k] * oss0_31[k]
                  - f_2 * pc_x[k] * oss1_31[k];

        t_32[k] = f_13 * nss0_32[k]
                  - f_14 * nss1_32[k]
                  + pa_x[k] * oss0_32[k]
                  - f_2 * pc_x[k] * oss1_32[k];
    }

#pragma omp simd aligned(t_33, t_34, pa_x, pa_y, pc_x, pc_y, nss0_33, nss1_33, oss0_27, \
                         oss0_33, oss1_27, oss1_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_13 * nss0_33[k]
                  - f_14 * nss1_33[k]
                  + pa_x[k] * oss0_33[k]
                  - f_2 * pc_x[k] * oss1_33[k];

        t_34[k] = pa_y[k] * oss0_27[k]
                  - f_2 * pc_y[k] * oss1_27[k];
    }

#pragma omp simd aligned(t_35, t_36, pa_x, pc_x, nss0_35, nss0_36, nss1_35, nss1_36, oss0_35, \
                         oss0_36, oss1_35, oss1_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_35[k] = f_13 * nss0_35[k]
                  - f_14 * nss1_35[k]
                  + pa_x[k] * oss0_35[k]
                  - f_2 * pc_x[k] * oss1_35[k];

        t_36[k] = f_15 * nss0_36[k]
                  - f_16 * nss1_36[k]
                  + pa_x[k] * oss0_36[k]
                  - f_2 * pc_x[k] * oss1_36[k];
    }

#pragma omp simd aligned(t_37, t_38, pa_x, pa_z, pc_x, pc_z, nss0_38, nss1_38, oss0_28, \
                         oss0_38, oss1_28, oss1_38 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = pa_z[k] * oss0_28[k]
                  - f_2 * pc_z[k] * oss1_28[k];

        t_38[k] = f_15 * nss0_38[k]
                  - f_16 * nss1_38[k]
                  + pa_x[k] * oss0_38[k]
                  - f_2 * pc_x[k] * oss1_38[k];
    }

#pragma omp simd aligned(t_39, t_40, pa_x, pc_x, nss0_39, nss0_40, nss1_39, nss1_40, oss0_39, \
                         oss0_40, oss1_39, oss1_40 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = f_15 * nss0_39[k]
                  - f_16 * nss1_39[k]
                  + pa_x[k] * oss0_39[k]
                  - f_2 * pc_x[k] * oss1_39[k];

        t_40[k] = f_15 * nss0_40[k]
                  - f_16 * nss1_40[k]
                  + pa_x[k] * oss0_40[k]
                  - f_2 * pc_x[k] * oss1_40[k];
    }

#pragma omp simd aligned(t_41, t_42, pa_x, pc_x, nss0_41, nss0_42, nss1_41, nss1_42, oss0_41, \
                         oss0_42, oss1_41, oss1_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_41[k] = f_15 * nss0_41[k]
                  - f_16 * nss1_41[k]
                  + pa_x[k] * oss0_41[k]
                  - f_2 * pc_x[k] * oss1_41[k];

        t_42[k] = f_15 * nss0_42[k]
                  - f_16 * nss1_42[k]
                  + pa_x[k] * oss0_42[k]
                  - f_2 * pc_x[k] * oss1_42[k];
    }

#pragma omp simd aligned(t_43, t_44, pa_x, pa_y, pc_x, pc_y, nss0_44, nss1_44, oss0_35, \
                         oss0_44, oss1_35, oss1_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = pa_y[k] * oss0_35[k]
                  - f_2 * pc_y[k] * oss1_35[k];

        t_44[k] = f_15 * nss0_44[k]
                  - f_16 * nss1_44[k]
                  + pa_x[k] * oss0_44[k]
                  - f_2 * pc_x[k] * oss1_44[k];
    }

#pragma omp simd aligned(t_45, t_46, pa_x, pa_z, pc_x, pc_z, nss0_45, nss1_45, oss0_36, \
                         oss0_45, oss1_36, oss1_45 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_45[k] = f_17 * nss0_45[k]
                  - f_18 * nss1_45[k]
                  + pa_x[k] * oss0_45[k]
                  - f_2 * pc_x[k] * oss1_45[k];

        t_46[k] = pa_z[k] * oss0_36[k]
                  - f_2 * pc_z[k] * oss1_36[k];
    }

#pragma omp simd aligned(t_47, t_48, pa_x, pc_x, nss0_47, nss0_48, nss1_47, nss1_48, oss0_47, \
                         oss0_48, oss1_47, oss1_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_47[k] = f_17 * nss0_47[k]
                  - f_18 * nss1_47[k]
                  + pa_x[k] * oss0_47[k]
                  - f_2 * pc_x[k] * oss1_47[k];

        t_48[k] = f_17 * nss0_48[k]
                  - f_18 * nss1_48[k]
                  + pa_x[k] * oss0_48[k]
                  - f_2 * pc_x[k] * oss1_48[k];
    }

#pragma omp simd aligned(t_49, t_50, pa_x, pc_x, nss0_49, nss0_50, nss1_49, nss1_50, oss0_49, \
                         oss0_50, oss1_49, oss1_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_17 * nss0_49[k]
                  - f_18 * nss1_49[k]
                  + pa_x[k] * oss0_49[k]
                  - f_2 * pc_x[k] * oss1_49[k];

        t_50[k] = f_17 * nss0_50[k]
                  - f_18 * nss1_50[k]
                  + pa_x[k] * oss0_50[k]
                  - f_2 * pc_x[k] * oss1_50[k];
    }

#pragma omp simd aligned(t_51, t_52, pa_x, pc_x, nss0_51, nss0_52, nss1_51, nss1_52, oss0_51, \
                         oss0_52, oss1_51, oss1_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_51[k] = f_17 * nss0_51[k]
                  - f_18 * nss1_51[k]
                  + pa_x[k] * oss0_51[k]
                  - f_2 * pc_x[k] * oss1_51[k];

        t_52[k] = f_17 * nss0_52[k]
                  - f_18 * nss1_52[k]
                  + pa_x[k] * oss0_52[k]
                  - f_2 * pc_x[k] * oss1_52[k];
    }

#pragma omp simd aligned(t_53, t_54, pa_x, pa_y, pc_x, pc_y, nss0_54, nss1_54, oss0_44, \
                         oss0_54, oss1_44, oss1_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = pa_y[k] * oss0_44[k]
                  - f_2 * pc_y[k] * oss1_44[k];

        t_54[k] = f_17 * nss0_54[k]
                  - f_18 * nss1_54[k]
                  + pa_x[k] * oss0_54[k]
                  - f_2 * pc_x[k] * oss1_54[k];
    }

#pragma omp simd aligned(t_55, t_56, pa_x, pa_z, pc_x, pc_z, nss0_55, nss1_55, oss0_45, \
                         oss0_55, oss1_45, oss1_55 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_19 * nss0_55[k]
                  - f_20 * nss1_55[k]
                  + pa_x[k] * oss0_55[k]
                  - f_2 * pc_x[k] * oss1_55[k];

        t_56[k] = pa_z[k] * oss0_45[k]
                  - f_2 * pc_z[k] * oss1_45[k];
    }

#pragma omp simd aligned(t_57, t_58, pa_x, pc_x, nss0_57, nss0_58, nss1_57, nss1_58, oss0_57, \
                         oss0_58, oss1_57, oss1_58 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_57[k] = f_19 * nss0_57[k]
                  - f_20 * nss1_57[k]
                  + pa_x[k] * oss0_57[k]
                  - f_2 * pc_x[k] * oss1_57[k];

        t_58[k] = f_19 * nss0_58[k]
                  - f_20 * nss1_58[k]
                  + pa_x[k] * oss0_58[k]
                  - f_2 * pc_x[k] * oss1_58[k];
    }

#pragma omp simd aligned(t_59, t_60, pa_x, pc_x, nss0_59, nss0_60, nss1_59, nss1_60, oss0_59, \
                         oss0_60, oss1_59, oss1_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = f_19 * nss0_59[k]
                  - f_20 * nss1_59[k]
                  + pa_x[k] * oss0_59[k]
                  - f_2 * pc_x[k] * oss1_59[k];

        t_60[k] = f_19 * nss0_60[k]
                  - f_20 * nss1_60[k]
                  + pa_x[k] * oss0_60[k]
                  - f_2 * pc_x[k] * oss1_60[k];
    }

#pragma omp simd aligned(t_61, t_62, pa_x, pc_x, nss0_61, nss0_62, nss1_61, nss1_62, oss0_61, \
                         oss0_62, oss1_61, oss1_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_19 * nss0_61[k]
                  - f_20 * nss1_61[k]
                  + pa_x[k] * oss0_61[k]
                  - f_2 * pc_x[k] * oss1_61[k];

        t_62[k] = f_19 * nss0_62[k]
                  - f_20 * nss1_62[k]
                  + pa_x[k] * oss0_62[k]
                  - f_2 * pc_x[k] * oss1_62[k];
    }

#pragma omp simd aligned(t_63, t_64, pa_x, pa_y, pc_x, pc_y, nss0_63, nss1_63, oss0_54, \
                         oss0_63, oss1_54, oss1_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_63[k] = f_19 * nss0_63[k]
                  - f_20 * nss1_63[k]
                  + pa_x[k] * oss0_63[k]
                  - f_2 * pc_x[k] * oss1_63[k];

        t_64[k] = pa_y[k] * oss0_54[k]
                  - f_2 * pc_y[k] * oss1_54[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, t_68, pa_x, pc_x, nss0_65, nss1_65, oss0_65, \
                         oss0_66, oss0_67, oss0_68, oss1_65, oss1_66, oss1_67, \
                         oss1_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_19 * nss0_65[k]
                  - f_20 * nss1_65[k]
                  + pa_x[k] * oss0_65[k]
                  - f_2 * pc_x[k] * oss1_65[k];

        t_66[k] = pa_x[k] * oss0_66[k]
                  - f_2 * pc_x[k] * oss1_66[k];

        t_67[k] = pa_x[k] * oss0_67[k]
                  - f_2 * pc_x[k] * oss1_67[k];

        t_68[k] = pa_x[k] * oss0_68[k]
                  - f_2 * pc_x[k] * oss1_68[k];
    }

#pragma omp simd aligned(t_69, t_70, t_71, t_72, pa_x, pc_x, oss0_69, oss0_70, oss0_71, \
                         oss0_72, oss1_69, oss1_70, oss1_71, oss1_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_69[k] = pa_x[k] * oss0_69[k]
                  - f_2 * pc_x[k] * oss1_69[k];

        t_70[k] = pa_x[k] * oss0_70[k]
                  - f_2 * pc_x[k] * oss1_70[k];

        t_71[k] = pa_x[k] * oss0_71[k]
                  - f_2 * pc_x[k] * oss1_71[k];

        t_72[k] = pa_x[k] * oss0_72[k]
                  - f_2 * pc_x[k] * oss1_72[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, t_76, pa_x, pc_x, oss0_73, oss0_74, oss0_75, \
                         oss0_76, oss1_73, oss1_74, oss1_75, oss1_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = pa_x[k] * oss0_73[k]
                  - f_2 * pc_x[k] * oss1_73[k];

        t_74[k] = pa_x[k] * oss0_74[k]
                  - f_2 * pc_x[k] * oss1_74[k];

        t_75[k] = pa_x[k] * oss0_75[k]
                  - f_2 * pc_x[k] * oss1_75[k];

        t_76[k] = pa_x[k] * oss0_76[k]
                  - f_2 * pc_x[k] * oss1_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, nss0_55, \
                         nss1_55, oss0_66, oss0_77, oss1_66, oss1_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = pa_x[k] * oss0_77[k]
                  - f_2 * pc_x[k] * oss1_77[k];

        t_78[k] = f_0 * nss0_55[k]
                  - f_1 * nss1_55[k]
                  + pa_y[k] * oss0_66[k]
                  - f_2 * pc_y[k] * oss1_66[k];

        t_79[k] = pa_z[k] * oss0_66[k]
                  - f_2 * pc_z[k] * oss1_66[k];
    }

#pragma omp simd aligned(t_80, t_81, pa_y, pc_y, nss0_57, nss0_58, nss1_57, nss1_58, oss0_68, \
                         oss0_69, oss1_68, oss1_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * nss0_57[k]
                  - f_4 * nss1_57[k]
                  + pa_y[k] * oss0_68[k]
                  - f_2 * pc_y[k] * oss1_68[k];

        t_81[k] = f_5 * nss0_58[k]
                  - f_6 * nss1_58[k]
                  + pa_y[k] * oss0_69[k]
                  - f_2 * pc_y[k] * oss1_69[k];
    }

#pragma omp simd aligned(t_82, t_83, pa_y, pc_y, nss0_59, nss0_60, nss1_59, nss1_60, oss0_70, \
                         oss0_71, oss1_70, oss1_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_7 * nss0_59[k]
                  - f_8 * nss1_59[k]
                  + pa_y[k] * oss0_70[k]
                  - f_2 * pc_y[k] * oss1_70[k];

        t_83[k] = f_9 * nss0_60[k]
                  - f_10 * nss1_60[k]
                  + pa_y[k] * oss0_71[k]
                  - f_2 * pc_y[k] * oss1_71[k];
    }

#pragma omp simd aligned(t_84, t_85, pa_y, pc_y, nss0_61, nss0_62, nss1_61, nss1_62, oss0_72, \
                         oss0_73, oss1_72, oss1_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_84[k] = f_11 * nss0_61[k]
                  - f_12 * nss1_61[k]
                  + pa_y[k] * oss0_72[k]
                  - f_2 * pc_y[k] * oss1_72[k];

        t_85[k] = f_13 * nss0_62[k]
                  - f_14 * nss1_62[k]
                  + pa_y[k] * oss0_73[k]
                  - f_2 * pc_y[k] * oss1_73[k];
    }

#pragma omp simd aligned(t_86, t_87, pa_y, pc_y, nss0_63, nss0_64, nss1_63, nss1_64, oss0_74, \
                         oss0_75, oss1_74, oss1_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = f_15 * nss0_63[k]
                  - f_16 * nss1_63[k]
                  + pa_y[k] * oss0_74[k]
                  - f_2 * pc_y[k] * oss1_74[k];

        t_87[k] = f_17 * nss0_64[k]
                  - f_18 * nss1_64[k]
                  + pa_y[k] * oss0_75[k]
                  - f_2 * pc_y[k] * oss1_75[k];
    }
}

static auto
compute_prim_qss_three_center_electron_repulsion_0_piece1(CSimdMatrix &buffer,
                                                          const size_t target, const size_t pa,
                                                          const size_t pc, const size_t nss0,
                                                          const size_t nss1, const size_t oss0,
                                                          const size_t oss1, const size_t ncols,
                                                          const double gamma, const double p,
                                                          const double q) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 5.5 / p;
    const auto f_1 = 5.5 * gamma / (p * q);
    const auto f_2 = gamma / q;
    const auto f_19 = 0.5 / p;
    const auto f_20 = 0.5 * gamma / (p * q);

    auto *t_88 = buffer.data(target + 88);
    auto *t_89 = buffer.data(target + 89);
    auto *t_90 = buffer.data(target + 90);

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pc_y = buffer.data(pc + 1);
    const auto *pc_z = buffer.data(pc + 2);

    const auto *nss0_65 = buffer.data(nss0 + 65);

    const auto *nss1_65 = buffer.data(nss1 + 65);

    const auto *oss0_76 = buffer.data(oss0 + 76);
    const auto *oss0_77 = buffer.data(oss0 + 77);

    const auto *oss1_76 = buffer.data(oss1 + 76);
    const auto *oss1_77 = buffer.data(oss1 + 77);

#pragma omp simd aligned(t_88, t_89, t_90, pa_y, pa_z, pc_y, pc_z, nss0_65, nss1_65, oss0_76, \
                         oss0_77, oss1_76, oss1_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_88[k] = f_19 * nss0_65[k]
                  - f_20 * nss1_65[k]
                  + pa_y[k] * oss0_76[k]
                  - f_2 * pc_y[k] * oss1_76[k];

        t_89[k] = pa_y[k] * oss0_77[k]
                  - f_2 * pc_y[k] * oss1_77[k];

        t_90[k] = f_0 * nss0_65[k]
                  - f_1 * nss1_65[k]
                  + pa_z[k] * oss0_77[k]
                  - f_2 * pc_z[k] * oss1_77[k];
    }
}

auto
compute_prim_qss_three_center_electron_repulsion_0(CSimdMatrix &buffer, const size_t target,
                                                   const size_t pa, const size_t pc,
                                                   const size_t nss0, const size_t nss1,
                                                   const size_t oss0, const size_t oss1,
                                                   const size_t ncols, const double gamma,
                                                   const double p, const double q) -> void
{
    compute_prim_qss_three_center_electron_repulsion_0_piece0(buffer, target, pa, pc, nss0,
                                                              nss1, oss0, oss1, ncols, gamma, p,
                                                              q);

    compute_prim_qss_three_center_electron_repulsion_0_piece1(buffer, target, pa, pc, nss0,
                                                              nss1, oss0, oss1, ncols, gamma, p,
                                                              q);
}

}  // namespace simdt3ceri
