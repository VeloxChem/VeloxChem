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


#include "SimdKineticEnergyVrrRecHF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_hf_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t ff_s,
                                        const size_t ff, const size_t gd, const size_t gf,
                                        const size_t hp_s, const size_t hf_s, const size_t hp,
                                        const size_t hd, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = 3.0 * beta / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = alpha / p;
    const auto f_9 = beta / p;
    const auto f_10 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff_s_0 = buffer.data(ff_s + 0);
    const auto *ff_s_10 = buffer.data(ff_s + 10);
    const auto *ff_s_16 = buffer.data(ff_s + 16);
    const auto *ff_s_20 = buffer.data(ff_s + 20);
    const auto *ff_s_29 = buffer.data(ff_s + 29);
    const auto *ff_s_36 = buffer.data(ff_s + 36);
    const auto *ff_s_59 = buffer.data(ff_s + 59);
    const auto *ff_s_66 = buffer.data(ff_s + 66);
    const auto *ff_s_79 = buffer.data(ff_s + 79);
    const auto *ff_s_86 = buffer.data(ff_s + 86);
    const auto *ff_s_99 = buffer.data(ff_s + 99);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_6 = buffer.data(gd + 6);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_11 = buffer.data(gd + 11);
    const auto *gd_12 = buffer.data(gd + 12);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_27 = buffer.data(gd + 27);
    const auto *gd_28 = buffer.data(gd + 28);
    const auto *gd_29 = buffer.data(gd + 29);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_35 = buffer.data(gd + 35);
    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_39 = buffer.data(gd + 39);
    const auto *gd_41 = buffer.data(gd + 41);
    const auto *gd_46 = buffer.data(gd + 46);
    const auto *gd_47 = buffer.data(gd + 47);
    const auto *gd_51 = buffer.data(gd + 51);
    const auto *gd_52 = buffer.data(gd + 52);
    const auto *gd_57 = buffer.data(gd + 57);
    const auto *gd_59 = buffer.data(gd + 59);
    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_63 = buffer.data(gd + 63);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);

    const auto *hp_s_0 = buffer.data(hp_s + 0);
    const auto *hp_s_1 = buffer.data(hp_s + 1);
    const auto *hp_s_2 = buffer.data(hp_s + 2);
    const auto *hp_s_8 = buffer.data(hp_s + 8);
    const auto *hp_s_11 = buffer.data(hp_s + 11);
    const auto *hp_s_16 = buffer.data(hp_s + 16);
    const auto *hp_s_17 = buffer.data(hp_s + 17);
    const auto *hp_s_20 = buffer.data(hp_s + 20);
    const auto *hp_s_28 = buffer.data(hp_s + 28);
    const auto *hp_s_29 = buffer.data(hp_s + 29);

    const auto *hf_s_0 = buffer.data(hf_s + 0);
    const auto *hf_s_1 = buffer.data(hf_s + 1);
    const auto *hf_s_2 = buffer.data(hf_s + 2);
    const auto *hf_s_3 = buffer.data(hf_s + 3);
    const auto *hf_s_4 = buffer.data(hf_s + 4);
    const auto *hf_s_5 = buffer.data(hf_s + 5);
    const auto *hf_s_6 = buffer.data(hf_s + 6);
    const auto *hf_s_7 = buffer.data(hf_s + 7);
    const auto *hf_s_8 = buffer.data(hf_s + 8);
    const auto *hf_s_9 = buffer.data(hf_s + 9);
    const auto *hf_s_10 = buffer.data(hf_s + 10);
    const auto *hf_s_11 = buffer.data(hf_s + 11);
    const auto *hf_s_12 = buffer.data(hf_s + 12);
    const auto *hf_s_13 = buffer.data(hf_s + 13);
    const auto *hf_s_14 = buffer.data(hf_s + 14);
    const auto *hf_s_15 = buffer.data(hf_s + 15);
    const auto *hf_s_16 = buffer.data(hf_s + 16);
    const auto *hf_s_17 = buffer.data(hf_s + 17);
    const auto *hf_s_18 = buffer.data(hf_s + 18);
    const auto *hf_s_19 = buffer.data(hf_s + 19);
    const auto *hf_s_20 = buffer.data(hf_s + 20);
    const auto *hf_s_21 = buffer.data(hf_s + 21);
    const auto *hf_s_22 = buffer.data(hf_s + 22);
    const auto *hf_s_23 = buffer.data(hf_s + 23);
    const auto *hf_s_24 = buffer.data(hf_s + 24);
    const auto *hf_s_25 = buffer.data(hf_s + 25);
    const auto *hf_s_26 = buffer.data(hf_s + 26);
    const auto *hf_s_27 = buffer.data(hf_s + 27);
    const auto *hf_s_28 = buffer.data(hf_s + 28);
    const auto *hf_s_29 = buffer.data(hf_s + 29);
    const auto *hf_s_30 = buffer.data(hf_s + 30);
    const auto *hf_s_31 = buffer.data(hf_s + 31);
    const auto *hf_s_32 = buffer.data(hf_s + 32);
    const auto *hf_s_33 = buffer.data(hf_s + 33);
    const auto *hf_s_34 = buffer.data(hf_s + 34);
    const auto *hf_s_35 = buffer.data(hf_s + 35);
    const auto *hf_s_36 = buffer.data(hf_s + 36);
    const auto *hf_s_37 = buffer.data(hf_s + 37);
    const auto *hf_s_38 = buffer.data(hf_s + 38);
    const auto *hf_s_39 = buffer.data(hf_s + 39);
    const auto *hf_s_40 = buffer.data(hf_s + 40);
    const auto *hf_s_41 = buffer.data(hf_s + 41);
    const auto *hf_s_42 = buffer.data(hf_s + 42);
    const auto *hf_s_43 = buffer.data(hf_s + 43);
    const auto *hf_s_44 = buffer.data(hf_s + 44);
    const auto *hf_s_45 = buffer.data(hf_s + 45);
    const auto *hf_s_46 = buffer.data(hf_s + 46);
    const auto *hf_s_47 = buffer.data(hf_s + 47);
    const auto *hf_s_48 = buffer.data(hf_s + 48);
    const auto *hf_s_49 = buffer.data(hf_s + 49);
    const auto *hf_s_50 = buffer.data(hf_s + 50);
    const auto *hf_s_51 = buffer.data(hf_s + 51);
    const auto *hf_s_52 = buffer.data(hf_s + 52);
    const auto *hf_s_53 = buffer.data(hf_s + 53);
    const auto *hf_s_54 = buffer.data(hf_s + 54);
    const auto *hf_s_55 = buffer.data(hf_s + 55);
    const auto *hf_s_56 = buffer.data(hf_s + 56);
    const auto *hf_s_57 = buffer.data(hf_s + 57);
    const auto *hf_s_58 = buffer.data(hf_s + 58);
    const auto *hf_s_59 = buffer.data(hf_s + 59);
    const auto *hf_s_60 = buffer.data(hf_s + 60);
    const auto *hf_s_61 = buffer.data(hf_s + 61);
    const auto *hf_s_62 = buffer.data(hf_s + 62);
    const auto *hf_s_63 = buffer.data(hf_s + 63);
    const auto *hf_s_64 = buffer.data(hf_s + 64);
    const auto *hf_s_65 = buffer.data(hf_s + 65);
    const auto *hf_s_66 = buffer.data(hf_s + 66);
    const auto *hf_s_67 = buffer.data(hf_s + 67);
    const auto *hf_s_68 = buffer.data(hf_s + 68);
    const auto *hf_s_69 = buffer.data(hf_s + 69);
    const auto *hf_s_70 = buffer.data(hf_s + 70);
    const auto *hf_s_71 = buffer.data(hf_s + 71);
    const auto *hf_s_72 = buffer.data(hf_s + 72);
    const auto *hf_s_73 = buffer.data(hf_s + 73);
    const auto *hf_s_74 = buffer.data(hf_s + 74);
    const auto *hf_s_75 = buffer.data(hf_s + 75);
    const auto *hf_s_76 = buffer.data(hf_s + 76);
    const auto *hf_s_77 = buffer.data(hf_s + 77);
    const auto *hf_s_78 = buffer.data(hf_s + 78);
    const auto *hf_s_79 = buffer.data(hf_s + 79);
    const auto *hf_s_80 = buffer.data(hf_s + 80);
    const auto *hf_s_81 = buffer.data(hf_s + 81);
    const auto *hf_s_82 = buffer.data(hf_s + 82);
    const auto *hf_s_83 = buffer.data(hf_s + 83);
    const auto *hf_s_84 = buffer.data(hf_s + 84);
    const auto *hf_s_85 = buffer.data(hf_s + 85);
    const auto *hf_s_86 = buffer.data(hf_s + 86);
    const auto *hf_s_87 = buffer.data(hf_s + 87);
    const auto *hf_s_88 = buffer.data(hf_s + 88);
    const auto *hf_s_89 = buffer.data(hf_s + 89);
    const auto *hf_s_90 = buffer.data(hf_s + 90);
    const auto *hf_s_91 = buffer.data(hf_s + 91);
    const auto *hf_s_92 = buffer.data(hf_s + 92);
    const auto *hf_s_93 = buffer.data(hf_s + 93);
    const auto *hf_s_94 = buffer.data(hf_s + 94);
    const auto *hf_s_95 = buffer.data(hf_s + 95);
    const auto *hf_s_96 = buffer.data(hf_s + 96);
    const auto *hf_s_97 = buffer.data(hf_s + 97);
    const auto *hf_s_98 = buffer.data(hf_s + 98);
    const auto *hf_s_99 = buffer.data(hf_s + 99);
    const auto *hf_s_100 = buffer.data(hf_s + 100);
    const auto *hf_s_101 = buffer.data(hf_s + 101);
    const auto *hf_s_102 = buffer.data(hf_s + 102);
    const auto *hf_s_103 = buffer.data(hf_s + 103);
    const auto *hf_s_104 = buffer.data(hf_s + 104);

    const auto *hp_0 = buffer.data(hp + 0);
    const auto *hp_1 = buffer.data(hp + 1);
    const auto *hp_2 = buffer.data(hp + 2);
    const auto *hp_8 = buffer.data(hp + 8);
    const auto *hp_11 = buffer.data(hp + 11);
    const auto *hp_16 = buffer.data(hp + 16);
    const auto *hp_17 = buffer.data(hp + 17);
    const auto *hp_20 = buffer.data(hp + 20);
    const auto *hp_28 = buffer.data(hp + 28);
    const auto *hp_29 = buffer.data(hp + 29);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_2 = buffer.data(hd + 2);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_7 = buffer.data(hd + 7);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_14 = buffer.data(hd + 14);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_19 = buffer.data(hd + 19);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_32 = buffer.data(hd + 32);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_37 = buffer.data(hd + 37);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_56 = buffer.data(hd + 56);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
    const auto *hd_61 = buffer.data(hd + 61);
    const auto *hd_63 = buffer.data(hd + 63);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, gd_0, hp_s_0, hf_s_0, hf_s_1, \
                         hf_s_2, hp_0, hd_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gd_0[k]
                 - f_1 * hp_s_0[k]
                 + f_2 * hf_s_0[k]
                 + f_3 * hp_0[k]
                 + pb_x[k] * hd_0[k];

        t_1[k] = f_2 * hf_s_1[k]
                 + pb_y[k] * hd_0[k];

        t_2[k] = f_2 * hf_s_2[k]
                 + pb_z[k] * hd_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, gd_3, gd_5, hf_s_3, hf_s_4, hf_s_5, hd_2, \
                         hd_3, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * gd_3[k]
                 + f_2 * hf_s_3[k]
                 + pb_x[k] * hd_3[k];

        t_4[k] = f_2 * hf_s_4[k]
                 + pb_y[k] * hd_2[k];

        t_5[k] = f_0 * gd_5[k]
                 + f_2 * hf_s_5[k]
                 + pb_x[k] * hd_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pb_z, hp_s_1, hp_s_2, hf_s_6, hf_s_7, \
                         hf_s_8, hf_s_9, hp_1, hp_2, hd_3, hd_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * hp_s_1[k]
                 + f_2 * hf_s_6[k]
                 + f_3 * hp_1[k]
                 + pb_y[k] * hd_3[k];

        t_7[k] = f_2 * hf_s_7[k]
                 + pb_z[k] * hd_3[k];

        t_8[k] = f_2 * hf_s_8[k]
                 + pb_y[k] * hd_5[k];

        t_9[k] = -f_1 * hp_s_2[k]
                 + f_2 * hf_s_9[k]
                 + f_3 * hp_2[k]
                 + pb_z[k] * hd_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, gd_0, gf_0, hf_s_10, hf_s_11, \
                         hf_s_12, hd_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * gf_0[k]
                  + f_2 * hf_s_10[k];

        t_11[k] = f_4 * gd_0[k]
                  + f_2 * hf_s_11[k]
                  + pb_y[k] * hd_6[k];

        t_12[k] = f_2 * hf_s_12[k]
                  + pb_z[k] * hd_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, gd_9, gf_5, hf_s_13, hf_s_14, \
                         hf_s_15, hd_7, hd_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * gd_9[k]
                  + f_2 * hf_s_13[k]
                  + pb_x[k] * hd_9[k];

        t_14[k] = f_2 * hf_s_14[k]
                  + pb_z[k] * hd_7[k];

        t_15[k] = pa_y[k] * gf_5[k]
                  + f_2 * hf_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, pb_z, ff_s_16, ff_16, gd_5, gf_16, \
                         hf_s_16, hf_s_17, hf_s_18, hd_9, hd_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * ff_s_16[k]
                  + f_7 * ff_16[k]
                  + pa_x[k] * gf_16[k]
                  + f_2 * hf_s_16[k];

        t_17[k] = f_2 * hf_s_17[k]
                  + pb_z[k] * hd_9[k];

        t_18[k] = f_4 * gd_5[k]
                  + f_2 * hf_s_18[k]
                  + pb_y[k] * hd_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, gd_0, gf_0, gf_9, \
                         hf_s_19, hf_s_20, hf_s_21, hf_s_22, hd_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * gf_9[k]
                  + f_2 * hf_s_19[k];

        t_20[k] = pa_z[k] * gf_0[k]
                  + f_2 * hf_s_20[k];

        t_21[k] = f_2 * hf_s_21[k]
                  + pb_y[k] * hd_12[k];

        t_22[k] = f_4 * gd_0[k]
                  + f_2 * hf_s_22[k]
                  + pb_z[k] * hd_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, gd_17, gf_3, gf_6, hf_s_23, \
                         hf_s_24, hf_s_25, hf_s_26, hd_14, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * gf_3[k]
                  + f_2 * hf_s_23[k];

        t_24[k] = f_2 * hf_s_24[k]
                  + pb_y[k] * hd_14[k];

        t_25[k] = f_5 * gd_17[k]
                  + f_2 * hf_s_25[k]
                  + pb_x[k] * hd_17[k];

        t_26[k] = pa_z[k] * gf_6[k]
                  + f_2 * hf_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, ff_s_29, ff_29, gf_29, hp_s_8, hf_s_27, \
                         hf_s_28, hf_s_29, hp_8, hd_16, hd_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_8 * hp_s_8[k]
                  + f_2 * hf_s_27[k]
                  + f_4 * hp_8[k]
                  + pb_y[k] * hd_16[k];

        t_28[k] = f_2 * hf_s_28[k]
                  + pb_y[k] * hd_17[k];

        t_29[k] = -f_6 * ff_s_29[k]
                  + f_7 * ff_29[k]
                  + pa_x[k] * gf_29[k]
                  + f_2 * hf_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_y, pb_z, ff_s_0, ff_0, gd_6, gf_10, \
                         hf_s_30, hf_s_31, hf_s_32, hd_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_9 * ff_s_0[k]
                  + f_4 * ff_0[k]
                  + pa_y[k] * gf_10[k]
                  + f_2 * hf_s_30[k];

        t_31[k] = f_3 * gd_6[k]
                  + f_2 * hf_s_31[k]
                  + pb_y[k] * hd_18[k];

        t_32[k] = f_2 * hf_s_32[k]
                  + pb_z[k] * hd_18[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_z, gd_21, gd_23, hf_s_33, hf_s_34, \
                         hf_s_35, hd_19, hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * gd_21[k]
                  + f_2 * hf_s_33[k]
                  + pb_x[k] * hd_21[k];

        t_34[k] = f_2 * hf_s_34[k]
                  + pb_z[k] * hd_19[k];

        t_35[k] = f_7 * gd_23[k]
                  + f_2 * hf_s_35[k]
                  + pb_x[k] * hd_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_y, pb_z, ff_s_36, ff_36, gd_11, gf_36, \
                         hf_s_36, hf_s_37, hf_s_38, hd_21, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_10 * ff_s_36[k]
                  + f_3 * ff_36[k]
                  + pa_x[k] * gf_36[k]
                  + f_2 * hf_s_36[k];

        t_37[k] = f_2 * hf_s_37[k]
                  + pb_z[k] * hd_21[k];

        t_38[k] = f_3 * gd_11[k]
                  + f_2 * hf_s_38[k]
                  + pb_y[k] * hd_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pa_z, pb_z, gf_11, gf_20, hp_s_11, hf_s_39, \
                         hf_s_40, hf_s_41, hp_11, hd_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_1 * hp_s_11[k]
                  + f_2 * hf_s_39[k]
                  + f_3 * hp_11[k]
                  + pb_z[k] * hd_23[k];

        t_40[k] = pa_y[k] * gf_20[k]
                  + f_2 * hf_s_40[k];

        t_41[k] = pa_z[k] * gf_11[k]
                  + f_2 * hf_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, gd_28, gf_13, gf_22, gf_25, \
                         hf_s_42, hf_s_43, hf_s_44, hf_s_45, hd_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * gf_22[k]
                  + f_2 * hf_s_42[k];

        t_43[k] = pa_z[k] * gf_13[k]
                  + f_2 * hf_s_43[k];

        t_44[k] = f_7 * gd_28[k]
                  + f_2 * hf_s_44[k]
                  + pb_x[k] * hd_28[k];

        t_45[k] = pa_y[k] * gf_25[k]
                  + f_2 * hf_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_z, pb_y, pb_z, gd_9, gd_17, gf_16, hf_s_46, \
                         hf_s_47, hf_s_48, hd_27, hd_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * gf_16[k]
                  + f_2 * hf_s_46[k];

        t_47[k] = f_4 * gd_9[k]
                  + f_2 * hf_s_47[k]
                  + pb_z[k] * hd_27[k];

        t_48[k] = f_4 * gd_17[k]
                  + f_2 * hf_s_48[k]
                  + pb_y[k] * hd_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_y, ff_s_0, ff_0, gf_20, gf_29, \
                         hf_s_49, hf_s_50, hf_s_51, hd_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_y[k] * gf_29[k]
                  + f_2 * hf_s_49[k];

        t_50[k] = -f_9 * ff_s_0[k]
                  + f_4 * ff_0[k]
                  + pa_z[k] * gf_20[k]
                  + f_2 * hf_s_50[k];

        t_51[k] = f_2 * hf_s_51[k]
                  + pb_y[k] * hd_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, pb_z, gd_12, gd_33, hf_s_52, hf_s_53, \
                         hf_s_54, hd_30, hd_32, hd_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * gd_12[k]
                  + f_2 * hf_s_52[k]
                  + pb_z[k] * hd_30[k];

        t_53[k] = f_7 * gd_33[k]
                  + f_2 * hf_s_53[k]
                  + pb_x[k] * hd_33[k];

        t_54[k] = f_2 * hf_s_54[k]
                  + pb_y[k] * hd_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, gd_35, hp_s_16, hp_s_17, hf_s_55, \
                         hf_s_56, hf_s_57, hp_16, hp_17, hd_33, hd_34, \
                         hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_7 * gd_35[k]
                  + f_2 * hf_s_55[k]
                  + pb_x[k] * hd_35[k];

        t_56[k] = -f_1 * hp_s_16[k]
                  + f_2 * hf_s_56[k]
                  + f_3 * hp_16[k]
                  + pb_y[k] * hd_33[k];

        t_57[k] = -f_8 * hp_s_17[k]
                  + f_2 * hf_s_57[k]
                  + f_4 * hp_17[k]
                  + pb_y[k] * hd_34[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, ff_s_10, ff_s_59, ff_10, ff_59, \
                         gf_30, gf_59, hf_s_58, hf_s_59, hf_s_60, \
                         hd_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * hf_s_58[k]
                  + pb_y[k] * hd_35[k];

        t_59[k] = -f_10 * ff_s_59[k]
                  + f_3 * ff_59[k]
                  + pa_x[k] * gf_59[k]
                  + f_2 * hf_s_59[k];

        t_60[k] = -f_10 * ff_s_10[k]
                  + f_3 * ff_10[k]
                  + pa_y[k] * gf_30[k]
                  + f_2 * hf_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_x, pb_y, pb_z, gd_18, gd_39, hf_s_61, \
                         hf_s_62, hf_s_63, hf_s_64, hd_36, hd_37, \
                         hd_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_7 * gd_18[k]
                  + f_2 * hf_s_61[k]
                  + pb_y[k] * hd_36[k];

        t_62[k] = f_2 * hf_s_62[k]
                  + pb_z[k] * hd_36[k];

        t_63[k] = f_3 * gd_39[k]
                  + f_2 * hf_s_63[k]
                  + pb_x[k] * hd_39[k];

        t_64[k] = f_2 * hf_s_64[k]
                  + pb_z[k] * hd_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_x, pb_x, pb_z, ff_s_66, ff_66, gd_41, gf_66, \
                         hf_s_65, hf_s_66, hf_s_67, hd_39, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_3 * gd_41[k]
                  + f_2 * hf_s_65[k]
                  + pb_x[k] * hd_41[k];

        t_66[k] = -f_9 * ff_s_66[k]
                  + f_4 * ff_66[k]
                  + pa_x[k] * gf_66[k]
                  + f_2 * hf_s_66[k];

        t_67[k] = f_2 * hf_s_67[k]
                  + pb_z[k] * hd_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_z, pb_y, pb_z, gd_23, gf_30, hp_s_20, hf_s_68, \
                         hf_s_69, hf_s_70, hp_20, hd_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_7 * gd_23[k]
                  + f_2 * hf_s_68[k]
                  + pb_y[k] * hd_41[k];

        t_69[k] = -f_1 * hp_s_20[k]
                  + f_2 * hf_s_69[k]
                  + f_3 * hp_20[k]
                  + pb_z[k] * hd_41[k];

        t_70[k] = pa_z[k] * gf_30[k]
                  + f_2 * hf_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_z, pb_z, gd_18, gf_31, gf_33, hf_s_71, hf_s_72, \
                         hf_s_73, hd_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_z[k] * gf_31[k]
                  + f_2 * hf_s_71[k];

        t_72[k] = f_4 * gd_18[k]
                  + f_2 * hf_s_72[k]
                  + pb_z[k] * hd_42[k];

        t_73[k] = pa_z[k] * gf_33[k]
                  + f_2 * hf_s_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_z, pb_x, gd_46, gd_47, gf_36, hf_s_74, hf_s_75, \
                         hf_s_76, hd_46, hd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_3 * gd_46[k]
                  + f_2 * hf_s_74[k]
                  + pb_x[k] * hd_46[k];

        t_75[k] = f_3 * gd_47[k]
                  + f_2 * hf_s_75[k]
                  + pb_x[k] * hd_47[k];

        t_76[k] = pa_z[k] * gf_36[k]
                  + f_2 * hf_s_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pb_y, pb_z, ff_s_79, ff_79, gd_21, gd_29, \
                         gf_79, hf_s_77, hf_s_78, hf_s_79, hd_45, \
                         hd_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_4 * gd_21[k]
                  + f_2 * hf_s_77[k]
                  + pb_z[k] * hd_45[k];

        t_78[k] = f_3 * gd_29[k]
                  + f_2 * hf_s_78[k]
                  + pb_y[k] * hd_47[k];

        t_79[k] = -f_9 * ff_s_79[k]
                  + f_4 * ff_79[k]
                  + pa_x[k] * gf_79[k]
                  + f_2 * hf_s_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_y, pb_y, gd_30, gf_50, gf_52, hf_s_80, hf_s_81, \
                         hf_s_82, hd_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * gf_50[k]
                  + f_2 * hf_s_80[k];

        t_81[k] = f_4 * gd_30[k]
                  + f_2 * hf_s_81[k]
                  + pb_y[k] * hd_48[k];

        t_82[k] = pa_y[k] * gf_52[k]
                  + f_2 * hf_s_82[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_y, pb_x, gd_51, gd_52, gf_55, hf_s_83, hf_s_84, \
                         hf_s_85, hd_51, hd_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_3 * gd_51[k]
                  + f_2 * hf_s_83[k]
                  + pb_x[k] * hd_51[k];

        t_84[k] = f_3 * gd_52[k]
                  + f_2 * hf_s_84[k]
                  + pb_x[k] * hd_52[k];

        t_85[k] = pa_y[k] * gf_55[k]
                  + f_2 * hf_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_x, pb_y, pb_z, ff_s_86, ff_86, gd_27, gd_35, \
                         gf_86, hf_s_86, hf_s_87, hf_s_88, hd_51, \
                         hd_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_9 * ff_s_86[k]
                  + f_4 * ff_86[k]
                  + pa_x[k] * gf_86[k]
                  + f_2 * hf_s_86[k];

        t_87[k] = f_3 * gd_27[k]
                  + f_2 * hf_s_87[k]
                  + pb_z[k] * hd_51[k];

        t_88[k] = f_4 * gd_35[k]
                  + f_2 * hf_s_88[k]
                  + pb_y[k] * hd_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_y, pa_z, pb_y, ff_s_20, ff_20, gf_50, gf_59, \
                         hf_s_89, hf_s_90, hf_s_91, hd_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * gf_59[k]
                  + f_2 * hf_s_89[k];

        t_90[k] = -f_10 * ff_s_20[k]
                  + f_3 * ff_20[k]
                  + pa_z[k] * gf_50[k]
                  + f_2 * hf_s_90[k];

        t_91[k] = f_2 * hf_s_91[k]
                  + pb_y[k] * hd_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pb_x, pb_y, pb_z, gd_30, gd_57, hf_s_92, hf_s_93, \
                         hf_s_94, hd_54, hd_56, hd_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_7 * gd_30[k]
                  + f_2 * hf_s_92[k]
                  + pb_z[k] * hd_54[k];

        t_93[k] = f_3 * gd_57[k]
                  + f_2 * hf_s_93[k]
                  + pb_x[k] * hd_57[k];

        t_94[k] = f_2 * hf_s_94[k]
                  + pb_y[k] * hd_56[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, pb_y, gd_59, hp_s_28, hp_s_29, hf_s_95, \
                         hf_s_96, hf_s_97, hp_28, hp_29, hd_57, hd_58, \
                         hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_3 * gd_59[k]
                  + f_2 * hf_s_95[k]
                  + pb_x[k] * hd_59[k];

        t_96[k] = -f_1 * hp_s_28[k]
                  + f_2 * hf_s_96[k]
                  + f_3 * hp_28[k]
                  + pb_y[k] * hd_57[k];

        t_97[k] = -f_8 * hp_s_29[k]
                  + f_2 * hf_s_97[k]
                  + f_4 * hp_29[k]
                  + pb_y[k] * hd_58[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_x, pb_y, ff_s_99, ff_99, gd_60, gf_99, gf_100, \
                         hf_s_98, hf_s_99, hf_s_100, hd_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_2 * hf_s_98[k]
                  + pb_y[k] * hd_59[k];

        t_99[k] = -f_9 * ff_s_99[k]
                  + f_4 * ff_99[k]
                  + pa_x[k] * gf_99[k]
                  + f_2 * hf_s_99[k];

        t_100[k] = f_7 * gd_60[k]
                   + pa_x[k] * gf_100[k]
                   + f_2 * hf_s_100[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, pb_y, pb_z, gd_36, gd_63, hf_s_101, \
                         hf_s_102, hf_s_103, hf_s_104, hd_60, hd_61, \
                         hd_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_5 * gd_36[k]
                   + f_2 * hf_s_101[k]
                   + pb_y[k] * hd_60[k];

        t_102[k] = f_2 * hf_s_102[k]
                   + pb_z[k] * hd_60[k];

        t_103[k] = f_4 * gd_63[k]
                   + f_2 * hf_s_103[k]
                   + pb_x[k] * hd_63[k];

        t_104[k] = f_2 * hf_s_104[k]
                   + pb_z[k] * hd_61[k];
    }
}

static auto
compute_prim_hf_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t ff_s,
                                        const size_t ff, const size_t gd, const size_t gf,
                                        const size_t hp_s, const size_t hf_s, const size_t hp,
                                        const size_t hd, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 / p;
    const auto f_6 = 3.0 * beta / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = alpha / p;
    const auto f_9 = beta / p;
    const auto f_10 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff_s_66 = buffer.data(ff_s + 66);
    const auto *ff_s_76 = buffer.data(ff_s + 76);
    const auto *ff_s_79 = buffer.data(ff_s + 79);
    const auto *ff_s_89 = buffer.data(ff_s + 89);
    const auto *ff_s_99 = buffer.data(ff_s + 99);

    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *gd_36 = buffer.data(gd + 36);
    const auto *gd_42 = buffer.data(gd + 42);
    const auto *gd_48 = buffer.data(gd + 48);
    const auto *gd_54 = buffer.data(gd + 54);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_69 = buffer.data(gd + 69);
    const auto *gd_70 = buffer.data(gd + 70);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_81 = buffer.data(gd + 81);
    const auto *gd_82 = buffer.data(gd + 82);
    const auto *gd_83 = buffer.data(gd + 83);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *hp_s_45 = buffer.data(hp_s + 45);
    const auto *hp_s_46 = buffer.data(hp_s + 46);
    const auto *hp_s_47 = buffer.data(hp_s + 47);
    const auto *hp_s_50 = buffer.data(hp_s + 50);
    const auto *hp_s_51 = buffer.data(hp_s + 51);
    const auto *hp_s_52 = buffer.data(hp_s + 52);
    const auto *hp_s_53 = buffer.data(hp_s + 53);
    const auto *hp_s_54 = buffer.data(hp_s + 54);
    const auto *hp_s_55 = buffer.data(hp_s + 55);
    const auto *hp_s_56 = buffer.data(hp_s + 56);
    const auto *hp_s_60 = buffer.data(hp_s + 60);
    const auto *hp_s_61 = buffer.data(hp_s + 61);
    const auto *hp_s_62 = buffer.data(hp_s + 62);

    const auto *hf_s_105 = buffer.data(hf_s + 105);
    const auto *hf_s_106 = buffer.data(hf_s + 106);
    const auto *hf_s_107 = buffer.data(hf_s + 107);
    const auto *hf_s_108 = buffer.data(hf_s + 108);
    const auto *hf_s_109 = buffer.data(hf_s + 109);
    const auto *hf_s_110 = buffer.data(hf_s + 110);
    const auto *hf_s_111 = buffer.data(hf_s + 111);
    const auto *hf_s_112 = buffer.data(hf_s + 112);
    const auto *hf_s_113 = buffer.data(hf_s + 113);
    const auto *hf_s_114 = buffer.data(hf_s + 114);
    const auto *hf_s_115 = buffer.data(hf_s + 115);
    const auto *hf_s_116 = buffer.data(hf_s + 116);
    const auto *hf_s_117 = buffer.data(hf_s + 117);
    const auto *hf_s_118 = buffer.data(hf_s + 118);
    const auto *hf_s_119 = buffer.data(hf_s + 119);
    const auto *hf_s_120 = buffer.data(hf_s + 120);
    const auto *hf_s_121 = buffer.data(hf_s + 121);
    const auto *hf_s_122 = buffer.data(hf_s + 122);
    const auto *hf_s_123 = buffer.data(hf_s + 123);
    const auto *hf_s_124 = buffer.data(hf_s + 124);
    const auto *hf_s_125 = buffer.data(hf_s + 125);
    const auto *hf_s_126 = buffer.data(hf_s + 126);
    const auto *hf_s_127 = buffer.data(hf_s + 127);
    const auto *hf_s_128 = buffer.data(hf_s + 128);
    const auto *hf_s_129 = buffer.data(hf_s + 129);
    const auto *hf_s_130 = buffer.data(hf_s + 130);
    const auto *hf_s_131 = buffer.data(hf_s + 131);
    const auto *hf_s_132 = buffer.data(hf_s + 132);
    const auto *hf_s_133 = buffer.data(hf_s + 133);
    const auto *hf_s_134 = buffer.data(hf_s + 134);
    const auto *hf_s_135 = buffer.data(hf_s + 135);
    const auto *hf_s_136 = buffer.data(hf_s + 136);
    const auto *hf_s_137 = buffer.data(hf_s + 137);
    const auto *hf_s_138 = buffer.data(hf_s + 138);
    const auto *hf_s_139 = buffer.data(hf_s + 139);
    const auto *hf_s_140 = buffer.data(hf_s + 140);
    const auto *hf_s_141 = buffer.data(hf_s + 141);
    const auto *hf_s_142 = buffer.data(hf_s + 142);
    const auto *hf_s_143 = buffer.data(hf_s + 143);
    const auto *hf_s_144 = buffer.data(hf_s + 144);
    const auto *hf_s_145 = buffer.data(hf_s + 145);
    const auto *hf_s_146 = buffer.data(hf_s + 146);
    const auto *hf_s_147 = buffer.data(hf_s + 147);
    const auto *hf_s_148 = buffer.data(hf_s + 148);
    const auto *hf_s_149 = buffer.data(hf_s + 149);
    const auto *hf_s_150 = buffer.data(hf_s + 150);
    const auto *hf_s_151 = buffer.data(hf_s + 151);
    const auto *hf_s_152 = buffer.data(hf_s + 152);
    const auto *hf_s_153 = buffer.data(hf_s + 153);
    const auto *hf_s_154 = buffer.data(hf_s + 154);
    const auto *hf_s_155 = buffer.data(hf_s + 155);
    const auto *hf_s_156 = buffer.data(hf_s + 156);
    const auto *hf_s_157 = buffer.data(hf_s + 157);
    const auto *hf_s_158 = buffer.data(hf_s + 158);
    const auto *hf_s_159 = buffer.data(hf_s + 159);
    const auto *hf_s_160 = buffer.data(hf_s + 160);
    const auto *hf_s_161 = buffer.data(hf_s + 161);
    const auto *hf_s_162 = buffer.data(hf_s + 162);
    const auto *hf_s_163 = buffer.data(hf_s + 163);
    const auto *hf_s_164 = buffer.data(hf_s + 164);
    const auto *hf_s_165 = buffer.data(hf_s + 165);
    const auto *hf_s_166 = buffer.data(hf_s + 166);
    const auto *hf_s_167 = buffer.data(hf_s + 167);
    const auto *hf_s_168 = buffer.data(hf_s + 168);
    const auto *hf_s_169 = buffer.data(hf_s + 169);
    const auto *hf_s_170 = buffer.data(hf_s + 170);
    const auto *hf_s_171 = buffer.data(hf_s + 171);
    const auto *hf_s_172 = buffer.data(hf_s + 172);
    const auto *hf_s_173 = buffer.data(hf_s + 173);
    const auto *hf_s_174 = buffer.data(hf_s + 174);
    const auto *hf_s_175 = buffer.data(hf_s + 175);
    const auto *hf_s_176 = buffer.data(hf_s + 176);
    const auto *hf_s_177 = buffer.data(hf_s + 177);
    const auto *hf_s_178 = buffer.data(hf_s + 178);
    const auto *hf_s_179 = buffer.data(hf_s + 179);
    const auto *hf_s_180 = buffer.data(hf_s + 180);
    const auto *hf_s_181 = buffer.data(hf_s + 181);
    const auto *hf_s_182 = buffer.data(hf_s + 182);
    const auto *hf_s_183 = buffer.data(hf_s + 183);
    const auto *hf_s_184 = buffer.data(hf_s + 184);
    const auto *hf_s_185 = buffer.data(hf_s + 185);
    const auto *hf_s_186 = buffer.data(hf_s + 186);
    const auto *hf_s_187 = buffer.data(hf_s + 187);
    const auto *hf_s_188 = buffer.data(hf_s + 188);
    const auto *hf_s_189 = buffer.data(hf_s + 189);
    const auto *hf_s_190 = buffer.data(hf_s + 190);
    const auto *hf_s_191 = buffer.data(hf_s + 191);
    const auto *hf_s_192 = buffer.data(hf_s + 192);
    const auto *hf_s_193 = buffer.data(hf_s + 193);
    const auto *hf_s_194 = buffer.data(hf_s + 194);
    const auto *hf_s_195 = buffer.data(hf_s + 195);
    const auto *hf_s_196 = buffer.data(hf_s + 196);
    const auto *hf_s_197 = buffer.data(hf_s + 197);
    const auto *hf_s_198 = buffer.data(hf_s + 198);
    const auto *hf_s_199 = buffer.data(hf_s + 199);
    const auto *hf_s_200 = buffer.data(hf_s + 200);
    const auto *hf_s_201 = buffer.data(hf_s + 201);
    const auto *hf_s_202 = buffer.data(hf_s + 202);
    const auto *hf_s_203 = buffer.data(hf_s + 203);
    const auto *hf_s_204 = buffer.data(hf_s + 204);
    const auto *hf_s_205 = buffer.data(hf_s + 205);
    const auto *hf_s_206 = buffer.data(hf_s + 206);
    const auto *hf_s_207 = buffer.data(hf_s + 207);
    const auto *hf_s_208 = buffer.data(hf_s + 208);
    const auto *hf_s_209 = buffer.data(hf_s + 209);

    const auto *hp_45 = buffer.data(hp + 45);
    const auto *hp_46 = buffer.data(hp + 46);
    const auto *hp_47 = buffer.data(hp + 47);
    const auto *hp_50 = buffer.data(hp + 50);
    const auto *hp_51 = buffer.data(hp + 51);
    const auto *hp_52 = buffer.data(hp + 52);
    const auto *hp_53 = buffer.data(hp + 53);
    const auto *hp_54 = buffer.data(hp + 54);
    const auto *hp_55 = buffer.data(hp + 55);
    const auto *hp_56 = buffer.data(hp + 56);
    const auto *hp_60 = buffer.data(hp + 60);
    const auto *hp_61 = buffer.data(hp + 61);
    const auto *hp_62 = buffer.data(hp + 62);

    const auto *hd_63 = buffer.data(hd + 63);
    const auto *hd_65 = buffer.data(hd + 65);
    const auto *hd_66 = buffer.data(hd + 66);
    const auto *hd_70 = buffer.data(hd + 70);
    const auto *hd_71 = buffer.data(hd + 71);
    const auto *hd_72 = buffer.data(hd + 72);
    const auto *hd_75 = buffer.data(hd + 75);
    const auto *hd_76 = buffer.data(hd + 76);
    const auto *hd_77 = buffer.data(hd + 77);
    const auto *hd_78 = buffer.data(hd + 78);
    const auto *hd_81 = buffer.data(hd + 81);
    const auto *hd_82 = buffer.data(hd + 82);
    const auto *hd_84 = buffer.data(hd + 84);
    const auto *hd_86 = buffer.data(hd + 86);
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_94 = buffer.data(hd + 94);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_98 = buffer.data(hd + 98);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_103 = buffer.data(hd + 103);
    const auto *hd_104 = buffer.data(hd + 104);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_109 = buffer.data(hd + 109);
    const auto *hd_110 = buffer.data(hd + 110);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

#pragma omp simd aligned(t_105, t_106, t_107, t_108, pa_x, pb_x, pb_z, gd_65, gf_106, gf_108, \
                         hf_s_105, hf_s_106, hf_s_107, hf_s_108, hd_63, \
                         hd_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_4 * gd_65[k]
                   + f_2 * hf_s_105[k]
                   + pb_x[k] * hd_65[k];

        t_106[k] = pa_x[k] * gf_106[k]
                   + f_2 * hf_s_106[k];

        t_107[k] = f_2 * hf_s_107[k]
                   + pb_z[k] * hd_63[k];

        t_108[k] = pa_x[k] * gf_108[k]
                   + f_2 * hf_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, t_112, pa_x, pa_z, pb_z, gd_36, gf_60, gf_61, \
                         gf_109, hf_s_109, hf_s_110, hf_s_111, hf_s_112, \
                         hd_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = pa_x[k] * gf_109[k]
                   + f_2 * hf_s_109[k];

        t_110[k] = pa_z[k] * gf_60[k]
                   + f_2 * hf_s_110[k];

        t_111[k] = pa_z[k] * gf_61[k]
                   + f_2 * hf_s_111[k];

        t_112[k] = f_4 * gd_36[k]
                   + f_2 * hf_s_112[k]
                   + pb_z[k] * hd_66[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, pa_z, pb_x, gd_70, gd_71, gf_63, hf_s_113, \
                         hf_s_114, hf_s_115, hd_70, hd_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = pa_z[k] * gf_63[k]
                   + f_2 * hf_s_113[k];

        t_114[k] = f_4 * gd_70[k]
                   + f_2 * hf_s_114[k]
                   + pb_x[k] * hd_70[k];

        t_115[k] = f_4 * gd_71[k]
                   + f_2 * hf_s_115[k]
                   + pb_x[k] * hd_71[k];
    }

#pragma omp simd aligned(t_116, t_117, t_118, t_119, pa_x, gf_116, gf_117, gf_118, gf_119, \
                         hf_s_116, hf_s_117, hf_s_118, hf_s_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_116[k] = pa_x[k] * gf_116[k]
                   + f_2 * hf_s_116[k];

        t_117[k] = pa_x[k] * gf_117[k]
                   + f_2 * hf_s_117[k];

        t_118[k] = pa_x[k] * gf_118[k]
                   + f_2 * hf_s_118[k];

        t_119[k] = pa_x[k] * gf_119[k]
                   + f_2 * hf_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_x, pb_y, pb_z, gd_42, gd_48, gd_72, gf_120, \
                         hf_s_120, hf_s_121, hf_s_122, hd_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = f_7 * gd_72[k]
                   + pa_x[k] * gf_120[k]
                   + f_2 * hf_s_120[k];

        t_121[k] = f_3 * gd_48[k]
                   + f_2 * hf_s_121[k]
                   + pb_y[k] * hd_72[k];

        t_122[k] = f_3 * gd_42[k]
                   + f_2 * hf_s_122[k]
                   + pb_z[k] * hd_72[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_x, gd_75, gd_76, gd_77, hf_s_123, hf_s_124, \
                         hf_s_125, hd_75, hd_76, hd_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_4 * gd_75[k]
                   + f_2 * hf_s_123[k]
                   + pb_x[k] * hd_75[k];

        t_124[k] = f_4 * gd_76[k]
                   + f_2 * hf_s_124[k]
                   + pb_x[k] * hd_76[k];

        t_125[k] = f_4 * gd_77[k]
                   + f_2 * hf_s_125[k]
                   + pb_x[k] * hd_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, t_129, pa_x, gf_126, gf_127, gf_128, gf_129, \
                         hf_s_126, hf_s_127, hf_s_128, hf_s_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = pa_x[k] * gf_126[k]
                   + f_2 * hf_s_126[k];

        t_127[k] = pa_x[k] * gf_127[k]
                   + f_2 * hf_s_127[k];

        t_128[k] = pa_x[k] * gf_128[k]
                   + f_2 * hf_s_128[k];

        t_129[k] = pa_x[k] * gf_129[k]
                   + f_2 * hf_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_y, pb_y, gd_54, gf_90, gf_92, hf_s_130, \
                         hf_s_131, hf_s_132, hd_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_y[k] * gf_90[k]
                   + f_2 * hf_s_130[k];

        t_131[k] = f_4 * gd_54[k]
                   + f_2 * hf_s_131[k]
                   + pb_y[k] * hd_78[k];

        t_132[k] = pa_y[k] * gf_92[k]
                   + f_2 * hf_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_y, pb_x, gd_81, gd_82, gf_95, hf_s_133, \
                         hf_s_134, hf_s_135, hd_81, hd_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_4 * gd_81[k]
                   + f_2 * hf_s_133[k]
                   + pb_x[k] * hd_81[k];

        t_134[k] = f_4 * gd_82[k]
                   + f_2 * hf_s_134[k]
                   + pb_x[k] * hd_82[k];

        t_135[k] = pa_y[k] * gf_95[k]
                   + f_2 * hf_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pa_x, gf_136, gf_137, gf_138, gf_139, \
                         hf_s_136, hf_s_137, hf_s_138, hf_s_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = pa_x[k] * gf_136[k]
                   + f_2 * hf_s_136[k];

        t_137[k] = pa_x[k] * gf_137[k]
                   + f_2 * hf_s_137[k];

        t_138[k] = pa_x[k] * gf_138[k]
                   + f_2 * hf_s_138[k];

        t_139[k] = pa_x[k] * gf_139[k]
                   + f_2 * hf_s_139[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pa_x, pb_y, pb_z, gd_54, gd_84, gf_140, \
                         hf_s_140, hf_s_141, hf_s_142, hd_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_7 * gd_84[k]
                   + pa_x[k] * gf_140[k]
                   + f_2 * hf_s_140[k];

        t_141[k] = f_2 * hf_s_141[k]
                   + pb_y[k] * hd_84[k];

        t_142[k] = f_5 * gd_54[k]
                   + f_2 * hf_s_142[k]
                   + pb_z[k] * hd_84[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, pb_y, gd_87, gd_89, hf_s_143, hf_s_144, \
                         hf_s_145, hd_86, hd_87, hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_4 * gd_87[k]
                   + f_2 * hf_s_143[k]
                   + pb_x[k] * hd_87[k];

        t_144[k] = f_2 * hf_s_144[k]
                   + pb_y[k] * hd_86[k];

        t_145[k] = f_4 * gd_89[k]
                   + f_2 * hf_s_145[k]
                   + pb_x[k] * hd_89[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, t_149, pa_x, pb_y, gf_146, gf_147, gf_149, \
                         hf_s_146, hf_s_147, hf_s_148, hf_s_149, \
                         hd_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = pa_x[k] * gf_146[k]
                   + f_2 * hf_s_146[k];

        t_147[k] = pa_x[k] * gf_147[k]
                   + f_2 * hf_s_147[k];

        t_148[k] = f_2 * hf_s_148[k]
                   + pb_y[k] * hd_89[k];

        t_149[k] = pa_x[k] * gf_149[k]
                   + f_2 * hf_s_149[k];
    }

#pragma omp simd aligned(t_150, t_151, t_152, pb_x, pb_z, hp_s_45, hp_s_46, hf_s_150, \
                         hf_s_151, hf_s_152, hp_45, hp_46, hd_90, \
                         hd_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_150[k] = -f_1 * hp_s_45[k]
                   + f_2 * hf_s_150[k]
                   + f_3 * hp_45[k]
                   + pb_x[k] * hd_90[k];

        t_151[k] = -f_8 * hp_s_46[k]
                   + f_2 * hf_s_151[k]
                   + f_4 * hp_46[k]
                   + pb_x[k] * hd_91[k];

        t_152[k] = f_2 * hf_s_152[k]
                   + pb_z[k] * hd_90[k];
    }

#pragma omp simd aligned(t_153, t_154, t_155, t_156, pb_x, pb_y, gd_63, hp_s_46, hf_s_153, \
                         hf_s_154, hf_s_155, hf_s_156, hp_46, hd_93, hd_94, \
                         hd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_153[k] = f_2 * hf_s_153[k]
                   + pb_x[k] * hd_93[k];

        t_154[k] = f_2 * hf_s_154[k]
                   + pb_x[k] * hd_94[k];

        t_155[k] = f_2 * hf_s_155[k]
                   + pb_x[k] * hd_95[k];

        t_156[k] = f_0 * gd_63[k]
                   - f_1 * hp_s_46[k]
                   + f_2 * hf_s_156[k]
                   + f_3 * hp_46[k]
                   + pb_y[k] * hd_93[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pb_y, pb_z, gd_65, hp_s_47, hf_s_157, hf_s_158, \
                         hf_s_159, hp_47, hd_93, hd_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_2 * hf_s_157[k]
                   + pb_z[k] * hd_93[k];

        t_158[k] = f_0 * gd_65[k]
                   + f_2 * hf_s_158[k]
                   + pb_y[k] * hd_95[k];

        t_159[k] = -f_1 * hp_s_47[k]
                   + f_2 * hf_s_159[k]
                   + f_3 * hp_47[k]
                   + pb_z[k] * hd_95[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_z, pb_x, gf_100, gf_101, hp_s_50, \
                         hf_s_160, hf_s_161, hf_s_162, hf_s_163, hp_50, hd_98, \
                         hd_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = pa_z[k] * gf_100[k]
                   + f_2 * hf_s_160[k];

        t_161[k] = pa_z[k] * gf_101[k]
                   + f_2 * hf_s_161[k];

        t_162[k] = -f_8 * hp_s_50[k]
                   + f_2 * hf_s_162[k]
                   + f_4 * hp_50[k]
                   + pb_x[k] * hd_98[k];

        t_163[k] = f_2 * hf_s_163[k]
                   + pb_x[k] * hd_99[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_z, pb_x, pb_z, gd_63, gf_106, \
                         hf_s_164, hf_s_165, hf_s_166, hf_s_167, hd_99, hd_100, \
                         hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = f_2 * hf_s_164[k]
                   + pb_x[k] * hd_100[k];

        t_165[k] = f_2 * hf_s_165[k]
                   + pb_x[k] * hd_101[k];

        t_166[k] = pa_z[k] * gf_106[k]
                   + f_2 * hf_s_166[k];

        t_167[k] = f_4 * gd_63[k]
                   + f_2 * hf_s_167[k]
                   + pb_z[k] * hd_99[k];
    }

#pragma omp simd aligned(t_168, t_169, pa_y, pb_y, ff_s_79, ff_79, gd_71, gf_119, hf_s_168, \
                         hf_s_169, hd_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = f_5 * gd_71[k]
                   + f_2 * hf_s_168[k]
                   + pb_y[k] * hd_101[k];

        t_169[k] = -f_6 * ff_s_79[k]
                   + f_7 * ff_79[k]
                   + pa_y[k] * gf_119[k]
                   + f_2 * hf_s_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pb_x, hp_s_51, hp_s_52, hp_s_53, hf_s_170, \
                         hf_s_171, hf_s_172, hp_51, hp_52, hp_53, hd_102, hd_103, \
                         hd_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -f_1 * hp_s_51[k]
                   + f_2 * hf_s_170[k]
                   + f_3 * hp_51[k]
                   + pb_x[k] * hd_102[k];

        t_171[k] = -f_8 * hp_s_52[k]
                   + f_2 * hf_s_171[k]
                   + f_4 * hp_52[k]
                   + pb_x[k] * hd_103[k];

        t_172[k] = -f_8 * hp_s_53[k]
                   + f_2 * hf_s_172[k]
                   + f_4 * hp_53[k]
                   + pb_x[k] * hd_104[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_z, pb_x, ff_s_66, ff_66, gf_116, \
                         hf_s_173, hf_s_174, hf_s_175, hf_s_176, hd_105, hd_106, \
                         hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_2 * hf_s_173[k]
                   + pb_x[k] * hd_105[k];

        t_174[k] = f_2 * hf_s_174[k]
                   + pb_x[k] * hd_106[k];

        t_175[k] = f_2 * hf_s_175[k]
                   + pb_x[k] * hd_107[k];

        t_176[k] = -f_9 * ff_s_66[k]
                   + f_4 * ff_66[k]
                   + pa_z[k] * gf_116[k]
                   + f_2 * hf_s_176[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pb_y, pb_z, ff_s_89, ff_89, gd_69, gd_77, \
                         gf_129, hf_s_177, hf_s_178, hf_s_179, hd_105, \
                         hd_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_3 * gd_69[k]
                   + f_2 * hf_s_177[k]
                   + pb_z[k] * hd_105[k];

        t_178[k] = f_7 * gd_77[k]
                   + f_2 * hf_s_178[k]
                   + pb_y[k] * hd_107[k];

        t_179[k] = -f_10 * ff_s_89[k]
                   + f_3 * ff_89[k]
                   + pa_y[k] * gf_129[k]
                   + f_2 * hf_s_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_x, hp_s_54, hp_s_55, hp_s_56, hf_s_180, \
                         hf_s_181, hf_s_182, hp_54, hp_55, hp_56, hd_108, hd_109, \
                         hd_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_1 * hp_s_54[k]
                   + f_2 * hf_s_180[k]
                   + f_3 * hp_54[k]
                   + pb_x[k] * hd_108[k];

        t_181[k] = -f_8 * hp_s_55[k]
                   + f_2 * hf_s_181[k]
                   + f_4 * hp_55[k]
                   + pb_x[k] * hd_109[k];

        t_182[k] = -f_8 * hp_s_56[k]
                   + f_2 * hf_s_182[k]
                   + f_4 * hp_56[k]
                   + pb_x[k] * hd_110[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, pa_z, pb_x, ff_s_76, ff_76, gf_126, \
                         hf_s_183, hf_s_184, hf_s_185, hf_s_186, hd_111, hd_112, \
                         hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_2 * hf_s_183[k]
                   + pb_x[k] * hd_111[k];

        t_184[k] = f_2 * hf_s_184[k]
                   + pb_x[k] * hd_112[k];

        t_185[k] = f_2 * hf_s_185[k]
                   + pb_x[k] * hd_113[k];

        t_186[k] = -f_10 * ff_s_76[k]
                   + f_3 * ff_76[k]
                   + pa_z[k] * gf_126[k]
                   + f_2 * hf_s_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pb_y, pb_z, ff_s_99, ff_99, gd_75, gd_83, \
                         gf_139, hf_s_187, hf_s_188, hf_s_189, hd_111, \
                         hd_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_7 * gd_75[k]
                   + f_2 * hf_s_187[k]
                   + pb_z[k] * hd_111[k];

        t_188[k] = f_3 * gd_83[k]
                   + f_2 * hf_s_188[k]
                   + pb_y[k] * hd_113[k];

        t_189[k] = -f_9 * ff_s_99[k]
                   + f_4 * ff_99[k]
                   + pa_y[k] * gf_139[k]
                   + f_2 * hf_s_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pa_y, pb_x, gd_84, gf_140, gf_141, \
                         gf_142, hf_s_190, hf_s_191, hf_s_192, hf_s_193, \
                         hd_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_y[k] * gf_140[k]
                   + f_2 * hf_s_190[k];

        t_191[k] = f_4 * gd_84[k]
                   + pa_y[k] * gf_141[k]
                   + f_2 * hf_s_191[k];

        t_192[k] = pa_y[k] * gf_142[k]
                   + f_2 * hf_s_192[k];

        t_193[k] = f_2 * hf_s_193[k]
                   + pb_x[k] * hd_117[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_y, pb_x, gd_87, gf_146, hf_s_194, hf_s_195, \
                         hf_s_196, hd_118, hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_2 * hf_s_194[k]
                   + pb_x[k] * hd_118[k];

        t_195[k] = f_2 * hf_s_195[k]
                   + pb_x[k] * hd_119[k];

        t_196[k] = f_7 * gd_87[k]
                   + pa_y[k] * gf_146[k]
                   + f_2 * hf_s_196[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, pa_y, pb_y, pb_z, gd_81, gd_89, gf_149, \
                         hf_s_197, hf_s_198, hf_s_199, hd_117, hd_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = f_5 * gd_81[k]
                   + f_2 * hf_s_197[k]
                   + pb_z[k] * hd_117[k];

        t_198[k] = f_4 * gd_89[k]
                   + f_2 * hf_s_198[k]
                   + pb_y[k] * hd_119[k];

        t_199[k] = pa_y[k] * gf_149[k]
                   + f_2 * hf_s_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pb_x, pb_y, hp_s_60, hp_s_62, hf_s_200, \
                         hf_s_201, hf_s_202, hp_60, hp_62, hd_120, \
                         hd_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = -f_1 * hp_s_60[k]
                   + f_2 * hf_s_200[k]
                   + f_3 * hp_60[k]
                   + pb_x[k] * hd_120[k];

        t_201[k] = f_2 * hf_s_201[k]
                   + pb_y[k] * hd_120[k];

        t_202[k] = -f_8 * hp_s_62[k]
                   + f_2 * hf_s_202[k]
                   + f_4 * hp_62[k]
                   + pb_x[k] * hd_122[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, t_206, pb_x, pb_y, hp_s_61, hf_s_203, hf_s_204, \
                         hf_s_205, hf_s_206, hp_61, hd_123, hd_124, \
                         hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_2 * hf_s_203[k]
                   + pb_x[k] * hd_123[k];

        t_204[k] = f_2 * hf_s_204[k]
                   + pb_x[k] * hd_124[k];

        t_205[k] = f_2 * hf_s_205[k]
                   + pb_x[k] * hd_125[k];

        t_206[k] = -f_1 * hp_s_61[k]
                   + f_2 * hf_s_206[k]
                   + f_3 * hp_61[k]
                   + pb_y[k] * hd_123[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pb_y, pb_z, gd_89, hp_s_62, hf_s_207, hf_s_208, \
                         hf_s_209, hp_62, hd_124, hd_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -f_8 * hp_s_62[k]
                   + f_2 * hf_s_207[k]
                   + f_4 * hp_62[k]
                   + pb_y[k] * hd_124[k];

        t_208[k] = f_2 * hf_s_208[k]
                   + pb_y[k] * hd_125[k];

        t_209[k] = f_0 * gd_89[k]
                   - f_1 * hp_s_62[k]
                   + f_2 * hf_s_209[k]
                   + f_3 * hp_62[k]
                   + pb_z[k] * hd_125[k];
    }
}

auto
compute_prim_hf_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t ff_s, const size_t ff,
                                 const size_t gd, const size_t gf, const size_t hp_s,
                                 const size_t hf_s, const size_t hp, const size_t hd,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_hf_kinetic_energy_0_piece0(buffer, target, pa, pb, ff_s, ff, gd, gf, hp_s,
                                            hf_s, hp, hd, ncols, alpha, beta, p);

    compute_prim_hf_kinetic_energy_0_piece1(buffer, target, pa, pb, ff_s, ff, gd, gf, hp_s,
                                            hf_s, hp, hd, ncols, alpha, beta, p);
}

}  // namespace simdkin
