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


#include "SimdKineticEnergyVrrRecHG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_hg_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fg_s,
                                        const size_t fg, const size_t gf, const size_t gg,
                                        const size_t hd_s, const size_t hg_s, const size_t hd,
                                        const size_t hf, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 3.0 * beta / p;
    const auto f_9 = 2.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_s_0 = buffer.data(fg_s + 0);
    const auto *fg_s_15 = buffer.data(fg_s + 15);
    const auto *fg_s_25 = buffer.data(fg_s + 25);
    const auto *fg_s_44 = buffer.data(fg_s + 44);
    const auto *fg_s_55 = buffer.data(fg_s + 55);
    const auto *fg_s_72 = buffer.data(fg_s + 72);
    const auto *fg_s_89 = buffer.data(fg_s + 89);
    const auto *fg_s_100 = buffer.data(fg_s + 100);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_100 = buffer.data(fg + 100);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_100 = buffer.data(gg + 100);

    const auto *hd_s_0 = buffer.data(hd_s + 0);
    const auto *hd_s_3 = buffer.data(hd_s + 3);
    const auto *hd_s_5 = buffer.data(hd_s + 5);
    const auto *hd_s_9 = buffer.data(hd_s + 9);
    const auto *hd_s_16 = buffer.data(hd_s + 16);
    const auto *hd_s_17 = buffer.data(hd_s + 17);
    const auto *hd_s_18 = buffer.data(hd_s + 18);
    const auto *hd_s_21 = buffer.data(hd_s + 21);
    const auto *hd_s_23 = buffer.data(hd_s + 23);
    const auto *hd_s_30 = buffer.data(hd_s + 30);
    const auto *hd_s_33 = buffer.data(hd_s + 33);
    const auto *hd_s_34 = buffer.data(hd_s + 34);
    const auto *hd_s_35 = buffer.data(hd_s + 35);
    const auto *hd_s_36 = buffer.data(hd_s + 36);
    const auto *hd_s_39 = buffer.data(hd_s + 39);

    const auto *hg_s_0 = buffer.data(hg_s + 0);
    const auto *hg_s_1 = buffer.data(hg_s + 1);
    const auto *hg_s_2 = buffer.data(hg_s + 2);
    const auto *hg_s_3 = buffer.data(hg_s + 3);
    const auto *hg_s_4 = buffer.data(hg_s + 4);
    const auto *hg_s_5 = buffer.data(hg_s + 5);
    const auto *hg_s_6 = buffer.data(hg_s + 6);
    const auto *hg_s_7 = buffer.data(hg_s + 7);
    const auto *hg_s_8 = buffer.data(hg_s + 8);
    const auto *hg_s_9 = buffer.data(hg_s + 9);
    const auto *hg_s_10 = buffer.data(hg_s + 10);
    const auto *hg_s_11 = buffer.data(hg_s + 11);
    const auto *hg_s_12 = buffer.data(hg_s + 12);
    const auto *hg_s_13 = buffer.data(hg_s + 13);
    const auto *hg_s_14 = buffer.data(hg_s + 14);
    const auto *hg_s_15 = buffer.data(hg_s + 15);
    const auto *hg_s_16 = buffer.data(hg_s + 16);
    const auto *hg_s_17 = buffer.data(hg_s + 17);
    const auto *hg_s_18 = buffer.data(hg_s + 18);
    const auto *hg_s_19 = buffer.data(hg_s + 19);
    const auto *hg_s_20 = buffer.data(hg_s + 20);
    const auto *hg_s_21 = buffer.data(hg_s + 21);
    const auto *hg_s_22 = buffer.data(hg_s + 22);
    const auto *hg_s_23 = buffer.data(hg_s + 23);
    const auto *hg_s_24 = buffer.data(hg_s + 24);
    const auto *hg_s_25 = buffer.data(hg_s + 25);
    const auto *hg_s_26 = buffer.data(hg_s + 26);
    const auto *hg_s_27 = buffer.data(hg_s + 27);
    const auto *hg_s_28 = buffer.data(hg_s + 28);
    const auto *hg_s_29 = buffer.data(hg_s + 29);
    const auto *hg_s_30 = buffer.data(hg_s + 30);
    const auto *hg_s_31 = buffer.data(hg_s + 31);
    const auto *hg_s_32 = buffer.data(hg_s + 32);
    const auto *hg_s_33 = buffer.data(hg_s + 33);
    const auto *hg_s_34 = buffer.data(hg_s + 34);
    const auto *hg_s_35 = buffer.data(hg_s + 35);
    const auto *hg_s_36 = buffer.data(hg_s + 36);
    const auto *hg_s_37 = buffer.data(hg_s + 37);
    const auto *hg_s_38 = buffer.data(hg_s + 38);
    const auto *hg_s_39 = buffer.data(hg_s + 39);
    const auto *hg_s_40 = buffer.data(hg_s + 40);
    const auto *hg_s_41 = buffer.data(hg_s + 41);
    const auto *hg_s_42 = buffer.data(hg_s + 42);
    const auto *hg_s_43 = buffer.data(hg_s + 43);
    const auto *hg_s_44 = buffer.data(hg_s + 44);
    const auto *hg_s_45 = buffer.data(hg_s + 45);
    const auto *hg_s_46 = buffer.data(hg_s + 46);
    const auto *hg_s_47 = buffer.data(hg_s + 47);
    const auto *hg_s_48 = buffer.data(hg_s + 48);
    const auto *hg_s_49 = buffer.data(hg_s + 49);
    const auto *hg_s_50 = buffer.data(hg_s + 50);
    const auto *hg_s_51 = buffer.data(hg_s + 51);
    const auto *hg_s_52 = buffer.data(hg_s + 52);
    const auto *hg_s_53 = buffer.data(hg_s + 53);
    const auto *hg_s_54 = buffer.data(hg_s + 54);
    const auto *hg_s_55 = buffer.data(hg_s + 55);
    const auto *hg_s_56 = buffer.data(hg_s + 56);
    const auto *hg_s_57 = buffer.data(hg_s + 57);
    const auto *hg_s_58 = buffer.data(hg_s + 58);
    const auto *hg_s_59 = buffer.data(hg_s + 59);
    const auto *hg_s_60 = buffer.data(hg_s + 60);
    const auto *hg_s_61 = buffer.data(hg_s + 61);
    const auto *hg_s_62 = buffer.data(hg_s + 62);
    const auto *hg_s_63 = buffer.data(hg_s + 63);
    const auto *hg_s_64 = buffer.data(hg_s + 64);
    const auto *hg_s_65 = buffer.data(hg_s + 65);
    const auto *hg_s_66 = buffer.data(hg_s + 66);
    const auto *hg_s_67 = buffer.data(hg_s + 67);
    const auto *hg_s_68 = buffer.data(hg_s + 68);
    const auto *hg_s_69 = buffer.data(hg_s + 69);
    const auto *hg_s_70 = buffer.data(hg_s + 70);
    const auto *hg_s_71 = buffer.data(hg_s + 71);
    const auto *hg_s_72 = buffer.data(hg_s + 72);
    const auto *hg_s_73 = buffer.data(hg_s + 73);
    const auto *hg_s_74 = buffer.data(hg_s + 74);
    const auto *hg_s_75 = buffer.data(hg_s + 75);
    const auto *hg_s_76 = buffer.data(hg_s + 76);
    const auto *hg_s_77 = buffer.data(hg_s + 77);
    const auto *hg_s_78 = buffer.data(hg_s + 78);
    const auto *hg_s_79 = buffer.data(hg_s + 79);
    const auto *hg_s_80 = buffer.data(hg_s + 80);
    const auto *hg_s_81 = buffer.data(hg_s + 81);
    const auto *hg_s_82 = buffer.data(hg_s + 82);
    const auto *hg_s_83 = buffer.data(hg_s + 83);
    const auto *hg_s_84 = buffer.data(hg_s + 84);
    const auto *hg_s_85 = buffer.data(hg_s + 85);
    const auto *hg_s_86 = buffer.data(hg_s + 86);
    const auto *hg_s_87 = buffer.data(hg_s + 87);
    const auto *hg_s_88 = buffer.data(hg_s + 88);
    const auto *hg_s_89 = buffer.data(hg_s + 89);
    const auto *hg_s_90 = buffer.data(hg_s + 90);
    const auto *hg_s_91 = buffer.data(hg_s + 91);
    const auto *hg_s_92 = buffer.data(hg_s + 92);
    const auto *hg_s_93 = buffer.data(hg_s + 93);
    const auto *hg_s_94 = buffer.data(hg_s + 94);
    const auto *hg_s_95 = buffer.data(hg_s + 95);
    const auto *hg_s_96 = buffer.data(hg_s + 96);
    const auto *hg_s_97 = buffer.data(hg_s + 97);
    const auto *hg_s_98 = buffer.data(hg_s + 98);
    const auto *hg_s_99 = buffer.data(hg_s + 99);
    const auto *hg_s_100 = buffer.data(hg_s + 100);
    const auto *hg_s_101 = buffer.data(hg_s + 101);
    const auto *hg_s_102 = buffer.data(hg_s + 102);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_16 = buffer.data(hd + 16);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_34 = buffer.data(hd + 34);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_39 = buffer.data(hd + 39);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gf_0, hd_s_0, hg_s_0, hg_s_1, \
                         hg_s_2, hg_s_3, hd_0, hf_0, hf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gf_0[k]
                 - f_1 * hd_s_0[k]
                 + f_2 * hg_s_0[k]
                 + f_3 * hd_0[k]
                 + pb_x[k] * hf_0[k];

        t_1[k] = f_2 * hg_s_1[k]
                 + pb_y[k] * hf_0[k];

        t_2[k] = f_2 * hg_s_2[k]
                 + pb_z[k] * hf_0[k];

        t_3[k] = -f_4 * hd_s_0[k]
                 + f_2 * hg_s_3[k]
                 + f_5 * hd_0[k]
                 + pb_y[k] * hf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, gf_6, hd_s_0, hg_s_4, hg_s_5, \
                         hg_s_6, hd_0, hf_2, hf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hg_s_4[k]
                 + pb_y[k] * hf_2[k];

        t_5[k] = -f_4 * hd_s_0[k]
                 + f_2 * hg_s_5[k]
                 + f_5 * hd_0[k]
                 + pb_z[k] * hf_2[k];

        t_6[k] = f_0 * gf_6[k]
                 + f_2 * hg_s_6[k]
                 + pb_x[k] * hf_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, gf_9, hg_s_7, hg_s_8, hg_s_9, hf_3, \
                         hf_5, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * hg_s_7[k]
                 + pb_z[k] * hf_3[k];

        t_8[k] = f_2 * hg_s_8[k]
                 + pb_y[k] * hf_5[k];

        t_9[k] = f_0 * gf_9[k]
                 + f_2 * hg_s_9[k]
                 + pb_x[k] * hf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, hd_s_3, hd_s_5, hg_s_10, hg_s_11, \
                         hg_s_12, hd_3, hd_5, hf_6, hf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * hd_s_3[k]
                  + f_2 * hg_s_10[k]
                  + f_3 * hd_3[k]
                  + pb_y[k] * hf_6[k];

        t_11[k] = f_2 * hg_s_11[k]
                  + pb_z[k] * hf_6[k];

        t_12[k] = -f_4 * hd_s_5[k]
                  + f_2 * hg_s_12[k]
                  + f_5 * hd_5[k]
                  + pb_y[k] * hf_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, gg_0, hd_s_5, hg_s_13, hg_s_14, \
                         hg_s_15, hd_5, hf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * hg_s_13[k]
                  + pb_y[k] * hf_9[k];

        t_14[k] = -f_1 * hd_s_5[k]
                  + f_2 * hg_s_14[k]
                  + f_3 * hd_5[k]
                  + pb_z[k] * hf_9[k];

        t_15[k] = pa_y[k] * gg_0[k]
                  + f_2 * hg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, gf_0, gf_1, gg_3, hg_s_16, \
                         hg_s_17, hg_s_18, hg_s_19, hf_10, hf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * gf_0[k]
                  + f_2 * hg_s_16[k]
                  + pb_y[k] * hf_10[k];

        t_17[k] = f_2 * hg_s_17[k]
                  + pb_z[k] * hf_10[k];

        t_18[k] = f_6 * gf_1[k]
                  + pa_y[k] * gg_3[k]
                  + f_2 * hg_s_18[k];

        t_19[k] = f_2 * hg_s_19[k]
                  + pb_z[k] * hf_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, gf_16, gg_5, hg_s_20, hg_s_21, \
                         hg_s_22, hf_13, hf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * gg_5[k]
                  + f_2 * hg_s_20[k];

        t_21[k] = f_7 * gf_16[k]
                  + f_2 * hg_s_21[k]
                  + pb_x[k] * hf_16[k];

        t_22[k] = f_2 * hg_s_22[k]
                  + pb_z[k] * hf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, fg_s_25, fg_25, gf_18, gg_9, \
                         gg_25, hg_s_23, hg_s_24, hg_s_25, hf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * gf_18[k]
                  + f_2 * hg_s_23[k]
                  + pb_x[k] * hf_18[k];

        t_24[k] = pa_y[k] * gg_9[k]
                  + f_2 * hg_s_24[k];

        t_25[k] = -f_8 * fg_s_25[k]
                  + f_3 * fg_25[k]
                  + pa_x[k] * gg_25[k]
                  + f_2 * hg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, gf_9, hd_s_9, hg_s_26, hg_s_27, \
                         hg_s_28, hd_9, hf_16, hf_17, hf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * hg_s_26[k]
                  + pb_z[k] * hf_16[k];

        t_27[k] = -f_4 * hd_s_9[k]
                  + f_2 * hg_s_27[k]
                  + f_5 * hd_9[k]
                  + pb_z[k] * hf_17[k];

        t_28[k] = f_5 * gf_9[k]
                  + f_2 * hg_s_28[k]
                  + pb_y[k] * hf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, gf_0, gg_0, gg_14, \
                         hg_s_29, hg_s_30, hg_s_31, hg_s_32, hf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * gg_14[k]
                  + f_2 * hg_s_29[k];

        t_30[k] = pa_z[k] * gg_0[k]
                  + f_2 * hg_s_30[k];

        t_31[k] = f_2 * hg_s_31[k]
                  + pb_y[k] * hf_20[k];

        t_32[k] = f_5 * gf_0[k]
                  + f_2 * hg_s_32[k]
                  + pb_z[k] * hf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, gf_2, gg_3, gg_5, gg_6, hg_s_33, \
                         hg_s_34, hg_s_35, hg_s_36, hf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * gg_3[k]
                  + f_2 * hg_s_33[k];

        t_34[k] = f_2 * hg_s_34[k]
                  + pb_y[k] * hf_22[k];

        t_35[k] = f_6 * gf_2[k]
                  + pa_z[k] * gg_5[k]
                  + f_2 * hg_s_35[k];

        t_36[k] = pa_z[k] * gg_6[k]
                  + f_2 * hg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, gf_27, gf_29, hg_s_37, hg_s_38, \
                         hg_s_39, hf_25, hf_27, hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * gf_27[k]
                  + f_2 * hg_s_37[k]
                  + pb_x[k] * hf_27[k];

        t_38[k] = f_2 * hg_s_38[k]
                  + pb_y[k] * hf_25[k];

        t_39[k] = f_7 * gf_29[k]
                  + f_2 * hg_s_39[k]
                  + pb_x[k] * hf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, gg_10, hd_s_16, hd_s_17, hg_s_40, \
                         hg_s_41, hg_s_42, hd_16, hd_17, hf_27, hf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * gg_10[k]
                  + f_2 * hg_s_40[k];

        t_41[k] = -f_9 * hd_s_16[k]
                  + f_2 * hg_s_41[k]
                  + f_6 * hd_16[k]
                  + pb_y[k] * hf_27[k];

        t_42[k] = -f_4 * hd_s_17[k]
                  + f_2 * hg_s_42[k]
                  + f_5 * hd_17[k]
                  + pb_y[k] * hf_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pb_y, fg_s_0, fg_s_44, fg_0, fg_44, \
                         gg_15, gg_44, hg_s_43, hg_s_44, hg_s_45, \
                         hf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * hg_s_43[k]
                  + pb_y[k] * hf_29[k];

        t_44[k] = -f_8 * fg_s_44[k]
                  + f_3 * fg_44[k]
                  + pa_x[k] * gg_44[k]
                  + f_2 * hg_s_44[k];

        t_45[k] = -f_10 * fg_s_0[k]
                  + f_5 * fg_0[k]
                  + pa_y[k] * gg_15[k]
                  + f_2 * hg_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, pb_z, gf_10, gf_33, hd_s_21, hg_s_46, \
                         hg_s_47, hg_s_48, hd_21, hf_30, hf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * gf_10[k]
                  + f_2 * hg_s_46[k]
                  + pb_y[k] * hf_30[k];

        t_47[k] = f_2 * hg_s_47[k]
                  + pb_z[k] * hf_30[k];

        t_48[k] = f_3 * gf_33[k]
                  - f_4 * hd_s_21[k]
                  + f_2 * hg_s_48[k]
                  + f_5 * hd_21[k]
                  + pb_x[k] * hf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_z, gf_36, hd_s_18, hg_s_49, hg_s_50, \
                         hg_s_51, hd_18, hf_31, hf_32, hf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * hg_s_49[k]
                  + pb_z[k] * hf_31[k];

        t_50[k] = -f_4 * hd_s_18[k]
                  + f_2 * hg_s_50[k]
                  + f_5 * hd_18[k]
                  + pb_z[k] * hf_32[k];

        t_51[k] = f_3 * gf_36[k]
                  + f_2 * hg_s_51[k]
                  + pb_x[k] * hf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_z, gf_38, gf_39, hg_s_52, hg_s_53, \
                         hg_s_54, hf_33, hf_38, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * hg_s_52[k]
                  + pb_z[k] * hf_33[k];

        t_53[k] = f_3 * gf_38[k]
                  + f_2 * hg_s_53[k]
                  + pb_x[k] * hf_38[k];

        t_54[k] = f_3 * gf_39[k]
                  + f_2 * hg_s_54[k]
                  + pb_x[k] * hf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_x, pb_z, fg_s_55, fg_55, gg_55, hd_s_21, \
                         hg_s_55, hg_s_56, hg_s_57, hd_21, hf_36, \
                         hf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_11 * fg_s_55[k]
                  + f_6 * fg_55[k]
                  + pa_x[k] * gg_55[k]
                  + f_2 * hg_s_55[k];

        t_56[k] = f_2 * hg_s_56[k]
                  + pb_z[k] * hf_36[k];

        t_57[k] = -f_4 * hd_s_21[k]
                  + f_2 * hg_s_57[k]
                  + f_5 * hd_21[k]
                  + pb_z[k] * hf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pb_y, pb_z, gf_19, gg_30, hd_s_23, hg_s_58, \
                         hg_s_59, hg_s_60, hd_23, hf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_6 * gf_19[k]
                  + f_2 * hg_s_58[k]
                  + pb_y[k] * hf_39[k];

        t_59[k] = -f_1 * hd_s_23[k]
                  + f_2 * hg_s_59[k]
                  + f_3 * hd_23[k]
                  + pb_z[k] * hf_39[k];

        t_60[k] = pa_y[k] * gg_30[k]
                  + f_2 * hg_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pa_z, pb_y, gf_22, gg_16, gg_18, gg_32, \
                         hg_s_61, hg_s_62, hg_s_63, hg_s_64, hf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_z[k] * gg_16[k]
                  + f_2 * hg_s_61[k];

        t_62[k] = pa_y[k] * gg_32[k]
                  + f_2 * hg_s_62[k];

        t_63[k] = pa_z[k] * gg_18[k]
                  + f_2 * hg_s_63[k];

        t_64[k] = f_5 * gf_22[k]
                  + f_2 * hg_s_64[k]
                  + pb_y[k] * hf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pa_z, pb_x, gf_47, gg_21, gg_35, hg_s_65, \
                         hg_s_66, hg_s_67, hf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * gg_35[k]
                  + f_2 * hg_s_65[k];

        t_66[k] = pa_z[k] * gg_21[k]
                  + f_2 * hg_s_66[k];

        t_67[k] = f_3 * gf_47[k]
                  + f_2 * hg_s_67[k]
                  + pb_x[k] * hf_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pa_z, pb_x, gf_48, gg_25, gg_39, hg_s_68, \
                         hg_s_69, hg_s_70, hf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_3 * gf_48[k]
                  + f_2 * hg_s_68[k]
                  + pb_x[k] * hf_48[k];

        t_69[k] = pa_y[k] * gg_39[k]
                  + f_2 * hg_s_69[k];

        t_70[k] = pa_z[k] * gg_25[k]
                  + f_2 * hg_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_x, pb_y, pb_z, fg_s_72, fg_72, gf_16, gf_29, \
                         gg_72, hg_s_71, hg_s_72, hg_s_73, hf_46, \
                         hf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * gf_16[k]
                  + f_2 * hg_s_71[k]
                  + pb_z[k] * hf_46[k];

        t_72[k] = -f_11 * fg_s_72[k]
                  + f_6 * fg_72[k]
                  + pa_x[k] * gg_72[k]
                  + f_2 * hg_s_72[k];

        t_73[k] = f_5 * gf_29[k]
                  + f_2 * hg_s_73[k]
                  + pb_y[k] * hf_49[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_y, fg_s_0, fg_0, gg_30, gg_44, \
                         hg_s_74, hg_s_75, hg_s_76, hf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_y[k] * gg_44[k]
                  + f_2 * hg_s_74[k];

        t_75[k] = -f_10 * fg_s_0[k]
                  + f_5 * fg_0[k]
                  + pa_z[k] * gg_30[k]
                  + f_2 * hg_s_75[k];

        t_76[k] = f_2 * hg_s_76[k]
                  + pb_y[k] * hf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pb_z, gf_20, hd_s_30, hg_s_77, hg_s_78, \
                         hg_s_79, hd_30, hf_50, hf_51, hf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_6 * gf_20[k]
                  + f_2 * hg_s_77[k]
                  + pb_z[k] * hf_50[k];

        t_78[k] = -f_4 * hd_s_30[k]
                  + f_2 * hg_s_78[k]
                  + f_5 * hd_30[k]
                  + pb_y[k] * hf_51[k];

        t_79[k] = f_2 * hg_s_79[k]
                  + pb_y[k] * hf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, gf_55, gf_56, gf_57, hd_s_35, hg_s_80, \
                         hg_s_81, hg_s_82, hd_35, hf_55, hf_56, hf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_3 * gf_55[k]
                  - f_4 * hd_s_35[k]
                  + f_2 * hg_s_80[k]
                  + f_5 * hd_35[k]
                  + pb_x[k] * hf_55[k];

        t_81[k] = f_3 * gf_56[k]
                  + f_2 * hg_s_81[k]
                  + pb_x[k] * hf_56[k];

        t_82[k] = f_3 * gf_57[k]
                  + f_2 * hg_s_82[k]
                  + pb_x[k] * hf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, pb_y, gf_59, hd_s_33, hg_s_83, hg_s_84, \
                         hg_s_85, hd_33, hf_55, hf_56, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * hg_s_83[k]
                  + pb_y[k] * hf_55[k];

        t_84[k] = f_3 * gf_59[k]
                  + f_2 * hg_s_84[k]
                  + pb_x[k] * hf_59[k];

        t_85[k] = -f_1 * hd_s_33[k]
                  + f_2 * hg_s_85[k]
                  + f_3 * hd_33[k]
                  + pb_y[k] * hf_56[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, hd_s_34, hd_s_35, hg_s_86, hg_s_87, hg_s_88, \
                         hd_34, hd_35, hf_57, hf_58, hf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_9 * hd_s_34[k]
                  + f_2 * hg_s_86[k]
                  + f_6 * hd_34[k]
                  + pb_y[k] * hf_57[k];

        t_87[k] = -f_4 * hd_s_35[k]
                  + f_2 * hg_s_87[k]
                  + f_5 * hd_35[k]
                  + pb_y[k] * hf_58[k];

        t_88[k] = f_2 * hg_s_88[k]
                  + pb_y[k] * hf_59[k];
    }

#pragma omp simd aligned(t_89, t_90, pa_x, pa_y, fg_s_15, fg_s_89, fg_15, fg_89, gg_45, gg_89, \
                         hg_s_89, hg_s_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -f_11 * fg_s_89[k]
                  + f_6 * fg_89[k]
                  + pa_x[k] * gg_89[k]
                  + f_2 * hg_s_89[k];

        t_90[k] = -f_11 * fg_s_15[k]
                  + f_6 * fg_15[k]
                  + pa_y[k] * gg_45[k]
                  + f_2 * hg_s_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_y, pb_z, gf_30, gf_63, hd_s_39, hg_s_91, \
                         hg_s_92, hg_s_93, hd_39, hf_60, hf_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * gf_30[k]
                  + f_2 * hg_s_91[k]
                  + pb_y[k] * hf_60[k];

        t_92[k] = f_2 * hg_s_92[k]
                  + pb_z[k] * hf_60[k];

        t_93[k] = f_6 * gf_63[k]
                  - f_4 * hd_s_39[k]
                  + f_2 * hg_s_93[k]
                  + f_5 * hd_39[k]
                  + pb_x[k] * hf_63[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, gf_66, hd_s_36, hg_s_94, hg_s_95, \
                         hg_s_96, hd_36, hf_61, hf_62, hf_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * hg_s_94[k]
                  + pb_z[k] * hf_61[k];

        t_95[k] = -f_4 * hd_s_36[k]
                  + f_2 * hg_s_95[k]
                  + f_5 * hd_36[k]
                  + pb_z[k] * hf_62[k];

        t_96[k] = f_6 * gf_66[k]
                  + f_2 * hg_s_96[k]
                  + pb_x[k] * hf_66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, gf_68, gf_69, hg_s_97, hg_s_98, \
                         hg_s_99, hf_63, hf_68, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * hg_s_97[k]
                  + pb_z[k] * hf_63[k];

        t_98[k] = f_6 * gf_68[k]
                  + f_2 * hg_s_98[k]
                  + pb_x[k] * hf_68[k];

        t_99[k] = f_6 * gf_69[k]
                  + f_2 * hg_s_99[k]
                  + pb_x[k] * hf_69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_x, pb_z, fg_s_100, fg_100, gg_100, hd_s_39, \
                         hg_s_100, hg_s_101, hg_s_102, hd_39, hf_66, \
                         hf_67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -f_10 * fg_s_100[k]
                   + f_5 * fg_100[k]
                   + pa_x[k] * gg_100[k]
                   + f_2 * hg_s_100[k];

        t_101[k] = f_2 * hg_s_101[k]
                   + pb_z[k] * hf_66[k];

        t_102[k] = -f_4 * hd_s_39[k]
                   + f_2 * hg_s_102[k]
                   + f_5 * hd_39[k]
                   + pb_z[k] * hf_67[k];
    }
}

static auto
compute_prim_hg_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fg_s,
                                        const size_t fg, const size_t gf, const size_t gg,
                                        const size_t hd_s, const size_t hg_s, const size_t hd,
                                        const size_t hf, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_9 = 2.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_s_30 = buffer.data(fg_s + 30);
    const auto *fg_s_35 = buffer.data(fg_s + 35);
    const auto *fg_s_117 = buffer.data(fg_s + 117);
    const auto *fg_s_119 = buffer.data(fg_s + 119);
    const auto *fg_s_130 = buffer.data(fg_s + 130);
    const auto *fg_s_132 = buffer.data(fg_s + 132);
    const auto *fg_s_149 = buffer.data(fg_s + 149);

    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_70 = buffer.data(gf + 70);
    const auto *gf_72 = buffer.data(gf + 72);
    const auto *gf_77 = buffer.data(gf + 77);
    const auto *gf_78 = buffer.data(gf + 78);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_80 = buffer.data(gf + 80);
    const auto *gf_82 = buffer.data(gf + 82);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_87 = buffer.data(gf + 87);
    const auto *gf_88 = buffer.data(gf + 88);
    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_92 = buffer.data(gf + 92);
    const auto *gf_95 = buffer.data(gf + 95);
    const auto *gf_96 = buffer.data(gf + 96);
    const auto *gf_97 = buffer.data(gf + 97);
    const auto *gf_99 = buffer.data(gf + 99);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_133 = buffer.data(gf + 133);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_140 = buffer.data(gf + 140);

    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);

    const auto *hd_s_41 = buffer.data(hd_s + 41);
    const auto *hd_s_54 = buffer.data(hd_s + 54);
    const auto *hd_s_57 = buffer.data(hd_s + 57);
    const auto *hd_s_58 = buffer.data(hd_s + 58);
    const auto *hd_s_59 = buffer.data(hd_s + 59);

    const auto *hg_s_103 = buffer.data(hg_s + 103);
    const auto *hg_s_104 = buffer.data(hg_s + 104);
    const auto *hg_s_105 = buffer.data(hg_s + 105);
    const auto *hg_s_106 = buffer.data(hg_s + 106);
    const auto *hg_s_107 = buffer.data(hg_s + 107);
    const auto *hg_s_108 = buffer.data(hg_s + 108);
    const auto *hg_s_109 = buffer.data(hg_s + 109);
    const auto *hg_s_110 = buffer.data(hg_s + 110);
    const auto *hg_s_111 = buffer.data(hg_s + 111);
    const auto *hg_s_112 = buffer.data(hg_s + 112);
    const auto *hg_s_113 = buffer.data(hg_s + 113);
    const auto *hg_s_114 = buffer.data(hg_s + 114);
    const auto *hg_s_115 = buffer.data(hg_s + 115);
    const auto *hg_s_116 = buffer.data(hg_s + 116);
    const auto *hg_s_117 = buffer.data(hg_s + 117);
    const auto *hg_s_118 = buffer.data(hg_s + 118);
    const auto *hg_s_119 = buffer.data(hg_s + 119);
    const auto *hg_s_120 = buffer.data(hg_s + 120);
    const auto *hg_s_121 = buffer.data(hg_s + 121);
    const auto *hg_s_122 = buffer.data(hg_s + 122);
    const auto *hg_s_123 = buffer.data(hg_s + 123);
    const auto *hg_s_124 = buffer.data(hg_s + 124);
    const auto *hg_s_125 = buffer.data(hg_s + 125);
    const auto *hg_s_126 = buffer.data(hg_s + 126);
    const auto *hg_s_127 = buffer.data(hg_s + 127);
    const auto *hg_s_128 = buffer.data(hg_s + 128);
    const auto *hg_s_129 = buffer.data(hg_s + 129);
    const auto *hg_s_130 = buffer.data(hg_s + 130);
    const auto *hg_s_131 = buffer.data(hg_s + 131);
    const auto *hg_s_132 = buffer.data(hg_s + 132);
    const auto *hg_s_133 = buffer.data(hg_s + 133);
    const auto *hg_s_134 = buffer.data(hg_s + 134);
    const auto *hg_s_135 = buffer.data(hg_s + 135);
    const auto *hg_s_136 = buffer.data(hg_s + 136);
    const auto *hg_s_137 = buffer.data(hg_s + 137);
    const auto *hg_s_138 = buffer.data(hg_s + 138);
    const auto *hg_s_139 = buffer.data(hg_s + 139);
    const auto *hg_s_140 = buffer.data(hg_s + 140);
    const auto *hg_s_141 = buffer.data(hg_s + 141);
    const auto *hg_s_142 = buffer.data(hg_s + 142);
    const auto *hg_s_143 = buffer.data(hg_s + 143);
    const auto *hg_s_144 = buffer.data(hg_s + 144);
    const auto *hg_s_145 = buffer.data(hg_s + 145);
    const auto *hg_s_146 = buffer.data(hg_s + 146);
    const auto *hg_s_147 = buffer.data(hg_s + 147);
    const auto *hg_s_148 = buffer.data(hg_s + 148);
    const auto *hg_s_149 = buffer.data(hg_s + 149);
    const auto *hg_s_150 = buffer.data(hg_s + 150);
    const auto *hg_s_151 = buffer.data(hg_s + 151);
    const auto *hg_s_152 = buffer.data(hg_s + 152);
    const auto *hg_s_153 = buffer.data(hg_s + 153);
    const auto *hg_s_154 = buffer.data(hg_s + 154);
    const auto *hg_s_155 = buffer.data(hg_s + 155);
    const auto *hg_s_156 = buffer.data(hg_s + 156);
    const auto *hg_s_157 = buffer.data(hg_s + 157);
    const auto *hg_s_158 = buffer.data(hg_s + 158);
    const auto *hg_s_159 = buffer.data(hg_s + 159);
    const auto *hg_s_160 = buffer.data(hg_s + 160);
    const auto *hg_s_161 = buffer.data(hg_s + 161);
    const auto *hg_s_162 = buffer.data(hg_s + 162);
    const auto *hg_s_163 = buffer.data(hg_s + 163);
    const auto *hg_s_164 = buffer.data(hg_s + 164);
    const auto *hg_s_165 = buffer.data(hg_s + 165);
    const auto *hg_s_166 = buffer.data(hg_s + 166);
    const auto *hg_s_167 = buffer.data(hg_s + 167);
    const auto *hg_s_168 = buffer.data(hg_s + 168);
    const auto *hg_s_169 = buffer.data(hg_s + 169);
    const auto *hg_s_170 = buffer.data(hg_s + 170);
    const auto *hg_s_171 = buffer.data(hg_s + 171);
    const auto *hg_s_172 = buffer.data(hg_s + 172);
    const auto *hg_s_173 = buffer.data(hg_s + 173);
    const auto *hg_s_174 = buffer.data(hg_s + 174);
    const auto *hg_s_175 = buffer.data(hg_s + 175);
    const auto *hg_s_176 = buffer.data(hg_s + 176);
    const auto *hg_s_177 = buffer.data(hg_s + 177);
    const auto *hg_s_178 = buffer.data(hg_s + 178);
    const auto *hg_s_179 = buffer.data(hg_s + 179);
    const auto *hg_s_180 = buffer.data(hg_s + 180);
    const auto *hg_s_181 = buffer.data(hg_s + 181);
    const auto *hg_s_182 = buffer.data(hg_s + 182);
    const auto *hg_s_183 = buffer.data(hg_s + 183);
    const auto *hg_s_184 = buffer.data(hg_s + 184);
    const auto *hg_s_185 = buffer.data(hg_s + 185);
    const auto *hg_s_186 = buffer.data(hg_s + 186);
    const auto *hg_s_187 = buffer.data(hg_s + 187);
    const auto *hg_s_188 = buffer.data(hg_s + 188);
    const auto *hg_s_189 = buffer.data(hg_s + 189);
    const auto *hg_s_190 = buffer.data(hg_s + 190);
    const auto *hg_s_191 = buffer.data(hg_s + 191);
    const auto *hg_s_192 = buffer.data(hg_s + 192);
    const auto *hg_s_193 = buffer.data(hg_s + 193);
    const auto *hg_s_194 = buffer.data(hg_s + 194);
    const auto *hg_s_195 = buffer.data(hg_s + 195);
    const auto *hg_s_196 = buffer.data(hg_s + 196);
    const auto *hg_s_197 = buffer.data(hg_s + 197);
    const auto *hg_s_198 = buffer.data(hg_s + 198);
    const auto *hg_s_199 = buffer.data(hg_s + 199);
    const auto *hg_s_200 = buffer.data(hg_s + 200);
    const auto *hg_s_201 = buffer.data(hg_s + 201);
    const auto *hg_s_202 = buffer.data(hg_s + 202);
    const auto *hg_s_203 = buffer.data(hg_s + 203);
    const auto *hg_s_204 = buffer.data(hg_s + 204);
    const auto *hg_s_205 = buffer.data(hg_s + 205);
    const auto *hg_s_206 = buffer.data(hg_s + 206);
    const auto *hg_s_207 = buffer.data(hg_s + 207);
    const auto *hg_s_208 = buffer.data(hg_s + 208);
    const auto *hg_s_209 = buffer.data(hg_s + 209);
    const auto *hg_s_210 = buffer.data(hg_s + 210);
    const auto *hg_s_211 = buffer.data(hg_s + 211);

    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_58 = buffer.data(hd + 58);
    const auto *hd_59 = buffer.data(hd + 59);

    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_70 = buffer.data(hf + 70);
    const auto *hf_72 = buffer.data(hf + 72);
    const auto *hf_76 = buffer.data(hf + 76);
    const auto *hf_77 = buffer.data(hf + 77);
    const auto *hf_78 = buffer.data(hf + 78);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_82 = buffer.data(hf + 82);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_87 = buffer.data(hf + 87);
    const auto *hf_88 = buffer.data(hf + 88);
    const auto *hf_89 = buffer.data(hf + 89);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);
    const auto *hf_140 = buffer.data(hf + 140);

#pragma omp simd aligned(t_103, t_104, t_105, pa_z, pb_y, pb_z, gf_39, gg_45, hd_s_41, \
                         hg_s_103, hg_s_104, hg_s_105, hd_41, hf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * gf_39[k]
                   + f_2 * hg_s_103[k]
                   + pb_y[k] * hf_69[k];

        t_104[k] = -f_1 * hd_s_41[k]
                   + f_2 * hg_s_104[k]
                   + f_3 * hd_41[k]
                   + pb_z[k] * hf_69[k];

        t_105[k] = pa_z[k] * gg_45[k]
                   + f_2 * hg_s_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pa_z, pb_z, gf_30, gg_46, gg_48, hg_s_106, \
                         hg_s_107, hg_s_108, hf_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_z[k] * gg_46[k]
                   + f_2 * hg_s_106[k];

        t_107[k] = f_5 * gf_30[k]
                   + f_2 * hg_s_107[k]
                   + pb_z[k] * hf_70[k];

        t_108[k] = pa_z[k] * gg_48[k]
                   + f_2 * hg_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pa_y, pa_z, pb_y, fg_s_35, fg_35, gf_42, gg_51, \
                         gg_65, hg_s_109, hg_s_110, hg_s_111, hf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * gf_42[k]
                   + f_2 * hg_s_109[k]
                   + pb_y[k] * hf_72[k];

        t_110[k] = -f_10 * fg_s_35[k]
                   + f_5 * fg_35[k]
                   + pa_y[k] * gg_65[k]
                   + f_2 * hg_s_110[k];

        t_111[k] = pa_z[k] * gg_51[k]
                   + f_2 * hg_s_111[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, gf_77, gf_78, gf_79, hg_s_112, hg_s_113, \
                         hg_s_114, hf_77, hf_78, hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_6 * gf_77[k]
                   + f_2 * hg_s_112[k]
                   + pb_x[k] * hf_77[k];

        t_113[k] = f_6 * gf_78[k]
                   + f_2 * hg_s_113[k]
                   + pb_x[k] * hf_78[k];

        t_114[k] = f_6 * gf_79[k]
                   + f_2 * hg_s_114[k]
                   + pb_x[k] * hf_79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_x, pa_z, pb_z, fg_s_117, fg_117, gf_36, \
                         gg_55, gg_117, hg_s_115, hg_s_116, hg_s_117, \
                         hf_76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * gg_55[k]
                   + f_2 * hg_s_115[k];

        t_116[k] = f_5 * gf_36[k]
                   + f_2 * hg_s_116[k]
                   + pb_z[k] * hf_76[k];

        t_117[k] = -f_10 * fg_s_117[k]
                   + f_5 * fg_117[k]
                   + pa_x[k] * gg_117[k]
                   + f_2 * hg_s_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_x, pa_y, pb_y, fg_s_119, fg_119, gf_49, \
                         gg_75, gg_119, hg_s_118, hg_s_119, hg_s_120, \
                         hf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_6 * gf_49[k]
                   + f_2 * hg_s_118[k]
                   + pb_y[k] * hf_79[k];

        t_119[k] = -f_10 * fg_s_119[k]
                   + f_5 * fg_119[k]
                   + pa_x[k] * gg_119[k]
                   + f_2 * hg_s_119[k];

        t_120[k] = pa_y[k] * gg_75[k]
                   + f_2 * hg_s_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pb_y, gf_50, gf_51, gg_77, gg_78, \
                         hg_s_121, hg_s_122, hg_s_123, hf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * gf_50[k]
                   + f_2 * hg_s_121[k]
                   + pb_y[k] * hf_80[k];

        t_122[k] = pa_y[k] * gg_77[k]
                   + f_2 * hg_s_122[k];

        t_123[k] = f_6 * gf_51[k]
                   + pa_y[k] * gg_78[k]
                   + f_2 * hg_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_x, pb_y, gf_52, gf_86, gg_80, hg_s_124, \
                         hg_s_125, hg_s_126, hf_82, hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * gf_52[k]
                   + f_2 * hg_s_124[k]
                   + pb_y[k] * hf_82[k];

        t_125[k] = pa_y[k] * gg_80[k]
                   + f_2 * hg_s_125[k];

        t_126[k] = f_6 * gf_86[k]
                   + f_2 * hg_s_126[k]
                   + pb_x[k] * hf_86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_y, pb_x, gf_87, gf_88, gg_84, hg_s_127, \
                         hg_s_128, hg_s_129, hf_87, hf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_6 * gf_87[k]
                   + f_2 * hg_s_127[k]
                   + pb_x[k] * hf_87[k];

        t_128[k] = f_6 * gf_88[k]
                   + f_2 * hg_s_128[k]
                   + pb_x[k] * hf_88[k];

        t_129[k] = pa_y[k] * gg_84[k]
                   + f_2 * hg_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_x, pb_z, fg_s_130, fg_s_132, fg_130, fg_132, \
                         gf_46, gg_130, gg_132, hg_s_130, hg_s_131, hg_s_132, \
                         hf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -f_10 * fg_s_130[k]
                   + f_5 * fg_130[k]
                   + pa_x[k] * gg_130[k]
                   + f_2 * hg_s_130[k];

        t_131[k] = f_6 * gf_46[k]
                   + f_2 * hg_s_131[k]
                   + pb_z[k] * hf_86[k];

        t_132[k] = -f_10 * fg_s_132[k]
                   + f_5 * fg_132[k]
                   + pa_x[k] * gg_132[k]
                   + f_2 * hg_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_y, pa_z, pb_y, fg_s_30, fg_30, gf_59, gg_75, \
                         gg_89, hg_s_133, hg_s_134, hg_s_135, hf_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * gf_59[k]
                   + f_2 * hg_s_133[k]
                   + pb_y[k] * hf_89[k];

        t_134[k] = pa_y[k] * gg_89[k]
                   + f_2 * hg_s_134[k];

        t_135[k] = -f_11 * fg_s_30[k]
                   + f_6 * fg_30[k]
                   + pa_z[k] * gg_75[k]
                   + f_2 * hg_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_y, pb_z, gf_50, hd_s_54, hg_s_136, \
                         hg_s_137, hg_s_138, hg_s_139, hd_54, hf_90, hf_91, \
                         hf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * hg_s_136[k]
                   + pb_y[k] * hf_90[k];

        t_137[k] = f_3 * gf_50[k]
                   + f_2 * hg_s_137[k]
                   + pb_z[k] * hf_90[k];

        t_138[k] = -f_4 * hd_s_54[k]
                   + f_2 * hg_s_138[k]
                   + f_5 * hd_54[k]
                   + pb_y[k] * hf_91[k];

        t_139[k] = f_2 * hg_s_139[k]
                   + pb_y[k] * hf_92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, gf_95, gf_96, gf_97, hd_s_59, hg_s_140, \
                         hg_s_141, hg_s_142, hd_59, hf_95, hf_96, \
                         hf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_6 * gf_95[k]
                   - f_4 * hd_s_59[k]
                   + f_2 * hg_s_140[k]
                   + f_5 * hd_59[k]
                   + pb_x[k] * hf_95[k];

        t_141[k] = f_6 * gf_96[k]
                   + f_2 * hg_s_141[k]
                   + pb_x[k] * hf_96[k];

        t_142[k] = f_6 * gf_97[k]
                   + f_2 * hg_s_142[k]
                   + pb_x[k] * hf_97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, pb_y, gf_99, hd_s_57, hg_s_143, hg_s_144, \
                         hg_s_145, hd_57, hf_95, hf_96, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_2 * hg_s_143[k]
                   + pb_y[k] * hf_95[k];

        t_144[k] = f_6 * gf_99[k]
                   + f_2 * hg_s_144[k]
                   + pb_x[k] * hf_99[k];

        t_145[k] = -f_1 * hd_s_57[k]
                   + f_2 * hg_s_145[k]
                   + f_3 * hd_57[k]
                   + pb_y[k] * hf_96[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pb_y, hd_s_58, hd_s_59, hg_s_146, hg_s_147, \
                         hg_s_148, hd_58, hd_59, hf_97, hf_98, hf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -f_9 * hd_s_58[k]
                   + f_2 * hg_s_146[k]
                   + f_6 * hd_58[k]
                   + pb_y[k] * hf_97[k];

        t_147[k] = -f_4 * hd_s_59[k]
                   + f_2 * hg_s_147[k]
                   + f_5 * hd_59[k]
                   + pb_y[k] * hf_98[k];

        t_148[k] = f_2 * hg_s_148[k]
                   + pb_y[k] * hf_99[k];
    }

#pragma omp simd aligned(t_149, t_150, t_151, pa_x, pb_y, fg_s_149, fg_149, gf_60, gf_100, \
                         gg_149, gg_150, hg_s_149, hg_s_150, hg_s_151, \
                         hf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -f_10 * fg_s_149[k]
                   + f_5 * fg_149[k]
                   + pa_x[k] * gg_149[k]
                   + f_2 * hg_s_149[k];

        t_150[k] = f_7 * gf_100[k]
                   + pa_x[k] * gg_150[k]
                   + f_2 * hg_s_150[k];

        t_151[k] = f_7 * gf_60[k]
                   + f_2 * hg_s_151[k]
                   + pb_y[k] * hf_100[k];
    }

#pragma omp simd aligned(t_152, t_153, t_154, t_155, pa_x, pb_z, gf_103, gf_105, gg_153, \
                         gg_155, hg_s_152, hg_s_153, hg_s_154, hg_s_155, hf_100, \
                         hf_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_152[k] = f_2 * hg_s_152[k]
                   + pb_z[k] * hf_100[k];

        t_153[k] = f_6 * gf_103[k]
                   + pa_x[k] * gg_153[k]
                   + f_2 * hg_s_153[k];

        t_154[k] = f_2 * hg_s_154[k]
                   + pb_z[k] * hf_101[k];

        t_155[k] = f_6 * gf_105[k]
                   + pa_x[k] * gg_155[k]
                   + f_2 * hg_s_155[k];
    }

#pragma omp simd aligned(t_156, t_157, t_158, pb_x, pb_z, gf_106, gf_108, hg_s_156, hg_s_157, \
                         hg_s_158, hf_103, hf_106, hf_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_156[k] = f_5 * gf_106[k]
                   + f_2 * hg_s_156[k]
                   + pb_x[k] * hf_106[k];

        t_157[k] = f_2 * hg_s_157[k]
                   + pb_z[k] * hf_103[k];

        t_158[k] = f_5 * gf_108[k]
                   + f_2 * hg_s_158[k]
                   + pb_x[k] * hf_108[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pb_x, pb_z, gf_109, gg_160, gg_162, \
                         hg_s_159, hg_s_160, hg_s_161, hg_s_162, hf_106, \
                         hf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = f_5 * gf_109[k]
                   + f_2 * hg_s_159[k]
                   + pb_x[k] * hf_109[k];

        t_160[k] = pa_x[k] * gg_160[k]
                   + f_2 * hg_s_160[k];

        t_161[k] = f_2 * hg_s_161[k]
                   + pb_z[k] * hf_106[k];

        t_162[k] = pa_x[k] * gg_162[k]
                   + f_2 * hg_s_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, t_166, pa_x, pa_z, gg_90, gg_91, gg_163, gg_164, \
                         hg_s_163, hg_s_164, hg_s_165, hg_s_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_x[k] * gg_163[k]
                   + f_2 * hg_s_163[k];

        t_164[k] = pa_x[k] * gg_164[k]
                   + f_2 * hg_s_164[k];

        t_165[k] = pa_z[k] * gg_90[k]
                   + f_2 * hg_s_165[k];

        t_166[k] = pa_z[k] * gg_91[k]
                   + f_2 * hg_s_166[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_z, pb_y, pb_z, gf_60, gf_72, gg_93, hg_s_167, \
                         hg_s_168, hg_s_169, hf_110, hf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = f_5 * gf_60[k]
                   + f_2 * hg_s_167[k]
                   + pb_z[k] * hf_110[k];

        t_168[k] = pa_z[k] * gg_93[k]
                   + f_2 * hg_s_168[k];

        t_169[k] = f_3 * gf_72[k]
                   + f_2 * hg_s_169[k]
                   + pb_y[k] * hf_112[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_x, pa_z, pb_x, gf_115, gf_117, gg_96, gg_170, \
                         hg_s_170, hg_s_171, hg_s_172, hf_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_6 * gf_115[k]
                   + pa_x[k] * gg_170[k]
                   + f_2 * hg_s_170[k];

        t_171[k] = pa_z[k] * gg_96[k]
                   + f_2 * hg_s_171[k];

        t_172[k] = f_5 * gf_117[k]
                   + f_2 * hg_s_172[k]
                   + pb_x[k] * hf_117[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, t_176, pa_x, pb_x, gf_118, gf_119, gg_175, \
                         gg_176, hg_s_173, hg_s_174, hg_s_175, hg_s_176, hf_118, \
                         hf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_5 * gf_118[k]
                   + f_2 * hg_s_173[k]
                   + pb_x[k] * hf_118[k];

        t_174[k] = f_5 * gf_119[k]
                   + f_2 * hg_s_174[k]
                   + pb_x[k] * hf_119[k];

        t_175[k] = pa_x[k] * gg_175[k]
                   + f_2 * hg_s_175[k];

        t_176[k] = pa_x[k] * gg_176[k]
                   + f_2 * hg_s_176[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, t_180, pa_x, gf_120, gg_177, gg_178, gg_179, \
                         gg_180, hg_s_177, hg_s_178, hg_s_179, \
                         hg_s_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_x[k] * gg_177[k]
                   + f_2 * hg_s_177[k];

        t_178[k] = pa_x[k] * gg_178[k]
                   + f_2 * hg_s_178[k];

        t_179[k] = pa_x[k] * gg_179[k]
                   + f_2 * hg_s_179[k];

        t_180[k] = f_7 * gf_120[k]
                   + pa_x[k] * gg_180[k]
                   + f_2 * hg_s_180[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pa_x, pb_y, pb_z, gf_70, gf_80, gf_123, gg_183, \
                         hg_s_181, hg_s_182, hg_s_183, hf_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_6 * gf_80[k]
                   + f_2 * hg_s_181[k]
                   + pb_y[k] * hf_120[k];

        t_182[k] = f_6 * gf_70[k]
                   + f_2 * hg_s_182[k]
                   + pb_z[k] * hf_120[k];

        t_183[k] = f_6 * gf_123[k]
                   + pa_x[k] * gg_183[k]
                   + f_2 * hg_s_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_x, pb_x, pb_y, gf_82, gf_125, gf_126, gg_185, \
                         hg_s_184, hg_s_185, hg_s_186, hf_122, hf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_6 * gf_82[k]
                   + f_2 * hg_s_184[k]
                   + pb_y[k] * hf_122[k];

        t_185[k] = f_6 * gf_125[k]
                   + pa_x[k] * gg_185[k]
                   + f_2 * hg_s_185[k];

        t_186[k] = f_5 * gf_126[k]
                   + f_2 * hg_s_186[k]
                   + pb_x[k] * hf_126[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pb_x, gf_127, gf_128, gf_129, hg_s_187, \
                         hg_s_188, hg_s_189, hf_127, hf_128, hf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_5 * gf_127[k]
                   + f_2 * hg_s_187[k]
                   + pb_x[k] * hf_127[k];

        t_188[k] = f_5 * gf_128[k]
                   + f_2 * hg_s_188[k]
                   + pb_x[k] * hf_128[k];

        t_189[k] = f_5 * gf_129[k]
                   + f_2 * hg_s_189[k]
                   + pb_x[k] * hf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, t_194, pa_x, gg_190, gg_191, gg_192, \
                         gg_193, gg_194, hg_s_190, hg_s_191, hg_s_192, hg_s_193, \
                         hg_s_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_x[k] * gg_190[k]
                   + f_2 * hg_s_190[k];

        t_191[k] = pa_x[k] * gg_191[k]
                   + f_2 * hg_s_191[k];

        t_192[k] = pa_x[k] * gg_192[k]
                   + f_2 * hg_s_192[k];

        t_193[k] = pa_x[k] * gg_193[k]
                   + f_2 * hg_s_193[k];

        t_194[k] = pa_x[k] * gg_194[k]
                   + f_2 * hg_s_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_y, pb_y, gf_90, gg_135, gg_137, hg_s_195, \
                         hg_s_196, hg_s_197, hf_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = pa_y[k] * gg_135[k]
                   + f_2 * hg_s_195[k];

        t_196[k] = f_5 * gf_90[k]
                   + f_2 * hg_s_196[k]
                   + pb_y[k] * hf_130[k];

        t_197[k] = pa_y[k] * gg_137[k]
                   + f_2 * hg_s_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_x, pa_y, pb_y, gf_92, gf_133, gg_140, gg_198, \
                         hg_s_198, hg_s_199, hg_s_200, hf_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_6 * gf_133[k]
                   + pa_x[k] * gg_198[k]
                   + f_2 * hg_s_198[k];

        t_199[k] = f_5 * gf_92[k]
                   + f_2 * hg_s_199[k]
                   + pb_y[k] * hf_132[k];

        t_200[k] = pa_y[k] * gg_140[k]
                   + f_2 * hg_s_200[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, gf_136, gf_137, gf_138, hg_s_201, \
                         hg_s_202, hg_s_203, hf_136, hf_137, hf_138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_5 * gf_136[k]
                   + f_2 * hg_s_201[k]
                   + pb_x[k] * hf_136[k];

        t_202[k] = f_5 * gf_137[k]
                   + f_2 * hg_s_202[k]
                   + pb_x[k] * hf_137[k];

        t_203[k] = f_5 * gf_138[k]
                   + f_2 * hg_s_203[k]
                   + pb_x[k] * hf_138[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, pa_y, gg_144, gg_205, gg_206, \
                         gg_207, hg_s_204, hg_s_205, hg_s_206, \
                         hg_s_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * gg_144[k]
                   + f_2 * hg_s_204[k];

        t_205[k] = pa_x[k] * gg_205[k]
                   + f_2 * hg_s_205[k];

        t_206[k] = pa_x[k] * gg_206[k]
                   + f_2 * hg_s_206[k];

        t_207[k] = pa_x[k] * gg_207[k]
                   + f_2 * hg_s_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, t_211, pa_x, pb_y, gf_140, gg_208, gg_209, \
                         gg_210, hg_s_208, hg_s_209, hg_s_210, hg_s_211, \
                         hf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = pa_x[k] * gg_208[k]
                   + f_2 * hg_s_208[k];

        t_209[k] = pa_x[k] * gg_209[k]
                   + f_2 * hg_s_209[k];

        t_210[k] = f_7 * gf_140[k]
                   + pa_x[k] * gg_210[k]
                   + f_2 * hg_s_210[k];

        t_211[k] = f_2 * hg_s_211[k]
                   + pb_y[k] * hf_140[k];
    }
}

static auto
compute_prim_hg_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fg_s,
                                        const size_t fg, const size_t gf, const size_t gg,
                                        const size_t hd_s, const size_t hg_s, const size_t hd,
                                        const size_t hf, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = 3.0 * beta / p;
    const auto f_9 = 2.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fg_s_100 = buffer.data(fg_s + 100);
    const auto *fg_s_115 = buffer.data(fg_s + 115);
    const auto *fg_s_119 = buffer.data(fg_s + 119);
    const auto *fg_s_134 = buffer.data(fg_s + 134);
    const auto *fg_s_149 = buffer.data(fg_s + 149);

    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *gf_90 = buffer.data(gf + 90);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_141 = buffer.data(gf + 141);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *hd_s_90 = buffer.data(hd_s + 90);
    const auto *hd_s_91 = buffer.data(hd_s + 91);
    const auto *hd_s_93 = buffer.data(hd_s + 93);
    const auto *hd_s_95 = buffer.data(hd_s + 95);
    const auto *hd_s_98 = buffer.data(hd_s + 98);
    const auto *hd_s_101 = buffer.data(hd_s + 101);
    const auto *hd_s_102 = buffer.data(hd_s + 102);
    const auto *hd_s_103 = buffer.data(hd_s + 103);
    const auto *hd_s_104 = buffer.data(hd_s + 104);
    const auto *hd_s_105 = buffer.data(hd_s + 105);
    const auto *hd_s_106 = buffer.data(hd_s + 106);
    const auto *hd_s_107 = buffer.data(hd_s + 107);
    const auto *hd_s_108 = buffer.data(hd_s + 108);
    const auto *hd_s_109 = buffer.data(hd_s + 109);
    const auto *hd_s_110 = buffer.data(hd_s + 110);
    const auto *hd_s_111 = buffer.data(hd_s + 111);
    const auto *hd_s_112 = buffer.data(hd_s + 112);
    const auto *hd_s_113 = buffer.data(hd_s + 113);
    const auto *hd_s_120 = buffer.data(hd_s + 120);
    const auto *hd_s_122 = buffer.data(hd_s + 122);
    const auto *hd_s_123 = buffer.data(hd_s + 123);
    const auto *hd_s_124 = buffer.data(hd_s + 124);
    const auto *hd_s_125 = buffer.data(hd_s + 125);

    const auto *hg_s_212 = buffer.data(hg_s + 212);
    const auto *hg_s_213 = buffer.data(hg_s + 213);
    const auto *hg_s_214 = buffer.data(hg_s + 214);
    const auto *hg_s_215 = buffer.data(hg_s + 215);
    const auto *hg_s_216 = buffer.data(hg_s + 216);
    const auto *hg_s_217 = buffer.data(hg_s + 217);
    const auto *hg_s_218 = buffer.data(hg_s + 218);
    const auto *hg_s_219 = buffer.data(hg_s + 219);
    const auto *hg_s_220 = buffer.data(hg_s + 220);
    const auto *hg_s_221 = buffer.data(hg_s + 221);
    const auto *hg_s_222 = buffer.data(hg_s + 222);
    const auto *hg_s_223 = buffer.data(hg_s + 223);
    const auto *hg_s_224 = buffer.data(hg_s + 224);
    const auto *hg_s_225 = buffer.data(hg_s + 225);
    const auto *hg_s_226 = buffer.data(hg_s + 226);
    const auto *hg_s_227 = buffer.data(hg_s + 227);
    const auto *hg_s_228 = buffer.data(hg_s + 228);
    const auto *hg_s_229 = buffer.data(hg_s + 229);
    const auto *hg_s_230 = buffer.data(hg_s + 230);
    const auto *hg_s_231 = buffer.data(hg_s + 231);
    const auto *hg_s_232 = buffer.data(hg_s + 232);
    const auto *hg_s_233 = buffer.data(hg_s + 233);
    const auto *hg_s_234 = buffer.data(hg_s + 234);
    const auto *hg_s_235 = buffer.data(hg_s + 235);
    const auto *hg_s_236 = buffer.data(hg_s + 236);
    const auto *hg_s_237 = buffer.data(hg_s + 237);
    const auto *hg_s_238 = buffer.data(hg_s + 238);
    const auto *hg_s_239 = buffer.data(hg_s + 239);
    const auto *hg_s_240 = buffer.data(hg_s + 240);
    const auto *hg_s_241 = buffer.data(hg_s + 241);
    const auto *hg_s_242 = buffer.data(hg_s + 242);
    const auto *hg_s_243 = buffer.data(hg_s + 243);
    const auto *hg_s_244 = buffer.data(hg_s + 244);
    const auto *hg_s_245 = buffer.data(hg_s + 245);
    const auto *hg_s_246 = buffer.data(hg_s + 246);
    const auto *hg_s_247 = buffer.data(hg_s + 247);
    const auto *hg_s_248 = buffer.data(hg_s + 248);
    const auto *hg_s_249 = buffer.data(hg_s + 249);
    const auto *hg_s_250 = buffer.data(hg_s + 250);
    const auto *hg_s_251 = buffer.data(hg_s + 251);
    const auto *hg_s_252 = buffer.data(hg_s + 252);
    const auto *hg_s_253 = buffer.data(hg_s + 253);
    const auto *hg_s_254 = buffer.data(hg_s + 254);
    const auto *hg_s_255 = buffer.data(hg_s + 255);
    const auto *hg_s_256 = buffer.data(hg_s + 256);
    const auto *hg_s_257 = buffer.data(hg_s + 257);
    const auto *hg_s_258 = buffer.data(hg_s + 258);
    const auto *hg_s_259 = buffer.data(hg_s + 259);
    const auto *hg_s_260 = buffer.data(hg_s + 260);
    const auto *hg_s_261 = buffer.data(hg_s + 261);
    const auto *hg_s_262 = buffer.data(hg_s + 262);
    const auto *hg_s_263 = buffer.data(hg_s + 263);
    const auto *hg_s_264 = buffer.data(hg_s + 264);
    const auto *hg_s_265 = buffer.data(hg_s + 265);
    const auto *hg_s_266 = buffer.data(hg_s + 266);
    const auto *hg_s_267 = buffer.data(hg_s + 267);
    const auto *hg_s_268 = buffer.data(hg_s + 268);
    const auto *hg_s_269 = buffer.data(hg_s + 269);
    const auto *hg_s_270 = buffer.data(hg_s + 270);
    const auto *hg_s_271 = buffer.data(hg_s + 271);
    const auto *hg_s_272 = buffer.data(hg_s + 272);
    const auto *hg_s_273 = buffer.data(hg_s + 273);
    const auto *hg_s_274 = buffer.data(hg_s + 274);
    const auto *hg_s_275 = buffer.data(hg_s + 275);
    const auto *hg_s_276 = buffer.data(hg_s + 276);
    const auto *hg_s_277 = buffer.data(hg_s + 277);
    const auto *hg_s_278 = buffer.data(hg_s + 278);
    const auto *hg_s_279 = buffer.data(hg_s + 279);
    const auto *hg_s_280 = buffer.data(hg_s + 280);
    const auto *hg_s_281 = buffer.data(hg_s + 281);
    const auto *hg_s_282 = buffer.data(hg_s + 282);
    const auto *hg_s_283 = buffer.data(hg_s + 283);
    const auto *hg_s_284 = buffer.data(hg_s + 284);
    const auto *hg_s_285 = buffer.data(hg_s + 285);
    const auto *hg_s_286 = buffer.data(hg_s + 286);
    const auto *hg_s_287 = buffer.data(hg_s + 287);
    const auto *hg_s_288 = buffer.data(hg_s + 288);
    const auto *hg_s_289 = buffer.data(hg_s + 289);
    const auto *hg_s_290 = buffer.data(hg_s + 290);
    const auto *hg_s_291 = buffer.data(hg_s + 291);
    const auto *hg_s_292 = buffer.data(hg_s + 292);
    const auto *hg_s_293 = buffer.data(hg_s + 293);
    const auto *hg_s_294 = buffer.data(hg_s + 294);
    const auto *hg_s_295 = buffer.data(hg_s + 295);
    const auto *hg_s_296 = buffer.data(hg_s + 296);
    const auto *hg_s_297 = buffer.data(hg_s + 297);
    const auto *hg_s_298 = buffer.data(hg_s + 298);
    const auto *hg_s_299 = buffer.data(hg_s + 299);
    const auto *hg_s_300 = buffer.data(hg_s + 300);
    const auto *hg_s_301 = buffer.data(hg_s + 301);
    const auto *hg_s_302 = buffer.data(hg_s + 302);
    const auto *hg_s_303 = buffer.data(hg_s + 303);
    const auto *hg_s_304 = buffer.data(hg_s + 304);
    const auto *hg_s_305 = buffer.data(hg_s + 305);
    const auto *hg_s_306 = buffer.data(hg_s + 306);
    const auto *hg_s_307 = buffer.data(hg_s + 307);
    const auto *hg_s_308 = buffer.data(hg_s + 308);
    const auto *hg_s_309 = buffer.data(hg_s + 309);
    const auto *hg_s_310 = buffer.data(hg_s + 310);
    const auto *hg_s_311 = buffer.data(hg_s + 311);
    const auto *hg_s_312 = buffer.data(hg_s + 312);
    const auto *hg_s_313 = buffer.data(hg_s + 313);
    const auto *hg_s_314 = buffer.data(hg_s + 314);

    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_91 = buffer.data(hd + 91);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_98 = buffer.data(hd + 98);
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
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_122 = buffer.data(hd + 122);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_124 = buffer.data(hd + 124);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_171 = buffer.data(hf + 171);
    const auto *hf_172 = buffer.data(hf + 172);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_174 = buffer.data(hf + 174);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_181 = buffer.data(hf + 181);
    const auto *hf_182 = buffer.data(hf + 182);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_184 = buffer.data(hf + 184);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

#pragma omp simd aligned(t_212, t_213, t_214, pa_x, pb_y, pb_z, gf_90, gf_143, gg_213, \
                         hg_s_212, hg_s_213, hg_s_214, hf_140, hf_142 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = f_7 * gf_90[k]
                   + f_2 * hg_s_212[k]
                   + pb_z[k] * hf_140[k];

        t_213[k] = f_6 * gf_143[k]
                   + pa_x[k] * gg_213[k]
                   + f_2 * hg_s_213[k];

        t_214[k] = f_2 * hg_s_214[k]
                   + pb_y[k] * hf_142[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pa_x, pb_x, gf_145, gf_146, gf_147, gg_215, \
                         hg_s_215, hg_s_216, hg_s_217, hf_146, hf_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_6 * gf_145[k]
                   + pa_x[k] * gg_215[k]
                   + f_2 * hg_s_215[k];

        t_216[k] = f_5 * gf_146[k]
                   + f_2 * hg_s_216[k]
                   + pb_x[k] * hf_146[k];

        t_217[k] = f_5 * gf_147[k]
                   + f_2 * hg_s_217[k]
                   + pb_x[k] * hf_147[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, t_221, pa_x, pb_x, pb_y, gf_149, gg_220, gg_221, \
                         hg_s_218, hg_s_219, hg_s_220, hg_s_221, hf_145, \
                         hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_2 * hg_s_218[k]
                   + pb_y[k] * hf_145[k];

        t_219[k] = f_5 * gf_149[k]
                   + f_2 * hg_s_219[k]
                   + pb_x[k] * hf_149[k];

        t_220[k] = pa_x[k] * gg_220[k]
                   + f_2 * hg_s_220[k];

        t_221[k] = pa_x[k] * gg_221[k]
                   + f_2 * hg_s_221[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_x, pb_y, gg_222, gg_224, hg_s_222, hg_s_223, \
                         hg_s_224, hf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = pa_x[k] * gg_222[k]
                   + f_2 * hg_s_222[k];

        t_223[k] = f_2 * hg_s_223[k]
                   + pb_y[k] * hf_149[k];

        t_224[k] = pa_x[k] * gg_224[k]
                   + f_2 * hg_s_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pb_x, pb_z, hd_s_90, hd_s_91, hg_s_225, \
                         hg_s_226, hg_s_227, hd_90, hd_91, hf_150, \
                         hf_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = -f_1 * hd_s_90[k]
                   + f_2 * hg_s_225[k]
                   + f_3 * hd_90[k]
                   + pb_x[k] * hf_150[k];

        t_226[k] = -f_9 * hd_s_91[k]
                   + f_2 * hg_s_226[k]
                   + f_6 * hd_91[k]
                   + pb_x[k] * hf_151[k];

        t_227[k] = f_2 * hg_s_227[k]
                   + pb_z[k] * hf_150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pb_x, pb_z, hd_s_93, hd_s_95, hg_s_228, \
                         hg_s_229, hg_s_230, hd_93, hd_95, hf_151, hf_153, \
                         hf_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = -f_4 * hd_s_93[k]
                   + f_2 * hg_s_228[k]
                   + f_5 * hd_93[k]
                   + pb_x[k] * hf_153[k];

        t_229[k] = f_2 * hg_s_229[k]
                   + pb_z[k] * hf_151[k];

        t_230[k] = -f_4 * hd_s_95[k]
                   + f_2 * hg_s_230[k]
                   + f_5 * hd_95[k]
                   + pb_x[k] * hf_155[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pb_x, hg_s_231, hg_s_232, hg_s_233, \
                         hg_s_234, hf_156, hf_157, hf_158, hf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_2 * hg_s_231[k]
                   + pb_x[k] * hf_156[k];

        t_232[k] = f_2 * hg_s_232[k]
                   + pb_x[k] * hf_157[k];

        t_233[k] = f_2 * hg_s_233[k]
                   + pb_x[k] * hf_158[k];

        t_234[k] = f_2 * hg_s_234[k]
                   + pb_x[k] * hf_159[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pb_y, pb_z, gf_106, hd_s_93, hg_s_235, hg_s_236, \
                         hg_s_237, hd_93, hf_156, hf_157 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_0 * gf_106[k]
                   - f_1 * hd_s_93[k]
                   + f_2 * hg_s_235[k]
                   + f_3 * hd_93[k]
                   + pb_y[k] * hf_156[k];

        t_236[k] = f_2 * hg_s_236[k]
                   + pb_z[k] * hf_156[k];

        t_237[k] = -f_4 * hd_s_93[k]
                   + f_2 * hg_s_237[k]
                   + f_5 * hd_93[k]
                   + pb_z[k] * hf_157[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_z, pb_y, pb_z, gf_109, gg_150, hd_s_95, \
                         hg_s_238, hg_s_239, hg_s_240, hd_95, hf_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_0 * gf_109[k]
                   + f_2 * hg_s_238[k]
                   + pb_y[k] * hf_159[k];

        t_239[k] = -f_1 * hd_s_95[k]
                   + f_2 * hg_s_239[k]
                   + f_3 * hd_95[k]
                   + pb_z[k] * hf_159[k];

        t_240[k] = pa_z[k] * gg_150[k]
                   + f_2 * hg_s_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, pa_z, pb_x, gg_151, gg_153, hd_s_98, hg_s_241, \
                         hg_s_242, hg_s_243, hd_98, hf_162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = pa_z[k] * gg_151[k]
                   + f_2 * hg_s_241[k];

        t_242[k] = -f_9 * hd_s_98[k]
                   + f_2 * hg_s_242[k]
                   + f_6 * hd_98[k]
                   + pb_x[k] * hf_162[k];

        t_243[k] = pa_z[k] * gg_153[k]
                   + f_2 * hg_s_243[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, pa_z, pb_x, gf_101, gg_154, hd_s_101, hg_s_244, \
                         hg_s_245, hg_s_246, hd_101, hf_165, hf_166 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_5 * gf_101[k]
                   + pa_z[k] * gg_154[k]
                   + f_2 * hg_s_244[k];

        t_245[k] = -f_4 * hd_s_101[k]
                   + f_2 * hg_s_245[k]
                   + f_5 * hd_101[k]
                   + pb_x[k] * hf_165[k];

        t_246[k] = f_2 * hg_s_246[k]
                   + pb_x[k] * hf_166[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, t_250, pa_z, pb_x, gg_160, hg_s_247, hg_s_248, \
                         hg_s_249, hg_s_250, hf_167, hf_168, hf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_2 * hg_s_247[k]
                   + pb_x[k] * hf_167[k];

        t_248[k] = f_2 * hg_s_248[k]
                   + pb_x[k] * hf_168[k];

        t_249[k] = f_2 * hg_s_249[k]
                   + pb_x[k] * hf_169[k];

        t_250[k] = pa_z[k] * gg_160[k]
                   + f_2 * hg_s_250[k];
    }

#pragma omp simd aligned(t_251, t_252, t_253, pa_z, pb_y, pb_z, gf_106, gf_107, gf_119, \
                         gg_162, hg_s_251, hg_s_252, hg_s_253, hf_166, \
                         hf_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_251[k] = f_5 * gf_106[k]
                   + f_2 * hg_s_251[k]
                   + pb_z[k] * hf_166[k];

        t_252[k] = f_6 * gf_107[k]
                   + pa_z[k] * gg_162[k]
                   + f_2 * hg_s_252[k];

        t_253[k] = f_7 * gf_119[k]
                   + f_2 * hg_s_253[k]
                   + pb_y[k] * hf_169[k];
    }

#pragma omp simd aligned(t_254, t_255, pa_y, pb_x, fg_s_119, fg_119, gg_179, hd_s_102, \
                         hg_s_254, hg_s_255, hd_102, hf_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_254[k] = -f_8 * fg_s_119[k]
                   + f_3 * fg_119[k]
                   + pa_y[k] * gg_179[k]
                   + f_2 * hg_s_254[k];

        t_255[k] = -f_1 * hd_s_102[k]
                   + f_2 * hg_s_255[k]
                   + f_3 * hd_102[k]
                   + pb_x[k] * hf_170[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pb_x, hd_s_103, hd_s_104, hd_s_105, hg_s_256, \
                         hg_s_257, hg_s_258, hd_103, hd_104, hd_105, hf_171, hf_172, \
                         hf_173 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = -f_9 * hd_s_103[k]
                   + f_2 * hg_s_256[k]
                   + f_6 * hd_103[k]
                   + pb_x[k] * hf_171[k];

        t_257[k] = -f_9 * hd_s_104[k]
                   + f_2 * hg_s_257[k]
                   + f_6 * hd_104[k]
                   + pb_x[k] * hf_172[k];

        t_258[k] = -f_4 * hd_s_105[k]
                   + f_2 * hg_s_258[k]
                   + f_5 * hd_105[k]
                   + pb_x[k] * hf_173[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pb_x, hd_s_106, hd_s_107, hg_s_259, hg_s_260, \
                         hg_s_261, hd_106, hd_107, hf_174, hf_175, \
                         hf_176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = -f_4 * hd_s_106[k]
                   + f_2 * hg_s_259[k]
                   + f_5 * hd_106[k]
                   + pb_x[k] * hf_174[k];

        t_260[k] = -f_4 * hd_s_107[k]
                   + f_2 * hg_s_260[k]
                   + f_5 * hd_107[k]
                   + pb_x[k] * hf_175[k];

        t_261[k] = f_2 * hg_s_261[k]
                   + pb_x[k] * hf_176[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, t_265, pa_z, pb_x, fg_s_100, fg_100, gg_175, \
                         hg_s_262, hg_s_263, hg_s_264, hg_s_265, hf_177, hf_178, \
                         hf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_2 * hg_s_262[k]
                   + pb_x[k] * hf_177[k];

        t_263[k] = f_2 * hg_s_263[k]
                   + pb_x[k] * hf_178[k];

        t_264[k] = f_2 * hg_s_264[k]
                   + pb_x[k] * hf_179[k];

        t_265[k] = -f_10 * fg_s_100[k]
                   + f_5 * fg_100[k]
                   + pa_z[k] * gg_175[k]
                   + f_2 * hg_s_265[k];
    }

#pragma omp simd aligned(t_266, t_267, t_268, pb_y, pb_z, gf_116, gf_128, gf_129, hd_s_107, \
                         hg_s_266, hg_s_267, hg_s_268, hd_107, hf_176, hf_178, \
                         hf_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_266[k] = f_6 * gf_116[k]
                   + f_2 * hg_s_266[k]
                   + pb_z[k] * hf_176[k];

        t_267[k] = f_3 * gf_128[k]
                   - f_4 * hd_s_107[k]
                   + f_2 * hg_s_267[k]
                   + f_5 * hd_107[k]
                   + pb_y[k] * hf_178[k];

        t_268[k] = f_3 * gf_129[k]
                   + f_2 * hg_s_268[k]
                   + pb_y[k] * hf_179[k];
    }

#pragma omp simd aligned(t_269, t_270, pa_y, pb_x, fg_s_134, fg_134, gg_194, hd_s_108, \
                         hg_s_269, hg_s_270, hd_108, hf_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_269[k] = -f_11 * fg_s_134[k]
                   + f_6 * fg_134[k]
                   + pa_y[k] * gg_194[k]
                   + f_2 * hg_s_269[k];

        t_270[k] = -f_1 * hd_s_108[k]
                   + f_2 * hg_s_270[k]
                   + f_3 * hd_108[k]
                   + pb_x[k] * hf_180[k];
    }

#pragma omp simd aligned(t_271, t_272, t_273, pb_x, hd_s_109, hd_s_110, hd_s_111, hg_s_271, \
                         hg_s_272, hg_s_273, hd_109, hd_110, hd_111, hf_181, hf_182, \
                         hf_183 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_271[k] = -f_9 * hd_s_109[k]
                   + f_2 * hg_s_271[k]
                   + f_6 * hd_109[k]
                   + pb_x[k] * hf_181[k];

        t_272[k] = -f_9 * hd_s_110[k]
                   + f_2 * hg_s_272[k]
                   + f_6 * hd_110[k]
                   + pb_x[k] * hf_182[k];

        t_273[k] = -f_4 * hd_s_111[k]
                   + f_2 * hg_s_273[k]
                   + f_5 * hd_111[k]
                   + pb_x[k] * hf_183[k];
    }

#pragma omp simd aligned(t_274, t_275, t_276, pb_x, hd_s_112, hd_s_113, hg_s_274, hg_s_275, \
                         hg_s_276, hd_112, hd_113, hf_184, hf_185, \
                         hf_186 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_274[k] = -f_4 * hd_s_112[k]
                   + f_2 * hg_s_274[k]
                   + f_5 * hd_112[k]
                   + pb_x[k] * hf_184[k];

        t_275[k] = -f_4 * hd_s_113[k]
                   + f_2 * hg_s_275[k]
                   + f_5 * hd_113[k]
                   + pb_x[k] * hf_185[k];

        t_276[k] = f_2 * hg_s_276[k]
                   + pb_x[k] * hf_186[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, t_280, pa_z, pb_x, fg_s_115, fg_115, gg_190, \
                         hg_s_277, hg_s_278, hg_s_279, hg_s_280, hf_187, hf_188, \
                         hf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = f_2 * hg_s_277[k]
                   + pb_x[k] * hf_187[k];

        t_278[k] = f_2 * hg_s_278[k]
                   + pb_x[k] * hf_188[k];

        t_279[k] = f_2 * hg_s_279[k]
                   + pb_x[k] * hf_189[k];

        t_280[k] = -f_11 * fg_s_115[k]
                   + f_6 * fg_115[k]
                   + pa_z[k] * gg_190[k]
                   + f_2 * hg_s_280[k];
    }

#pragma omp simd aligned(t_281, t_282, t_283, pb_y, pb_z, gf_126, gf_138, gf_139, hd_s_113, \
                         hg_s_281, hg_s_282, hg_s_283, hd_113, hf_186, hf_188, \
                         hf_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_281[k] = f_3 * gf_126[k]
                   + f_2 * hg_s_281[k]
                   + pb_z[k] * hf_186[k];

        t_282[k] = f_6 * gf_138[k]
                   - f_4 * hd_s_113[k]
                   + f_2 * hg_s_282[k]
                   + f_5 * hd_113[k]
                   + pb_y[k] * hf_188[k];

        t_283[k] = f_6 * gf_139[k]
                   + f_2 * hg_s_283[k]
                   + pb_y[k] * hf_189[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pa_y, fg_s_149, fg_149, gf_140, gg_209, \
                         gg_210, gg_211, gg_212, hg_s_284, hg_s_285, hg_s_286, \
                         hg_s_287 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = -f_10 * fg_s_149[k]
                   + f_5 * fg_149[k]
                   + pa_y[k] * gg_209[k]
                   + f_2 * hg_s_284[k];

        t_285[k] = pa_y[k] * gg_210[k]
                   + f_2 * hg_s_285[k];

        t_286[k] = f_5 * gf_140[k]
                   + pa_y[k] * gg_211[k]
                   + f_2 * hg_s_286[k];

        t_287[k] = pa_y[k] * gg_212[k]
                   + f_2 * hg_s_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, pa_y, pb_x, gf_141, gf_142, gg_213, \
                         gg_214, gg_215, hg_s_288, hg_s_289, hg_s_290, hg_s_291, \
                         hf_196 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_6 * gf_141[k]
                   + pa_y[k] * gg_213[k]
                   + f_2 * hg_s_288[k];

        t_289[k] = f_5 * gf_142[k]
                   + pa_y[k] * gg_214[k]
                   + f_2 * hg_s_289[k];

        t_290[k] = pa_y[k] * gg_215[k]
                   + f_2 * hg_s_290[k];

        t_291[k] = f_2 * hg_s_291[k]
                   + pb_x[k] * hf_196[k];
    }

#pragma omp simd aligned(t_292, t_293, t_294, t_295, pa_y, pb_x, gf_146, gg_220, hg_s_292, \
                         hg_s_293, hg_s_294, hg_s_295, hf_197, hf_198, \
                         hf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_292[k] = f_2 * hg_s_292[k]
                   + pb_x[k] * hf_197[k];

        t_293[k] = f_2 * hg_s_293[k]
                   + pb_x[k] * hf_198[k];

        t_294[k] = f_2 * hg_s_294[k]
                   + pb_x[k] * hf_199[k];

        t_295[k] = f_7 * gf_146[k]
                   + pa_y[k] * gg_220[k]
                   + f_2 * hg_s_295[k];
    }

#pragma omp simd aligned(t_296, t_297, t_298, pa_y, pb_y, pb_z, gf_136, gf_148, gf_149, \
                         gg_222, hg_s_296, hg_s_297, hg_s_298, hf_196, \
                         hf_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_296[k] = f_7 * gf_136[k]
                   + f_2 * hg_s_296[k]
                   + pb_z[k] * hf_196[k];

        t_297[k] = f_6 * gf_148[k]
                   + pa_y[k] * gg_222[k]
                   + f_2 * hg_s_297[k];

        t_298[k] = f_5 * gf_149[k]
                   + f_2 * hg_s_298[k]
                   + pb_y[k] * hf_199[k];
    }

#pragma omp simd aligned(t_299, t_300, t_301, pa_y, pb_x, pb_y, gg_224, hd_s_120, hg_s_299, \
                         hg_s_300, hg_s_301, hd_120, hf_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_299[k] = pa_y[k] * gg_224[k]
                   + f_2 * hg_s_299[k];

        t_300[k] = -f_1 * hd_s_120[k]
                   + f_2 * hg_s_300[k]
                   + f_3 * hd_120[k]
                   + pb_x[k] * hf_200[k];

        t_301[k] = f_2 * hg_s_301[k]
                   + pb_y[k] * hf_200[k];
    }

#pragma omp simd aligned(t_302, t_303, t_304, pb_x, pb_y, hd_s_122, hd_s_123, hg_s_302, \
                         hg_s_303, hg_s_304, hd_122, hd_123, hf_202, \
                         hf_203 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_302[k] = -f_9 * hd_s_122[k]
                   + f_2 * hg_s_302[k]
                   + f_6 * hd_122[k]
                   + pb_x[k] * hf_202[k];

        t_303[k] = -f_4 * hd_s_123[k]
                   + f_2 * hg_s_303[k]
                   + f_5 * hd_123[k]
                   + pb_x[k] * hf_203[k];

        t_304[k] = f_2 * hg_s_304[k]
                   + pb_y[k] * hf_202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, t_308, pb_x, hd_s_125, hg_s_305, hg_s_306, \
                         hg_s_307, hg_s_308, hd_125, hf_205, hf_206, hf_207, \
                         hf_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = -f_4 * hd_s_125[k]
                   + f_2 * hg_s_305[k]
                   + f_5 * hd_125[k]
                   + pb_x[k] * hf_205[k];

        t_306[k] = f_2 * hg_s_306[k]
                   + pb_x[k] * hf_206[k];

        t_307[k] = f_2 * hg_s_307[k]
                   + pb_x[k] * hf_207[k];

        t_308[k] = f_2 * hg_s_308[k]
                   + pb_x[k] * hf_208[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, pb_x, pb_y, hd_s_123, hd_s_124, hg_s_309, \
                         hg_s_310, hg_s_311, hd_123, hd_124, hf_206, hf_207, \
                         hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = f_2 * hg_s_309[k]
                   + pb_x[k] * hf_209[k];

        t_310[k] = -f_1 * hd_s_123[k]
                   + f_2 * hg_s_310[k]
                   + f_3 * hd_123[k]
                   + pb_y[k] * hf_206[k];

        t_311[k] = -f_9 * hd_s_124[k]
                   + f_2 * hg_s_311[k]
                   + f_6 * hd_124[k]
                   + pb_y[k] * hf_207[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pb_y, pb_z, gf_149, hd_s_125, hg_s_312, \
                         hg_s_313, hg_s_314, hd_125, hf_208, hf_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = -f_4 * hd_s_125[k]
                   + f_2 * hg_s_312[k]
                   + f_5 * hd_125[k]
                   + pb_y[k] * hf_208[k];

        t_313[k] = f_2 * hg_s_313[k]
                   + pb_y[k] * hf_209[k];

        t_314[k] = f_0 * gf_149[k]
                   - f_1 * hd_s_125[k]
                   + f_2 * hg_s_314[k]
                   + f_3 * hd_125[k]
                   + pb_z[k] * hf_209[k];
    }
}

auto
compute_prim_hg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fg_s, const size_t fg,
                                 const size_t gf, const size_t gg, const size_t hd_s,
                                 const size_t hg_s, const size_t hd, const size_t hf,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_hg_kinetic_energy_0_piece0(buffer, target, pa, pb, fg_s, fg, gf, gg, hd_s,
                                            hg_s, hd, hf, ncols, alpha, beta, p);

    compute_prim_hg_kinetic_energy_0_piece1(buffer, target, pa, pb, fg_s, fg, gf, gg, hd_s,
                                            hg_s, hd, hf, ncols, alpha, beta, p);

    compute_prim_hg_kinetic_energy_0_piece2(buffer, target, pa, pb, fg_s, fg, gf, gg, hd_s,
                                            hg_s, hd, hf, ncols, alpha, beta, p);
}

}  // namespace simdkin
