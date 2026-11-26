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


#include "SimdKineticEnergyVrrRecGG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_gg_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t dg_s,
                                        const size_t dg, const size_t ff, const size_t fg,
                                        const size_t gd_s, const size_t gg_s, const size_t gd,
                                        const size_t gf, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 2.0 * alpha / p;
    const auto f_9 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_s_0 = buffer.data(dg_s + 0);
    const auto *dg_s_25 = buffer.data(dg_s + 25);
    const auto *dg_s_44 = buffer.data(dg_s + 44);
    const auto *dg_s_55 = buffer.data(dg_s + 55);
    const auto *dg_s_72 = buffer.data(dg_s + 72);
    const auto *dg_s_89 = buffer.data(dg_s + 89);

    const auto *dg_0 = buffer.data(dg + 0);
    const auto *dg_25 = buffer.data(dg + 25);
    const auto *dg_44 = buffer.data(dg + 44);
    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_72 = buffer.data(dg + 72);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *ff_0 = buffer.data(ff + 0);
    const auto *ff_1 = buffer.data(ff + 1);
    const auto *ff_2 = buffer.data(ff + 2);
    const auto *ff_6 = buffer.data(ff + 6);
    const auto *ff_9 = buffer.data(ff + 9);
    const auto *ff_10 = buffer.data(ff + 10);
    const auto *ff_16 = buffer.data(ff + 16);
    const auto *ff_18 = buffer.data(ff + 18);
    const auto *ff_19 = buffer.data(ff + 19);
    const auto *ff_20 = buffer.data(ff + 20);
    const auto *ff_22 = buffer.data(ff + 22);
    const auto *ff_27 = buffer.data(ff + 27);
    const auto *ff_29 = buffer.data(ff + 29);
    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_33 = buffer.data(ff + 33);
    const auto *ff_36 = buffer.data(ff + 36);
    const auto *ff_38 = buffer.data(ff + 38);
    const auto *ff_39 = buffer.data(ff + 39);
    const auto *ff_47 = buffer.data(ff + 47);
    const auto *ff_48 = buffer.data(ff + 48);
    const auto *ff_55 = buffer.data(ff + 55);
    const auto *ff_56 = buffer.data(ff + 56);
    const auto *ff_57 = buffer.data(ff + 57);
    const auto *ff_59 = buffer.data(ff + 59);
    const auto *ff_60 = buffer.data(ff + 60);
    const auto *ff_63 = buffer.data(ff + 63);
    const auto *ff_65 = buffer.data(ff + 65);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_68 = buffer.data(ff + 68);
    const auto *ff_69 = buffer.data(ff + 69);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_6 = buffer.data(fg + 6);
    const auto *fg_9 = buffer.data(fg + 9);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_16 = buffer.data(fg + 16);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_21 = buffer.data(fg + 21);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_39 = buffer.data(fg + 39);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_46 = buffer.data(fg + 46);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);

    const auto *gd_s_0 = buffer.data(gd_s + 0);
    const auto *gd_s_3 = buffer.data(gd_s + 3);
    const auto *gd_s_5 = buffer.data(gd_s + 5);
    const auto *gd_s_9 = buffer.data(gd_s + 9);
    const auto *gd_s_16 = buffer.data(gd_s + 16);
    const auto *gd_s_17 = buffer.data(gd_s + 17);
    const auto *gd_s_18 = buffer.data(gd_s + 18);
    const auto *gd_s_21 = buffer.data(gd_s + 21);
    const auto *gd_s_23 = buffer.data(gd_s + 23);
    const auto *gd_s_30 = buffer.data(gd_s + 30);
    const auto *gd_s_33 = buffer.data(gd_s + 33);
    const auto *gd_s_34 = buffer.data(gd_s + 34);
    const auto *gd_s_35 = buffer.data(gd_s + 35);

    const auto *gg_s_0 = buffer.data(gg_s + 0);
    const auto *gg_s_1 = buffer.data(gg_s + 1);
    const auto *gg_s_2 = buffer.data(gg_s + 2);
    const auto *gg_s_3 = buffer.data(gg_s + 3);
    const auto *gg_s_4 = buffer.data(gg_s + 4);
    const auto *gg_s_5 = buffer.data(gg_s + 5);
    const auto *gg_s_6 = buffer.data(gg_s + 6);
    const auto *gg_s_7 = buffer.data(gg_s + 7);
    const auto *gg_s_8 = buffer.data(gg_s + 8);
    const auto *gg_s_9 = buffer.data(gg_s + 9);
    const auto *gg_s_10 = buffer.data(gg_s + 10);
    const auto *gg_s_11 = buffer.data(gg_s + 11);
    const auto *gg_s_12 = buffer.data(gg_s + 12);
    const auto *gg_s_13 = buffer.data(gg_s + 13);
    const auto *gg_s_14 = buffer.data(gg_s + 14);
    const auto *gg_s_15 = buffer.data(gg_s + 15);
    const auto *gg_s_16 = buffer.data(gg_s + 16);
    const auto *gg_s_17 = buffer.data(gg_s + 17);
    const auto *gg_s_18 = buffer.data(gg_s + 18);
    const auto *gg_s_19 = buffer.data(gg_s + 19);
    const auto *gg_s_20 = buffer.data(gg_s + 20);
    const auto *gg_s_21 = buffer.data(gg_s + 21);
    const auto *gg_s_22 = buffer.data(gg_s + 22);
    const auto *gg_s_23 = buffer.data(gg_s + 23);
    const auto *gg_s_24 = buffer.data(gg_s + 24);
    const auto *gg_s_25 = buffer.data(gg_s + 25);
    const auto *gg_s_26 = buffer.data(gg_s + 26);
    const auto *gg_s_27 = buffer.data(gg_s + 27);
    const auto *gg_s_28 = buffer.data(gg_s + 28);
    const auto *gg_s_29 = buffer.data(gg_s + 29);
    const auto *gg_s_30 = buffer.data(gg_s + 30);
    const auto *gg_s_31 = buffer.data(gg_s + 31);
    const auto *gg_s_32 = buffer.data(gg_s + 32);
    const auto *gg_s_33 = buffer.data(gg_s + 33);
    const auto *gg_s_34 = buffer.data(gg_s + 34);
    const auto *gg_s_35 = buffer.data(gg_s + 35);
    const auto *gg_s_36 = buffer.data(gg_s + 36);
    const auto *gg_s_37 = buffer.data(gg_s + 37);
    const auto *gg_s_38 = buffer.data(gg_s + 38);
    const auto *gg_s_39 = buffer.data(gg_s + 39);
    const auto *gg_s_40 = buffer.data(gg_s + 40);
    const auto *gg_s_41 = buffer.data(gg_s + 41);
    const auto *gg_s_42 = buffer.data(gg_s + 42);
    const auto *gg_s_43 = buffer.data(gg_s + 43);
    const auto *gg_s_44 = buffer.data(gg_s + 44);
    const auto *gg_s_45 = buffer.data(gg_s + 45);
    const auto *gg_s_46 = buffer.data(gg_s + 46);
    const auto *gg_s_47 = buffer.data(gg_s + 47);
    const auto *gg_s_48 = buffer.data(gg_s + 48);
    const auto *gg_s_49 = buffer.data(gg_s + 49);
    const auto *gg_s_50 = buffer.data(gg_s + 50);
    const auto *gg_s_51 = buffer.data(gg_s + 51);
    const auto *gg_s_52 = buffer.data(gg_s + 52);
    const auto *gg_s_53 = buffer.data(gg_s + 53);
    const auto *gg_s_54 = buffer.data(gg_s + 54);
    const auto *gg_s_55 = buffer.data(gg_s + 55);
    const auto *gg_s_56 = buffer.data(gg_s + 56);
    const auto *gg_s_57 = buffer.data(gg_s + 57);
    const auto *gg_s_58 = buffer.data(gg_s + 58);
    const auto *gg_s_59 = buffer.data(gg_s + 59);
    const auto *gg_s_60 = buffer.data(gg_s + 60);
    const auto *gg_s_61 = buffer.data(gg_s + 61);
    const auto *gg_s_62 = buffer.data(gg_s + 62);
    const auto *gg_s_63 = buffer.data(gg_s + 63);
    const auto *gg_s_64 = buffer.data(gg_s + 64);
    const auto *gg_s_65 = buffer.data(gg_s + 65);
    const auto *gg_s_66 = buffer.data(gg_s + 66);
    const auto *gg_s_67 = buffer.data(gg_s + 67);
    const auto *gg_s_68 = buffer.data(gg_s + 68);
    const auto *gg_s_69 = buffer.data(gg_s + 69);
    const auto *gg_s_70 = buffer.data(gg_s + 70);
    const auto *gg_s_71 = buffer.data(gg_s + 71);
    const auto *gg_s_72 = buffer.data(gg_s + 72);
    const auto *gg_s_73 = buffer.data(gg_s + 73);
    const auto *gg_s_74 = buffer.data(gg_s + 74);
    const auto *gg_s_75 = buffer.data(gg_s + 75);
    const auto *gg_s_76 = buffer.data(gg_s + 76);
    const auto *gg_s_77 = buffer.data(gg_s + 77);
    const auto *gg_s_78 = buffer.data(gg_s + 78);
    const auto *gg_s_79 = buffer.data(gg_s + 79);
    const auto *gg_s_80 = buffer.data(gg_s + 80);
    const auto *gg_s_81 = buffer.data(gg_s + 81);
    const auto *gg_s_82 = buffer.data(gg_s + 82);
    const auto *gg_s_83 = buffer.data(gg_s + 83);
    const auto *gg_s_84 = buffer.data(gg_s + 84);
    const auto *gg_s_85 = buffer.data(gg_s + 85);
    const auto *gg_s_86 = buffer.data(gg_s + 86);
    const auto *gg_s_87 = buffer.data(gg_s + 87);
    const auto *gg_s_88 = buffer.data(gg_s + 88);
    const auto *gg_s_89 = buffer.data(gg_s + 89);
    const auto *gg_s_90 = buffer.data(gg_s + 90);
    const auto *gg_s_91 = buffer.data(gg_s + 91);
    const auto *gg_s_92 = buffer.data(gg_s + 92);
    const auto *gg_s_93 = buffer.data(gg_s + 93);
    const auto *gg_s_94 = buffer.data(gg_s + 94);
    const auto *gg_s_95 = buffer.data(gg_s + 95);
    const auto *gg_s_96 = buffer.data(gg_s + 96);
    const auto *gg_s_97 = buffer.data(gg_s + 97);
    const auto *gg_s_98 = buffer.data(gg_s + 98);
    const auto *gg_s_99 = buffer.data(gg_s + 99);
    const auto *gg_s_100 = buffer.data(gg_s + 100);
    const auto *gg_s_101 = buffer.data(gg_s + 101);
    const auto *gg_s_102 = buffer.data(gg_s + 102);
    const auto *gg_s_103 = buffer.data(gg_s + 103);
    const auto *gg_s_104 = buffer.data(gg_s + 104);
    const auto *gg_s_105 = buffer.data(gg_s + 105);
    const auto *gg_s_106 = buffer.data(gg_s + 106);

    const auto *gd_0 = buffer.data(gd + 0);
    const auto *gd_3 = buffer.data(gd + 3);
    const auto *gd_5 = buffer.data(gd + 5);
    const auto *gd_9 = buffer.data(gd + 9);
    const auto *gd_16 = buffer.data(gd + 16);
    const auto *gd_17 = buffer.data(gd + 17);
    const auto *gd_18 = buffer.data(gd + 18);
    const auto *gd_21 = buffer.data(gd + 21);
    const auto *gd_23 = buffer.data(gd + 23);
    const auto *gd_30 = buffer.data(gd + 30);
    const auto *gd_33 = buffer.data(gd + 33);
    const auto *gd_34 = buffer.data(gd + 34);
    const auto *gd_35 = buffer.data(gd + 35);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_3 = buffer.data(gf + 3);
    const auto *gf_5 = buffer.data(gf + 5);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_11 = buffer.data(gf + 11);
    const auto *gf_13 = buffer.data(gf + 13);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_18 = buffer.data(gf + 18);
    const auto *gf_19 = buffer.data(gf + 19);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_22 = buffer.data(gf + 22);
    const auto *gf_25 = buffer.data(gf + 25);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_31 = buffer.data(gf + 31);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_38 = buffer.data(gf + 38);
    const auto *gf_39 = buffer.data(gf + 39);
    const auto *gf_42 = buffer.data(gf + 42);
    const auto *gf_46 = buffer.data(gf + 46);
    const auto *gf_47 = buffer.data(gf + 47);
    const auto *gf_48 = buffer.data(gf + 48);
    const auto *gf_49 = buffer.data(gf + 49);
    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_60 = buffer.data(gf + 60);
    const auto *gf_61 = buffer.data(gf + 61);
    const auto *gf_63 = buffer.data(gf + 63);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_68 = buffer.data(gf + 68);
    const auto *gf_69 = buffer.data(gf + 69);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, ff_0, gd_s_0, gg_s_0, gg_s_1, \
                         gg_s_2, gg_s_3, gd_0, gf_0, gf_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * ff_0[k]
                 - f_1 * gd_s_0[k]
                 + f_2 * gg_s_0[k]
                 + f_3 * gd_0[k]
                 + pb_x[k] * gf_0[k];

        t_1[k] = f_2 * gg_s_1[k]
                 + pb_y[k] * gf_0[k];

        t_2[k] = f_2 * gg_s_2[k]
                 + pb_z[k] * gf_0[k];

        t_3[k] = -f_4 * gd_s_0[k]
                 + f_2 * gg_s_3[k]
                 + f_5 * gd_0[k]
                 + pb_y[k] * gf_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, ff_6, gd_s_0, gg_s_4, gg_s_5, \
                         gg_s_6, gd_0, gf_2, gf_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gg_s_4[k]
                 + pb_y[k] * gf_2[k];

        t_5[k] = -f_4 * gd_s_0[k]
                 + f_2 * gg_s_5[k]
                 + f_5 * gd_0[k]
                 + pb_z[k] * gf_2[k];

        t_6[k] = f_0 * ff_6[k]
                 + f_2 * gg_s_6[k]
                 + pb_x[k] * gf_6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, ff_9, gg_s_7, gg_s_8, gg_s_9, gf_3, \
                         gf_5, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * gg_s_7[k]
                 + pb_z[k] * gf_3[k];

        t_8[k] = f_2 * gg_s_8[k]
                 + pb_y[k] * gf_5[k];

        t_9[k] = f_0 * ff_9[k]
                 + f_2 * gg_s_9[k]
                 + pb_x[k] * gf_9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, gd_s_3, gd_s_5, gg_s_10, gg_s_11, \
                         gg_s_12, gd_3, gd_5, gf_6, gf_8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * gd_s_3[k]
                  + f_2 * gg_s_10[k]
                  + f_3 * gd_3[k]
                  + pb_y[k] * gf_6[k];

        t_11[k] = f_2 * gg_s_11[k]
                  + pb_z[k] * gf_6[k];

        t_12[k] = -f_4 * gd_s_5[k]
                  + f_2 * gg_s_12[k]
                  + f_5 * gd_5[k]
                  + pb_y[k] * gf_8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, fg_0, gd_s_5, gg_s_13, gg_s_14, \
                         gg_s_15, gd_5, gf_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * gg_s_13[k]
                  + pb_y[k] * gf_9[k];

        t_14[k] = -f_1 * gd_s_5[k]
                  + f_2 * gg_s_14[k]
                  + f_3 * gd_5[k]
                  + pb_z[k] * gf_9[k];

        t_15[k] = pa_y[k] * fg_0[k]
                  + f_2 * gg_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, ff_0, ff_1, fg_3, gg_s_16, \
                         gg_s_17, gg_s_18, gg_s_19, gf_10, gf_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * ff_0[k]
                  + f_2 * gg_s_16[k]
                  + pb_y[k] * gf_10[k];

        t_17[k] = f_2 * gg_s_17[k]
                  + pb_z[k] * gf_10[k];

        t_18[k] = f_6 * ff_1[k]
                  + pa_y[k] * fg_3[k]
                  + f_2 * gg_s_18[k];

        t_19[k] = f_2 * gg_s_19[k]
                  + pb_z[k] * gf_11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, ff_16, fg_5, gg_s_20, gg_s_21, \
                         gg_s_22, gf_13, gf_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * fg_5[k]
                  + f_2 * gg_s_20[k];

        t_21[k] = f_3 * ff_16[k]
                  + f_2 * gg_s_21[k]
                  + pb_x[k] * gf_16[k];

        t_22[k] = f_2 * gg_s_22[k]
                  + pb_z[k] * gf_13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, dg_s_25, dg_25, ff_18, fg_9, \
                         fg_25, gg_s_23, gg_s_24, gg_s_25, gf_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_3 * ff_18[k]
                  + f_2 * gg_s_23[k]
                  + pb_x[k] * gf_18[k];

        t_24[k] = pa_y[k] * fg_9[k]
                  + f_2 * gg_s_24[k];

        t_25[k] = -f_7 * dg_s_25[k]
                  + f_6 * dg_25[k]
                  + pa_x[k] * fg_25[k]
                  + f_2 * gg_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, ff_9, gd_s_9, gg_s_26, gg_s_27, \
                         gg_s_28, gd_9, gf_16, gf_17, gf_19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * gg_s_26[k]
                  + pb_z[k] * gf_16[k];

        t_27[k] = -f_4 * gd_s_9[k]
                  + f_2 * gg_s_27[k]
                  + f_5 * gd_9[k]
                  + pb_z[k] * gf_17[k];

        t_28[k] = f_5 * ff_9[k]
                  + f_2 * gg_s_28[k]
                  + pb_y[k] * gf_19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, ff_0, fg_0, fg_14, \
                         gg_s_29, gg_s_30, gg_s_31, gg_s_32, gf_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * fg_14[k]
                  + f_2 * gg_s_29[k];

        t_30[k] = pa_z[k] * fg_0[k]
                  + f_2 * gg_s_30[k];

        t_31[k] = f_2 * gg_s_31[k]
                  + pb_y[k] * gf_20[k];

        t_32[k] = f_5 * ff_0[k]
                  + f_2 * gg_s_32[k]
                  + pb_z[k] * gf_20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, ff_2, fg_3, fg_5, fg_6, gg_s_33, \
                         gg_s_34, gg_s_35, gg_s_36, gf_22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * fg_3[k]
                  + f_2 * gg_s_33[k];

        t_34[k] = f_2 * gg_s_34[k]
                  + pb_y[k] * gf_22[k];

        t_35[k] = f_6 * ff_2[k]
                  + pa_z[k] * fg_5[k]
                  + f_2 * gg_s_35[k];

        t_36[k] = pa_z[k] * fg_6[k]
                  + f_2 * gg_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, ff_27, ff_29, gg_s_37, gg_s_38, \
                         gg_s_39, gf_25, gf_27, gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_3 * ff_27[k]
                  + f_2 * gg_s_37[k]
                  + pb_x[k] * gf_27[k];

        t_38[k] = f_2 * gg_s_38[k]
                  + pb_y[k] * gf_25[k];

        t_39[k] = f_3 * ff_29[k]
                  + f_2 * gg_s_39[k]
                  + pb_x[k] * gf_29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, fg_10, gd_s_16, gd_s_17, gg_s_40, \
                         gg_s_41, gg_s_42, gd_16, gd_17, gf_27, gf_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * fg_10[k]
                  + f_2 * gg_s_40[k];

        t_41[k] = -f_8 * gd_s_16[k]
                  + f_2 * gg_s_41[k]
                  + f_6 * gd_16[k]
                  + pb_y[k] * gf_27[k];

        t_42[k] = -f_4 * gd_s_17[k]
                  + f_2 * gg_s_42[k]
                  + f_5 * gd_17[k]
                  + pb_y[k] * gf_28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pb_y, dg_s_0, dg_s_44, dg_0, dg_44, \
                         fg_15, fg_44, gg_s_43, gg_s_44, gg_s_45, \
                         gf_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * gg_s_43[k]
                  + pb_y[k] * gf_29[k];

        t_44[k] = -f_7 * dg_s_44[k]
                  + f_6 * dg_44[k]
                  + pa_x[k] * fg_44[k]
                  + f_2 * gg_s_44[k];

        t_45[k] = -f_9 * dg_s_0[k]
                  + f_5 * dg_0[k]
                  + pa_y[k] * fg_15[k]
                  + f_2 * gg_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, pb_z, ff_10, ff_33, gd_s_21, gg_s_46, \
                         gg_s_47, gg_s_48, gd_21, gf_30, gf_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * ff_10[k]
                  + f_2 * gg_s_46[k]
                  + pb_y[k] * gf_30[k];

        t_47[k] = f_2 * gg_s_47[k]
                  + pb_z[k] * gf_30[k];

        t_48[k] = f_6 * ff_33[k]
                  - f_4 * gd_s_21[k]
                  + f_2 * gg_s_48[k]
                  + f_5 * gd_21[k]
                  + pb_x[k] * gf_33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_z, ff_36, gd_s_18, gg_s_49, gg_s_50, \
                         gg_s_51, gd_18, gf_31, gf_32, gf_36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * gg_s_49[k]
                  + pb_z[k] * gf_31[k];

        t_50[k] = -f_4 * gd_s_18[k]
                  + f_2 * gg_s_50[k]
                  + f_5 * gd_18[k]
                  + pb_z[k] * gf_32[k];

        t_51[k] = f_6 * ff_36[k]
                  + f_2 * gg_s_51[k]
                  + pb_x[k] * gf_36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_z, ff_38, ff_39, gg_s_52, gg_s_53, \
                         gg_s_54, gf_33, gf_38, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * gg_s_52[k]
                  + pb_z[k] * gf_33[k];

        t_53[k] = f_6 * ff_38[k]
                  + f_2 * gg_s_53[k]
                  + pb_x[k] * gf_38[k];

        t_54[k] = f_6 * ff_39[k]
                  + f_2 * gg_s_54[k]
                  + pb_x[k] * gf_39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_x, pb_z, dg_s_55, dg_55, fg_55, gd_s_21, \
                         gg_s_55, gg_s_56, gg_s_57, gd_21, gf_36, \
                         gf_37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_9 * dg_s_55[k]
                  + f_5 * dg_55[k]
                  + pa_x[k] * fg_55[k]
                  + f_2 * gg_s_55[k];

        t_56[k] = f_2 * gg_s_56[k]
                  + pb_z[k] * gf_36[k];

        t_57[k] = -f_4 * gd_s_21[k]
                  + f_2 * gg_s_57[k]
                  + f_5 * gd_21[k]
                  + pb_z[k] * gf_37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pb_y, pb_z, ff_19, fg_30, gd_s_23, gg_s_58, \
                         gg_s_59, gg_s_60, gd_23, gf_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_6 * ff_19[k]
                  + f_2 * gg_s_58[k]
                  + pb_y[k] * gf_39[k];

        t_59[k] = -f_1 * gd_s_23[k]
                  + f_2 * gg_s_59[k]
                  + f_3 * gd_23[k]
                  + pb_z[k] * gf_39[k];

        t_60[k] = pa_y[k] * fg_30[k]
                  + f_2 * gg_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pa_z, pb_y, ff_22, fg_16, fg_18, fg_32, \
                         gg_s_61, gg_s_62, gg_s_63, gg_s_64, gf_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_z[k] * fg_16[k]
                  + f_2 * gg_s_61[k];

        t_62[k] = pa_y[k] * fg_32[k]
                  + f_2 * gg_s_62[k];

        t_63[k] = pa_z[k] * fg_18[k]
                  + f_2 * gg_s_63[k];

        t_64[k] = f_5 * ff_22[k]
                  + f_2 * gg_s_64[k]
                  + pb_y[k] * gf_42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pa_z, pb_x, ff_47, fg_21, fg_35, gg_s_65, \
                         gg_s_66, gg_s_67, gf_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * fg_35[k]
                  + f_2 * gg_s_65[k];

        t_66[k] = pa_z[k] * fg_21[k]
                  + f_2 * gg_s_66[k];

        t_67[k] = f_6 * ff_47[k]
                  + f_2 * gg_s_67[k]
                  + pb_x[k] * gf_47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pa_z, pb_x, ff_48, fg_25, fg_39, gg_s_68, \
                         gg_s_69, gg_s_70, gf_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_6 * ff_48[k]
                  + f_2 * gg_s_68[k]
                  + pb_x[k] * gf_48[k];

        t_69[k] = pa_y[k] * fg_39[k]
                  + f_2 * gg_s_69[k];

        t_70[k] = pa_z[k] * fg_25[k]
                  + f_2 * gg_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_x, pb_y, pb_z, dg_s_72, dg_72, ff_16, ff_29, \
                         fg_72, gg_s_71, gg_s_72, gg_s_73, gf_46, \
                         gf_49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * ff_16[k]
                  + f_2 * gg_s_71[k]
                  + pb_z[k] * gf_46[k];

        t_72[k] = -f_9 * dg_s_72[k]
                  + f_5 * dg_72[k]
                  + pa_x[k] * fg_72[k]
                  + f_2 * gg_s_72[k];

        t_73[k] = f_5 * ff_29[k]
                  + f_2 * gg_s_73[k]
                  + pb_y[k] * gf_49[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_y, dg_s_0, dg_0, fg_30, fg_44, \
                         gg_s_74, gg_s_75, gg_s_76, gf_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_y[k] * fg_44[k]
                  + f_2 * gg_s_74[k];

        t_75[k] = -f_9 * dg_s_0[k]
                  + f_5 * dg_0[k]
                  + pa_z[k] * fg_30[k]
                  + f_2 * gg_s_75[k];

        t_76[k] = f_2 * gg_s_76[k]
                  + pb_y[k] * gf_50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pb_z, ff_20, gd_s_30, gg_s_77, gg_s_78, \
                         gg_s_79, gd_30, gf_50, gf_51, gf_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_6 * ff_20[k]
                  + f_2 * gg_s_77[k]
                  + pb_z[k] * gf_50[k];

        t_78[k] = -f_4 * gd_s_30[k]
                  + f_2 * gg_s_78[k]
                  + f_5 * gd_30[k]
                  + pb_y[k] * gf_51[k];

        t_79[k] = f_2 * gg_s_79[k]
                  + pb_y[k] * gf_52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, ff_55, ff_56, ff_57, gd_s_35, gg_s_80, \
                         gg_s_81, gg_s_82, gd_35, gf_55, gf_56, gf_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_6 * ff_55[k]
                  - f_4 * gd_s_35[k]
                  + f_2 * gg_s_80[k]
                  + f_5 * gd_35[k]
                  + pb_x[k] * gf_55[k];

        t_81[k] = f_6 * ff_56[k]
                  + f_2 * gg_s_81[k]
                  + pb_x[k] * gf_56[k];

        t_82[k] = f_6 * ff_57[k]
                  + f_2 * gg_s_82[k]
                  + pb_x[k] * gf_57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, pb_y, ff_59, gd_s_33, gg_s_83, gg_s_84, \
                         gg_s_85, gd_33, gf_55, gf_56, gf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * gg_s_83[k]
                  + pb_y[k] * gf_55[k];

        t_84[k] = f_6 * ff_59[k]
                  + f_2 * gg_s_84[k]
                  + pb_x[k] * gf_59[k];

        t_85[k] = -f_1 * gd_s_33[k]
                  + f_2 * gg_s_85[k]
                  + f_3 * gd_33[k]
                  + pb_y[k] * gf_56[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, gd_s_34, gd_s_35, gg_s_86, gg_s_87, gg_s_88, \
                         gd_34, gd_35, gf_57, gf_58, gf_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_8 * gd_s_34[k]
                  + f_2 * gg_s_86[k]
                  + f_6 * gd_34[k]
                  + pb_y[k] * gf_57[k];

        t_87[k] = -f_4 * gd_s_35[k]
                  + f_2 * gg_s_87[k]
                  + f_5 * gd_35[k]
                  + pb_y[k] * gf_58[k];

        t_88[k] = f_2 * gg_s_88[k]
                  + pb_y[k] * gf_59[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_x, pb_y, dg_s_89, dg_89, ff_30, ff_60, fg_89, \
                         fg_90, gg_s_89, gg_s_90, gg_s_91, gf_60 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -f_9 * dg_s_89[k]
                  + f_5 * dg_89[k]
                  + pa_x[k] * fg_89[k]
                  + f_2 * gg_s_89[k];

        t_90[k] = f_0 * ff_60[k]
                  + pa_x[k] * fg_90[k]
                  + f_2 * gg_s_90[k];

        t_91[k] = f_3 * ff_30[k]
                  + f_2 * gg_s_91[k]
                  + pb_y[k] * gf_60[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, t_95, pa_x, pb_z, ff_63, ff_65, fg_93, fg_95, \
                         gg_s_92, gg_s_93, gg_s_94, gg_s_95, gf_60, \
                         gf_61 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_2 * gg_s_92[k]
                  + pb_z[k] * gf_60[k];

        t_93[k] = f_6 * ff_63[k]
                  + pa_x[k] * fg_93[k]
                  + f_2 * gg_s_93[k];

        t_94[k] = f_2 * gg_s_94[k]
                  + pb_z[k] * gf_61[k];

        t_95[k] = f_6 * ff_65[k]
                  + pa_x[k] * fg_95[k]
                  + f_2 * gg_s_95[k];
    }

#pragma omp simd aligned(t_96, t_97, t_98, pb_x, pb_z, ff_66, ff_68, gg_s_96, gg_s_97, \
                         gg_s_98, gf_63, gf_66, gf_68 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_96[k] = f_5 * ff_66[k]
                  + f_2 * gg_s_96[k]
                  + pb_x[k] * gf_66[k];

        t_97[k] = f_2 * gg_s_97[k]
                  + pb_z[k] * gf_63[k];

        t_98[k] = f_5 * ff_68[k]
                  + f_2 * gg_s_98[k]
                  + pb_x[k] * gf_68[k];
    }

#pragma omp simd aligned(t_99, t_100, t_101, t_102, pa_x, pb_x, pb_z, ff_69, fg_100, fg_102, \
                         gg_s_99, gg_s_100, gg_s_101, gg_s_102, gf_66, \
                         gf_69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_99[k] = f_5 * ff_69[k]
                  + f_2 * gg_s_99[k]
                  + pb_x[k] * gf_69[k];

        t_100[k] = pa_x[k] * fg_100[k]
                   + f_2 * gg_s_100[k];

        t_101[k] = f_2 * gg_s_101[k]
                   + pb_z[k] * gf_66[k];

        t_102[k] = pa_x[k] * fg_102[k]
                   + f_2 * gg_s_102[k];
    }

#pragma omp simd aligned(t_103, t_104, t_105, t_106, pa_x, pa_z, fg_45, fg_46, fg_103, fg_104, \
                         gg_s_103, gg_s_104, gg_s_105, gg_s_106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = pa_x[k] * fg_103[k]
                   + f_2 * gg_s_103[k];

        t_104[k] = pa_x[k] * fg_104[k]
                   + f_2 * gg_s_104[k];

        t_105[k] = pa_z[k] * fg_45[k]
                   + f_2 * gg_s_105[k];

        t_106[k] = pa_z[k] * fg_46[k]
                   + f_2 * gg_s_106[k];
    }
}

static auto
compute_prim_gg_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t dg_s,
                                        const size_t dg, const size_t ff, const size_t fg,
                                        const size_t gd_s, const size_t gg_s, const size_t gd,
                                        const size_t gf, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.0 * beta / p;
    const auto f_8 = 2.0 * alpha / p;
    const auto f_9 = beta / p;

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
    auto *t_212 = buffer.data(target + 212);
    auto *t_213 = buffer.data(target + 213);
    auto *t_214 = buffer.data(target + 214);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dg_s_55 = buffer.data(dg_s + 55);
    const auto *dg_s_74 = buffer.data(dg_s + 74);
    const auto *dg_s_89 = buffer.data(dg_s + 89);

    const auto *dg_55 = buffer.data(dg + 55);
    const auto *dg_74 = buffer.data(dg + 74);
    const auto *dg_89 = buffer.data(dg + 89);

    const auto *ff_30 = buffer.data(ff + 30);
    const auto *ff_42 = buffer.data(ff + 42);
    const auto *ff_50 = buffer.data(ff + 50);
    const auto *ff_52 = buffer.data(ff + 52);
    const auto *ff_61 = buffer.data(ff + 61);
    const auto *ff_66 = buffer.data(ff + 66);
    const auto *ff_67 = buffer.data(ff + 67);
    const auto *ff_69 = buffer.data(ff + 69);
    const auto *ff_75 = buffer.data(ff + 75);
    const auto *ff_76 = buffer.data(ff + 76);
    const auto *ff_77 = buffer.data(ff + 77);
    const auto *ff_78 = buffer.data(ff + 78);
    const auto *ff_79 = buffer.data(ff + 79);
    const auto *ff_83 = buffer.data(ff + 83);
    const auto *ff_86 = buffer.data(ff + 86);
    const auto *ff_87 = buffer.data(ff + 87);
    const auto *ff_88 = buffer.data(ff + 88);
    const auto *ff_89 = buffer.data(ff + 89);
    const auto *ff_90 = buffer.data(ff + 90);
    const auto *ff_91 = buffer.data(ff + 91);
    const auto *ff_92 = buffer.data(ff + 92);
    const auto *ff_93 = buffer.data(ff + 93);
    const auto *ff_95 = buffer.data(ff + 95);
    const auto *ff_96 = buffer.data(ff + 96);
    const auto *ff_97 = buffer.data(ff + 97);
    const auto *ff_98 = buffer.data(ff + 98);
    const auto *ff_99 = buffer.data(ff + 99);

    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_134 = buffer.data(fg + 134);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_136 = buffer.data(fg + 136);
    const auto *fg_137 = buffer.data(fg + 137);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_139 = buffer.data(fg + 139);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *gd_s_60 = buffer.data(gd_s + 60);
    const auto *gd_s_61 = buffer.data(gd_s + 61);
    const auto *gd_s_63 = buffer.data(gd_s + 63);
    const auto *gd_s_65 = buffer.data(gd_s + 65);
    const auto *gd_s_68 = buffer.data(gd_s + 68);
    const auto *gd_s_71 = buffer.data(gd_s + 71);
    const auto *gd_s_72 = buffer.data(gd_s + 72);
    const auto *gd_s_73 = buffer.data(gd_s + 73);
    const auto *gd_s_74 = buffer.data(gd_s + 74);
    const auto *gd_s_75 = buffer.data(gd_s + 75);
    const auto *gd_s_76 = buffer.data(gd_s + 76);
    const auto *gd_s_77 = buffer.data(gd_s + 77);
    const auto *gd_s_84 = buffer.data(gd_s + 84);
    const auto *gd_s_86 = buffer.data(gd_s + 86);
    const auto *gd_s_87 = buffer.data(gd_s + 87);

    const auto *gg_s_107 = buffer.data(gg_s + 107);
    const auto *gg_s_108 = buffer.data(gg_s + 108);
    const auto *gg_s_109 = buffer.data(gg_s + 109);
    const auto *gg_s_110 = buffer.data(gg_s + 110);
    const auto *gg_s_111 = buffer.data(gg_s + 111);
    const auto *gg_s_112 = buffer.data(gg_s + 112);
    const auto *gg_s_113 = buffer.data(gg_s + 113);
    const auto *gg_s_114 = buffer.data(gg_s + 114);
    const auto *gg_s_115 = buffer.data(gg_s + 115);
    const auto *gg_s_116 = buffer.data(gg_s + 116);
    const auto *gg_s_117 = buffer.data(gg_s + 117);
    const auto *gg_s_118 = buffer.data(gg_s + 118);
    const auto *gg_s_119 = buffer.data(gg_s + 119);
    const auto *gg_s_120 = buffer.data(gg_s + 120);
    const auto *gg_s_121 = buffer.data(gg_s + 121);
    const auto *gg_s_122 = buffer.data(gg_s + 122);
    const auto *gg_s_123 = buffer.data(gg_s + 123);
    const auto *gg_s_124 = buffer.data(gg_s + 124);
    const auto *gg_s_125 = buffer.data(gg_s + 125);
    const auto *gg_s_126 = buffer.data(gg_s + 126);
    const auto *gg_s_127 = buffer.data(gg_s + 127);
    const auto *gg_s_128 = buffer.data(gg_s + 128);
    const auto *gg_s_129 = buffer.data(gg_s + 129);
    const auto *gg_s_130 = buffer.data(gg_s + 130);
    const auto *gg_s_131 = buffer.data(gg_s + 131);
    const auto *gg_s_132 = buffer.data(gg_s + 132);
    const auto *gg_s_133 = buffer.data(gg_s + 133);
    const auto *gg_s_134 = buffer.data(gg_s + 134);
    const auto *gg_s_135 = buffer.data(gg_s + 135);
    const auto *gg_s_136 = buffer.data(gg_s + 136);
    const auto *gg_s_137 = buffer.data(gg_s + 137);
    const auto *gg_s_138 = buffer.data(gg_s + 138);
    const auto *gg_s_139 = buffer.data(gg_s + 139);
    const auto *gg_s_140 = buffer.data(gg_s + 140);
    const auto *gg_s_141 = buffer.data(gg_s + 141);
    const auto *gg_s_142 = buffer.data(gg_s + 142);
    const auto *gg_s_143 = buffer.data(gg_s + 143);
    const auto *gg_s_144 = buffer.data(gg_s + 144);
    const auto *gg_s_145 = buffer.data(gg_s + 145);
    const auto *gg_s_146 = buffer.data(gg_s + 146);
    const auto *gg_s_147 = buffer.data(gg_s + 147);
    const auto *gg_s_148 = buffer.data(gg_s + 148);
    const auto *gg_s_149 = buffer.data(gg_s + 149);
    const auto *gg_s_150 = buffer.data(gg_s + 150);
    const auto *gg_s_151 = buffer.data(gg_s + 151);
    const auto *gg_s_152 = buffer.data(gg_s + 152);
    const auto *gg_s_153 = buffer.data(gg_s + 153);
    const auto *gg_s_154 = buffer.data(gg_s + 154);
    const auto *gg_s_155 = buffer.data(gg_s + 155);
    const auto *gg_s_156 = buffer.data(gg_s + 156);
    const auto *gg_s_157 = buffer.data(gg_s + 157);
    const auto *gg_s_158 = buffer.data(gg_s + 158);
    const auto *gg_s_159 = buffer.data(gg_s + 159);
    const auto *gg_s_160 = buffer.data(gg_s + 160);
    const auto *gg_s_161 = buffer.data(gg_s + 161);
    const auto *gg_s_162 = buffer.data(gg_s + 162);
    const auto *gg_s_163 = buffer.data(gg_s + 163);
    const auto *gg_s_164 = buffer.data(gg_s + 164);
    const auto *gg_s_165 = buffer.data(gg_s + 165);
    const auto *gg_s_166 = buffer.data(gg_s + 166);
    const auto *gg_s_167 = buffer.data(gg_s + 167);
    const auto *gg_s_168 = buffer.data(gg_s + 168);
    const auto *gg_s_169 = buffer.data(gg_s + 169);
    const auto *gg_s_170 = buffer.data(gg_s + 170);
    const auto *gg_s_171 = buffer.data(gg_s + 171);
    const auto *gg_s_172 = buffer.data(gg_s + 172);
    const auto *gg_s_173 = buffer.data(gg_s + 173);
    const auto *gg_s_174 = buffer.data(gg_s + 174);
    const auto *gg_s_175 = buffer.data(gg_s + 175);
    const auto *gg_s_176 = buffer.data(gg_s + 176);
    const auto *gg_s_177 = buffer.data(gg_s + 177);
    const auto *gg_s_178 = buffer.data(gg_s + 178);
    const auto *gg_s_179 = buffer.data(gg_s + 179);
    const auto *gg_s_180 = buffer.data(gg_s + 180);
    const auto *gg_s_181 = buffer.data(gg_s + 181);
    const auto *gg_s_182 = buffer.data(gg_s + 182);
    const auto *gg_s_183 = buffer.data(gg_s + 183);
    const auto *gg_s_184 = buffer.data(gg_s + 184);
    const auto *gg_s_185 = buffer.data(gg_s + 185);
    const auto *gg_s_186 = buffer.data(gg_s + 186);
    const auto *gg_s_187 = buffer.data(gg_s + 187);
    const auto *gg_s_188 = buffer.data(gg_s + 188);
    const auto *gg_s_189 = buffer.data(gg_s + 189);
    const auto *gg_s_190 = buffer.data(gg_s + 190);
    const auto *gg_s_191 = buffer.data(gg_s + 191);
    const auto *gg_s_192 = buffer.data(gg_s + 192);
    const auto *gg_s_193 = buffer.data(gg_s + 193);
    const auto *gg_s_194 = buffer.data(gg_s + 194);
    const auto *gg_s_195 = buffer.data(gg_s + 195);
    const auto *gg_s_196 = buffer.data(gg_s + 196);
    const auto *gg_s_197 = buffer.data(gg_s + 197);
    const auto *gg_s_198 = buffer.data(gg_s + 198);
    const auto *gg_s_199 = buffer.data(gg_s + 199);
    const auto *gg_s_200 = buffer.data(gg_s + 200);
    const auto *gg_s_201 = buffer.data(gg_s + 201);
    const auto *gg_s_202 = buffer.data(gg_s + 202);
    const auto *gg_s_203 = buffer.data(gg_s + 203);
    const auto *gg_s_204 = buffer.data(gg_s + 204);
    const auto *gg_s_205 = buffer.data(gg_s + 205);
    const auto *gg_s_206 = buffer.data(gg_s + 206);
    const auto *gg_s_207 = buffer.data(gg_s + 207);
    const auto *gg_s_208 = buffer.data(gg_s + 208);
    const auto *gg_s_209 = buffer.data(gg_s + 209);
    const auto *gg_s_210 = buffer.data(gg_s + 210);
    const auto *gg_s_211 = buffer.data(gg_s + 211);
    const auto *gg_s_212 = buffer.data(gg_s + 212);
    const auto *gg_s_213 = buffer.data(gg_s + 213);
    const auto *gg_s_214 = buffer.data(gg_s + 214);

    const auto *gd_60 = buffer.data(gd + 60);
    const auto *gd_61 = buffer.data(gd + 61);
    const auto *gd_63 = buffer.data(gd + 63);
    const auto *gd_65 = buffer.data(gd + 65);
    const auto *gd_68 = buffer.data(gd + 68);
    const auto *gd_71 = buffer.data(gd + 71);
    const auto *gd_72 = buffer.data(gd + 72);
    const auto *gd_73 = buffer.data(gd + 73);
    const auto *gd_74 = buffer.data(gd + 74);
    const auto *gd_75 = buffer.data(gd + 75);
    const auto *gd_76 = buffer.data(gd + 76);
    const auto *gd_77 = buffer.data(gd + 77);
    const auto *gd_84 = buffer.data(gd + 84);
    const auto *gd_86 = buffer.data(gd + 86);
    const auto *gd_87 = buffer.data(gd + 87);

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
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_115 = buffer.data(gf + 115);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_117 = buffer.data(gf + 117);
    const auto *gf_118 = buffer.data(gf + 118);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_120 = buffer.data(gf + 120);
    const auto *gf_121 = buffer.data(gf + 121);
    const auto *gf_122 = buffer.data(gf + 122);
    const auto *gf_123 = buffer.data(gf + 123);
    const auto *gf_124 = buffer.data(gf + 124);
    const auto *gf_125 = buffer.data(gf + 125);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_127 = buffer.data(gf + 127);
    const auto *gf_128 = buffer.data(gf + 128);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_137 = buffer.data(gf + 137);
    const auto *gf_138 = buffer.data(gf + 138);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);

#pragma omp simd aligned(t_107, t_108, t_109, pa_z, pb_y, pb_z, ff_30, ff_42, fg_48, gg_s_107, \
                         gg_s_108, gg_s_109, gf_70, gf_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_5 * ff_30[k]
                   + f_2 * gg_s_107[k]
                   + pb_z[k] * gf_70[k];

        t_108[k] = pa_z[k] * fg_48[k]
                   + f_2 * gg_s_108[k];

        t_109[k] = f_6 * ff_42[k]
                   + f_2 * gg_s_109[k]
                   + pb_y[k] * gf_72[k];
    }

#pragma omp simd aligned(t_110, t_111, t_112, pa_x, pa_z, pb_x, ff_75, ff_77, fg_51, fg_110, \
                         gg_s_110, gg_s_111, gg_s_112, gf_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_6 * ff_75[k]
                   + pa_x[k] * fg_110[k]
                   + f_2 * gg_s_110[k];

        t_111[k] = pa_z[k] * fg_51[k]
                   + f_2 * gg_s_111[k];

        t_112[k] = f_5 * ff_77[k]
                   + f_2 * gg_s_112[k]
                   + pb_x[k] * gf_77[k];
    }

#pragma omp simd aligned(t_113, t_114, t_115, t_116, pa_x, pb_x, ff_78, ff_79, fg_115, fg_116, \
                         gg_s_113, gg_s_114, gg_s_115, gg_s_116, gf_78, \
                         gf_79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_113[k] = f_5 * ff_78[k]
                   + f_2 * gg_s_113[k]
                   + pb_x[k] * gf_78[k];

        t_114[k] = f_5 * ff_79[k]
                   + f_2 * gg_s_114[k]
                   + pb_x[k] * gf_79[k];

        t_115[k] = pa_x[k] * fg_115[k]
                   + f_2 * gg_s_115[k];

        t_116[k] = pa_x[k] * fg_116[k]
                   + f_2 * gg_s_116[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, t_120, pa_x, pa_y, fg_75, fg_117, fg_118, \
                         fg_119, gg_s_117, gg_s_118, gg_s_119, \
                         gg_s_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = pa_x[k] * fg_117[k]
                   + f_2 * gg_s_117[k];

        t_118[k] = pa_x[k] * fg_118[k]
                   + f_2 * gg_s_118[k];

        t_119[k] = pa_x[k] * fg_119[k]
                   + f_2 * gg_s_119[k];

        t_120[k] = pa_y[k] * fg_75[k]
                   + f_2 * gg_s_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_x, pa_y, pb_y, ff_50, ff_83, fg_77, fg_123, \
                         gg_s_121, gg_s_122, gg_s_123, gf_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * ff_50[k]
                   + f_2 * gg_s_121[k]
                   + pb_y[k] * gf_80[k];

        t_122[k] = pa_y[k] * fg_77[k]
                   + f_2 * gg_s_122[k];

        t_123[k] = f_6 * ff_83[k]
                   + pa_x[k] * fg_123[k]
                   + f_2 * gg_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_x, pb_y, ff_52, ff_86, fg_80, gg_s_124, \
                         gg_s_125, gg_s_126, gf_82, gf_86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * ff_52[k]
                   + f_2 * gg_s_124[k]
                   + pb_y[k] * gf_82[k];

        t_125[k] = pa_y[k] * fg_80[k]
                   + f_2 * gg_s_125[k];

        t_126[k] = f_5 * ff_86[k]
                   + f_2 * gg_s_126[k]
                   + pb_x[k] * gf_86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_y, pb_x, ff_87, ff_88, fg_84, gg_s_127, \
                         gg_s_128, gg_s_129, gf_87, gf_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_5 * ff_87[k]
                   + f_2 * gg_s_127[k]
                   + pb_x[k] * gf_87[k];

        t_128[k] = f_5 * ff_88[k]
                   + f_2 * gg_s_128[k]
                   + pb_x[k] * gf_88[k];

        t_129[k] = pa_y[k] * fg_84[k]
                   + f_2 * gg_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, t_133, t_134, pa_x, fg_130, fg_131, fg_132, \
                         fg_133, fg_134, gg_s_130, gg_s_131, gg_s_132, gg_s_133, \
                         gg_s_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = pa_x[k] * fg_130[k]
                   + f_2 * gg_s_130[k];

        t_131[k] = pa_x[k] * fg_131[k]
                   + f_2 * gg_s_131[k];

        t_132[k] = pa_x[k] * fg_132[k]
                   + f_2 * gg_s_132[k];

        t_133[k] = pa_x[k] * fg_133[k]
                   + f_2 * gg_s_133[k];

        t_134[k] = pa_x[k] * fg_134[k]
                   + f_2 * gg_s_134[k];
    }

#pragma omp simd aligned(t_135, t_136, t_137, pa_x, pb_y, pb_z, ff_50, ff_90, fg_135, \
                         gg_s_135, gg_s_136, gg_s_137, gf_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_135[k] = f_0 * ff_90[k]
                   + pa_x[k] * fg_135[k]
                   + f_2 * gg_s_135[k];

        t_136[k] = f_2 * gg_s_136[k]
                   + pb_y[k] * gf_90[k];

        t_137[k] = f_3 * ff_50[k]
                   + f_2 * gg_s_137[k]
                   + pb_z[k] * gf_90[k];
    }

#pragma omp simd aligned(t_138, t_139, t_140, pa_x, pb_y, ff_93, ff_95, fg_138, fg_140, \
                         gg_s_138, gg_s_139, gg_s_140, gf_92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_138[k] = f_6 * ff_93[k]
                   + pa_x[k] * fg_138[k]
                   + f_2 * gg_s_138[k];

        t_139[k] = f_2 * gg_s_139[k]
                   + pb_y[k] * gf_92[k];

        t_140[k] = f_6 * ff_95[k]
                   + pa_x[k] * fg_140[k]
                   + f_2 * gg_s_140[k];
    }

#pragma omp simd aligned(t_141, t_142, t_143, pb_x, pb_y, ff_96, ff_97, gg_s_141, gg_s_142, \
                         gg_s_143, gf_95, gf_96, gf_97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_141[k] = f_5 * ff_96[k]
                   + f_2 * gg_s_141[k]
                   + pb_x[k] * gf_96[k];

        t_142[k] = f_5 * ff_97[k]
                   + f_2 * gg_s_142[k]
                   + pb_x[k] * gf_97[k];

        t_143[k] = f_2 * gg_s_143[k]
                   + pb_y[k] * gf_95[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pb_x, ff_99, fg_145, fg_146, \
                         fg_147, gg_s_144, gg_s_145, gg_s_146, gg_s_147, \
                         gf_99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = f_5 * ff_99[k]
                   + f_2 * gg_s_144[k]
                   + pb_x[k] * gf_99[k];

        t_145[k] = pa_x[k] * fg_145[k]
                   + f_2 * gg_s_145[k];

        t_146[k] = pa_x[k] * fg_146[k]
                   + f_2 * gg_s_146[k];

        t_147[k] = pa_x[k] * fg_147[k]
                   + f_2 * gg_s_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_x, pb_x, pb_y, fg_149, gd_s_60, gg_s_148, \
                         gg_s_149, gg_s_150, gd_60, gf_99, gf_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_2 * gg_s_148[k]
                   + pb_y[k] * gf_99[k];

        t_149[k] = pa_x[k] * fg_149[k]
                   + f_2 * gg_s_149[k];

        t_150[k] = -f_1 * gd_s_60[k]
                   + f_2 * gg_s_150[k]
                   + f_3 * gd_60[k]
                   + pb_x[k] * gf_100[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pb_x, pb_z, gd_s_61, gd_s_63, gg_s_151, \
                         gg_s_152, gg_s_153, gd_61, gd_63, gf_100, gf_101, \
                         gf_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = -f_8 * gd_s_61[k]
                   + f_2 * gg_s_151[k]
                   + f_6 * gd_61[k]
                   + pb_x[k] * gf_101[k];

        t_152[k] = f_2 * gg_s_152[k]
                   + pb_z[k] * gf_100[k];

        t_153[k] = -f_4 * gd_s_63[k]
                   + f_2 * gg_s_153[k]
                   + f_5 * gd_63[k]
                   + pb_x[k] * gf_103[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, t_157, pb_x, pb_z, gd_s_65, gg_s_154, gg_s_155, \
                         gg_s_156, gg_s_157, gd_65, gf_101, gf_105, gf_106, \
                         gf_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_2 * gg_s_154[k]
                   + pb_z[k] * gf_101[k];

        t_155[k] = -f_4 * gd_s_65[k]
                   + f_2 * gg_s_155[k]
                   + f_5 * gd_65[k]
                   + pb_x[k] * gf_105[k];

        t_156[k] = f_2 * gg_s_156[k]
                   + pb_x[k] * gf_106[k];

        t_157[k] = f_2 * gg_s_157[k]
                   + pb_x[k] * gf_107[k];
    }

#pragma omp simd aligned(t_158, t_159, t_160, pb_x, pb_y, ff_66, gd_s_63, gg_s_158, gg_s_159, \
                         gg_s_160, gd_63, gf_106, gf_108, gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_158[k] = f_2 * gg_s_158[k]
                   + pb_x[k] * gf_108[k];

        t_159[k] = f_2 * gg_s_159[k]
                   + pb_x[k] * gf_109[k];

        t_160[k] = f_0 * ff_66[k]
                   - f_1 * gd_s_63[k]
                   + f_2 * gg_s_160[k]
                   + f_3 * gd_63[k]
                   + pb_y[k] * gf_106[k];
    }

#pragma omp simd aligned(t_161, t_162, t_163, pb_y, pb_z, ff_69, gd_s_63, gg_s_161, gg_s_162, \
                         gg_s_163, gd_63, gf_106, gf_107, gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_161[k] = f_2 * gg_s_161[k]
                   + pb_z[k] * gf_106[k];

        t_162[k] = -f_4 * gd_s_63[k]
                   + f_2 * gg_s_162[k]
                   + f_5 * gd_63[k]
                   + pb_z[k] * gf_107[k];

        t_163[k] = f_0 * ff_69[k]
                   + f_2 * gg_s_163[k]
                   + pb_y[k] * gf_109[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, pa_z, pb_z, fg_90, fg_91, gd_s_65, gg_s_164, \
                         gg_s_165, gg_s_166, gd_65, gf_109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = -f_1 * gd_s_65[k]
                   + f_2 * gg_s_164[k]
                   + f_3 * gd_65[k]
                   + pb_z[k] * gf_109[k];

        t_165[k] = pa_z[k] * fg_90[k]
                   + f_2 * gg_s_165[k];

        t_166[k] = pa_z[k] * fg_91[k]
                   + f_2 * gg_s_166[k];
    }

#pragma omp simd aligned(t_167, t_168, t_169, pa_z, pb_x, ff_61, fg_93, fg_94, gd_s_68, \
                         gg_s_167, gg_s_168, gg_s_169, gd_68, gf_112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_167[k] = -f_8 * gd_s_68[k]
                   + f_2 * gg_s_167[k]
                   + f_6 * gd_68[k]
                   + pb_x[k] * gf_112[k];

        t_168[k] = pa_z[k] * fg_93[k]
                   + f_2 * gg_s_168[k];

        t_169[k] = f_5 * ff_61[k]
                   + pa_z[k] * fg_94[k]
                   + f_2 * gg_s_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, t_173, pb_x, gd_s_71, gg_s_170, gg_s_171, \
                         gg_s_172, gg_s_173, gd_71, gf_115, gf_116, gf_117, \
                         gf_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = -f_4 * gd_s_71[k]
                   + f_2 * gg_s_170[k]
                   + f_5 * gd_71[k]
                   + pb_x[k] * gf_115[k];

        t_171[k] = f_2 * gg_s_171[k]
                   + pb_x[k] * gf_116[k];

        t_172[k] = f_2 * gg_s_172[k]
                   + pb_x[k] * gf_117[k];

        t_173[k] = f_2 * gg_s_173[k]
                   + pb_x[k] * gf_118[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_z, pb_x, pb_z, ff_66, fg_100, gg_s_174, \
                         gg_s_175, gg_s_176, gf_116, gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_2 * gg_s_174[k]
                   + pb_x[k] * gf_119[k];

        t_175[k] = pa_z[k] * fg_100[k]
                   + f_2 * gg_s_175[k];

        t_176[k] = f_5 * ff_66[k]
                   + f_2 * gg_s_176[k]
                   + pb_z[k] * gf_116[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pa_z, pb_y, dg_s_74, dg_74, ff_67, ff_79, \
                         fg_102, fg_119, gg_s_177, gg_s_178, gg_s_179, \
                         gf_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = f_6 * ff_67[k]
                   + pa_z[k] * fg_102[k]
                   + f_2 * gg_s_177[k];

        t_178[k] = f_3 * ff_79[k]
                   + f_2 * gg_s_178[k]
                   + pb_y[k] * gf_119[k];

        t_179[k] = -f_7 * dg_s_74[k]
                   + f_6 * dg_74[k]
                   + pa_y[k] * fg_119[k]
                   + f_2 * gg_s_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pb_x, gd_s_72, gd_s_73, gd_s_74, gg_s_180, \
                         gg_s_181, gg_s_182, gd_72, gd_73, gd_74, gf_120, gf_121, \
                         gf_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_1 * gd_s_72[k]
                   + f_2 * gg_s_180[k]
                   + f_3 * gd_72[k]
                   + pb_x[k] * gf_120[k];

        t_181[k] = -f_8 * gd_s_73[k]
                   + f_2 * gg_s_181[k]
                   + f_6 * gd_73[k]
                   + pb_x[k] * gf_121[k];

        t_182[k] = -f_8 * gd_s_74[k]
                   + f_2 * gg_s_182[k]
                   + f_6 * gd_74[k]
                   + pb_x[k] * gf_122[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_x, gd_s_75, gd_s_76, gd_s_77, gg_s_183, \
                         gg_s_184, gg_s_185, gd_75, gd_76, gd_77, gf_123, gf_124, \
                         gf_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = -f_4 * gd_s_75[k]
                   + f_2 * gg_s_183[k]
                   + f_5 * gd_75[k]
                   + pb_x[k] * gf_123[k];

        t_184[k] = -f_4 * gd_s_76[k]
                   + f_2 * gg_s_184[k]
                   + f_5 * gd_76[k]
                   + pb_x[k] * gf_124[k];

        t_185[k] = -f_4 * gd_s_77[k]
                   + f_2 * gg_s_185[k]
                   + f_5 * gd_77[k]
                   + pb_x[k] * gf_125[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pb_x, gg_s_186, gg_s_187, gg_s_188, \
                         gg_s_189, gf_126, gf_127, gf_128, gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = f_2 * gg_s_186[k]
                   + pb_x[k] * gf_126[k];

        t_187[k] = f_2 * gg_s_187[k]
                   + pb_x[k] * gf_127[k];

        t_188[k] = f_2 * gg_s_188[k]
                   + pb_x[k] * gf_128[k];

        t_189[k] = f_2 * gg_s_189[k]
                   + pb_x[k] * gf_129[k];
    }

#pragma omp simd aligned(t_190, t_191, pa_z, pb_z, dg_s_55, dg_55, ff_76, fg_115, gg_s_190, \
                         gg_s_191, gf_126 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = -f_9 * dg_s_55[k]
                   + f_5 * dg_55[k]
                   + pa_z[k] * fg_115[k]
                   + f_2 * gg_s_190[k];

        t_191[k] = f_6 * ff_76[k]
                   + f_2 * gg_s_191[k]
                   + pb_z[k] * gf_126[k];
    }

#pragma omp simd aligned(t_192, t_193, pb_y, ff_88, ff_89, gd_s_77, gg_s_192, gg_s_193, gd_77, \
                         gf_128, gf_129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_6 * ff_88[k]
                   - f_4 * gd_s_77[k]
                   + f_2 * gg_s_192[k]
                   + f_5 * gd_77[k]
                   + pb_y[k] * gf_128[k];

        t_193[k] = f_6 * ff_89[k]
                   + f_2 * gg_s_193[k]
                   + pb_y[k] * gf_129[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, t_197, pa_y, dg_s_89, dg_89, ff_90, fg_134, \
                         fg_135, fg_136, fg_137, gg_s_194, gg_s_195, gg_s_196, \
                         gg_s_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -f_9 * dg_s_89[k]
                   + f_5 * dg_89[k]
                   + pa_y[k] * fg_134[k]
                   + f_2 * gg_s_194[k];

        t_195[k] = pa_y[k] * fg_135[k]
                   + f_2 * gg_s_195[k];

        t_196[k] = f_5 * ff_90[k]
                   + pa_y[k] * fg_136[k]
                   + f_2 * gg_s_196[k];

        t_197[k] = pa_y[k] * fg_137[k]
                   + f_2 * gg_s_197[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, t_201, pa_y, pb_x, ff_91, ff_92, fg_138, fg_139, \
                         fg_140, gg_s_198, gg_s_199, gg_s_200, gg_s_201, \
                         gf_136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_6 * ff_91[k]
                   + pa_y[k] * fg_138[k]
                   + f_2 * gg_s_198[k];

        t_199[k] = f_5 * ff_92[k]
                   + pa_y[k] * fg_139[k]
                   + f_2 * gg_s_199[k];

        t_200[k] = pa_y[k] * fg_140[k]
                   + f_2 * gg_s_200[k];

        t_201[k] = f_2 * gg_s_201[k]
                   + pb_x[k] * gf_136[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, t_205, pa_y, pb_x, ff_96, fg_145, gg_s_202, \
                         gg_s_203, gg_s_204, gg_s_205, gf_137, gf_138, \
                         gf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_2 * gg_s_202[k]
                   + pb_x[k] * gf_137[k];

        t_203[k] = f_2 * gg_s_203[k]
                   + pb_x[k] * gf_138[k];

        t_204[k] = f_2 * gg_s_204[k]
                   + pb_x[k] * gf_139[k];

        t_205[k] = f_0 * ff_96[k]
                   + pa_y[k] * fg_145[k]
                   + f_2 * gg_s_205[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, pa_y, pb_y, pb_z, ff_86, ff_98, ff_99, fg_147, \
                         gg_s_206, gg_s_207, gg_s_208, gf_136, gf_139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = f_3 * ff_86[k]
                   + f_2 * gg_s_206[k]
                   + pb_z[k] * gf_136[k];

        t_207[k] = f_6 * ff_98[k]
                   + pa_y[k] * fg_147[k]
                   + f_2 * gg_s_207[k];

        t_208[k] = f_5 * ff_99[k]
                   + f_2 * gg_s_208[k]
                   + pb_y[k] * gf_139[k];
    }

#pragma omp simd aligned(t_209, t_210, t_211, pa_y, pb_x, pb_y, fg_149, gd_s_84, gg_s_209, \
                         gg_s_210, gg_s_211, gd_84, gf_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_209[k] = pa_y[k] * fg_149[k]
                   + f_2 * gg_s_209[k];

        t_210[k] = -f_1 * gd_s_84[k]
                   + f_2 * gg_s_210[k]
                   + f_3 * gd_84[k]
                   + pb_x[k] * gf_140[k];

        t_211[k] = f_2 * gg_s_211[k]
                   + pb_y[k] * gf_140[k];
    }

#pragma omp simd aligned(t_212, t_213, t_214, pb_x, pb_y, gd_s_86, gd_s_87, gg_s_212, \
                         gg_s_213, gg_s_214, gd_86, gd_87, gf_142, \
                         gf_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_212[k] = -f_8 * gd_s_86[k]
                   + f_2 * gg_s_212[k]
                   + f_6 * gd_86[k]
                   + pb_x[k] * gf_142[k];

        t_213[k] = -f_4 * gd_s_87[k]
                   + f_2 * gg_s_213[k]
                   + f_5 * gd_87[k]
                   + pb_x[k] * gf_143[k];

        t_214[k] = f_2 * gg_s_214[k]
                   + pb_y[k] * gf_142[k];
    }
}

static auto
compute_prim_gg_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pb, const size_t ff, const size_t gd_s,
                                        const size_t gg_s, const size_t gd, const size_t gf,
                                        const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_8 = 2.0 * alpha / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *ff_99 = buffer.data(ff + 99);

    const auto *gd_s_87 = buffer.data(gd_s + 87);
    const auto *gd_s_88 = buffer.data(gd_s + 88);
    const auto *gd_s_89 = buffer.data(gd_s + 89);

    const auto *gg_s_215 = buffer.data(gg_s + 215);
    const auto *gg_s_216 = buffer.data(gg_s + 216);
    const auto *gg_s_217 = buffer.data(gg_s + 217);
    const auto *gg_s_218 = buffer.data(gg_s + 218);
    const auto *gg_s_219 = buffer.data(gg_s + 219);
    const auto *gg_s_220 = buffer.data(gg_s + 220);
    const auto *gg_s_221 = buffer.data(gg_s + 221);
    const auto *gg_s_222 = buffer.data(gg_s + 222);
    const auto *gg_s_223 = buffer.data(gg_s + 223);
    const auto *gg_s_224 = buffer.data(gg_s + 224);

    const auto *gd_87 = buffer.data(gd + 87);
    const auto *gd_88 = buffer.data(gd + 88);
    const auto *gd_89 = buffer.data(gd + 89);

    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

#pragma omp simd aligned(t_215, t_216, t_217, t_218, pb_x, gd_s_89, gg_s_215, gg_s_216, \
                         gg_s_217, gg_s_218, gd_89, gf_145, gf_146, gf_147, \
                         gf_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = -f_4 * gd_s_89[k]
                   + f_2 * gg_s_215[k]
                   + f_5 * gd_89[k]
                   + pb_x[k] * gf_145[k];

        t_216[k] = f_2 * gg_s_216[k]
                   + pb_x[k] * gf_146[k];

        t_217[k] = f_2 * gg_s_217[k]
                   + pb_x[k] * gf_147[k];

        t_218[k] = f_2 * gg_s_218[k]
                   + pb_x[k] * gf_148[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_y, gd_s_87, gd_s_88, gg_s_219, \
                         gg_s_220, gg_s_221, gd_87, gd_88, gf_146, gf_147, \
                         gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_2 * gg_s_219[k]
                   + pb_x[k] * gf_149[k];

        t_220[k] = -f_1 * gd_s_87[k]
                   + f_2 * gg_s_220[k]
                   + f_3 * gd_87[k]
                   + pb_y[k] * gf_146[k];

        t_221[k] = -f_8 * gd_s_88[k]
                   + f_2 * gg_s_221[k]
                   + f_6 * gd_88[k]
                   + pb_y[k] * gf_147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pb_y, pb_z, ff_99, gd_s_89, gg_s_222, gg_s_223, \
                         gg_s_224, gd_89, gf_148, gf_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -f_4 * gd_s_89[k]
                   + f_2 * gg_s_222[k]
                   + f_5 * gd_89[k]
                   + pb_y[k] * gf_148[k];

        t_223[k] = f_2 * gg_s_223[k]
                   + pb_y[k] * gf_149[k];

        t_224[k] = f_0 * ff_99[k]
                   - f_1 * gd_s_89[k]
                   + f_2 * gg_s_224[k]
                   + f_3 * gd_89[k]
                   + pb_z[k] * gf_149[k];
    }
}

auto
compute_prim_gg_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dg_s, const size_t dg,
                                 const size_t ff, const size_t fg, const size_t gd_s,
                                 const size_t gg_s, const size_t gd, const size_t gf,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_gg_kinetic_energy_0_piece0(buffer, target, pa, pb, dg_s, dg, ff, fg, gd_s,
                                            gg_s, gd, gf, ncols, alpha, beta, p);

    compute_prim_gg_kinetic_energy_0_piece1(buffer, target, pa, pb, dg_s, dg, ff, fg, gd_s,
                                            gg_s, gd, gf, ncols, alpha, beta, p);

    compute_prim_gg_kinetic_energy_0_piece2(buffer, target, pb, ff, gd_s, gg_s, gd, gf, ncols,
                                            alpha, beta, p);
}

}  // namespace simdkin
