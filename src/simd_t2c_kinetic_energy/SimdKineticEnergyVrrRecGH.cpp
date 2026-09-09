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


#include "SimdKineticEnergyVrrRecGH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_gh_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t dh_s,
                                        const size_t dh, const size_t fg, const size_t fh,
                                        const size_t gf_s, const size_t gh_s, const size_t gf,
                                        const size_t gg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;

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

    const auto *dh_s_0 = buffer.data(dh_s + 0);
    const auto *dh_s_36 = buffer.data(dh_s + 36);
    const auto *dh_s_62 = buffer.data(dh_s + 62);
    const auto *dh_s_78 = buffer.data(dh_s + 78);
    const auto *dh_s_101 = buffer.data(dh_s + 101);
    const auto *dh_s_102 = buffer.data(dh_s + 102);

    const auto *dh_0 = buffer.data(dh + 0);
    const auto *dh_36 = buffer.data(dh + 36);
    const auto *dh_62 = buffer.data(dh + 62);
    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_101 = buffer.data(dh + 101);
    const auto *dh_102 = buffer.data(dh + 102);

    const auto *fg_0 = buffer.data(fg + 0);
    const auto *fg_1 = buffer.data(fg + 1);
    const auto *fg_2 = buffer.data(fg + 2);
    const auto *fg_3 = buffer.data(fg + 3);
    const auto *fg_5 = buffer.data(fg + 5);
    const auto *fg_10 = buffer.data(fg + 10);
    const auto *fg_12 = buffer.data(fg + 12);
    const auto *fg_14 = buffer.data(fg + 14);
    const auto *fg_15 = buffer.data(fg + 15);
    const auto *fg_18 = buffer.data(fg + 18);
    const auto *fg_20 = buffer.data(fg + 20);
    const auto *fg_25 = buffer.data(fg + 25);
    const auto *fg_27 = buffer.data(fg + 27);
    const auto *fg_28 = buffer.data(fg + 28);
    const auto *fg_29 = buffer.data(fg + 29);
    const auto *fg_32 = buffer.data(fg + 32);
    const auto *fg_35 = buffer.data(fg + 35);
    const auto *fg_41 = buffer.data(fg + 41);
    const auto *fg_42 = buffer.data(fg + 42);
    const auto *fg_44 = buffer.data(fg + 44);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_51 = buffer.data(fg + 51);
    const auto *fg_55 = buffer.data(fg + 55);
    const auto *fg_57 = buffer.data(fg + 57);
    const auto *fg_58 = buffer.data(fg + 58);
    const auto *fg_59 = buffer.data(fg + 59);
    const auto *fg_71 = buffer.data(fg + 71);
    const auto *fg_72 = buffer.data(fg + 72);
    const auto *fg_73 = buffer.data(fg + 73);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_3 = buffer.data(fh + 3);
    const auto *fh_5 = buffer.data(fh + 5);
    const auto *fh_6 = buffer.data(fh + 6);
    const auto *fh_7 = buffer.data(fh + 7);
    const auto *fh_9 = buffer.data(fh + 9);
    const auto *fh_10 = buffer.data(fh + 10);
    const auto *fh_14 = buffer.data(fh + 14);
    const auto *fh_15 = buffer.data(fh + 15);
    const auto *fh_20 = buffer.data(fh + 20);
    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_22 = buffer.data(fh + 22);
    const auto *fh_24 = buffer.data(fh + 24);
    const auto *fh_27 = buffer.data(fh + 27);
    const auto *fh_31 = buffer.data(fh + 31);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_44 = buffer.data(fh + 44);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_56 = buffer.data(fh + 56);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_1 = buffer.data(gf_s + 1);
    const auto *gf_s_2 = buffer.data(gf_s + 2);
    const auto *gf_s_6 = buffer.data(gf_s + 6);
    const auto *gf_s_8 = buffer.data(gf_s + 8);
    const auto *gf_s_9 = buffer.data(gf_s + 9);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_17 = buffer.data(gf_s + 17);
    const auto *gf_s_27 = buffer.data(gf_s + 27);
    const auto *gf_s_28 = buffer.data(gf_s + 28);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_32 = buffer.data(gf_s + 32);
    const auto *gf_s_33 = buffer.data(gf_s + 33);
    const auto *gf_s_36 = buffer.data(gf_s + 36);
    const auto *gf_s_37 = buffer.data(gf_s + 37);
    const auto *gf_s_39 = buffer.data(gf_s + 39);

    const auto *gh_s_0 = buffer.data(gh_s + 0);
    const auto *gh_s_1 = buffer.data(gh_s + 1);
    const auto *gh_s_2 = buffer.data(gh_s + 2);
    const auto *gh_s_3 = buffer.data(gh_s + 3);
    const auto *gh_s_4 = buffer.data(gh_s + 4);
    const auto *gh_s_5 = buffer.data(gh_s + 5);
    const auto *gh_s_6 = buffer.data(gh_s + 6);
    const auto *gh_s_7 = buffer.data(gh_s + 7);
    const auto *gh_s_8 = buffer.data(gh_s + 8);
    const auto *gh_s_9 = buffer.data(gh_s + 9);
    const auto *gh_s_10 = buffer.data(gh_s + 10);
    const auto *gh_s_11 = buffer.data(gh_s + 11);
    const auto *gh_s_12 = buffer.data(gh_s + 12);
    const auto *gh_s_13 = buffer.data(gh_s + 13);
    const auto *gh_s_14 = buffer.data(gh_s + 14);
    const auto *gh_s_15 = buffer.data(gh_s + 15);
    const auto *gh_s_16 = buffer.data(gh_s + 16);
    const auto *gh_s_17 = buffer.data(gh_s + 17);
    const auto *gh_s_18 = buffer.data(gh_s + 18);
    const auto *gh_s_19 = buffer.data(gh_s + 19);
    const auto *gh_s_20 = buffer.data(gh_s + 20);
    const auto *gh_s_21 = buffer.data(gh_s + 21);
    const auto *gh_s_22 = buffer.data(gh_s + 22);
    const auto *gh_s_23 = buffer.data(gh_s + 23);
    const auto *gh_s_24 = buffer.data(gh_s + 24);
    const auto *gh_s_25 = buffer.data(gh_s + 25);
    const auto *gh_s_26 = buffer.data(gh_s + 26);
    const auto *gh_s_27 = buffer.data(gh_s + 27);
    const auto *gh_s_28 = buffer.data(gh_s + 28);
    const auto *gh_s_29 = buffer.data(gh_s + 29);
    const auto *gh_s_30 = buffer.data(gh_s + 30);
    const auto *gh_s_31 = buffer.data(gh_s + 31);
    const auto *gh_s_32 = buffer.data(gh_s + 32);
    const auto *gh_s_33 = buffer.data(gh_s + 33);
    const auto *gh_s_34 = buffer.data(gh_s + 34);
    const auto *gh_s_35 = buffer.data(gh_s + 35);
    const auto *gh_s_36 = buffer.data(gh_s + 36);
    const auto *gh_s_37 = buffer.data(gh_s + 37);
    const auto *gh_s_38 = buffer.data(gh_s + 38);
    const auto *gh_s_39 = buffer.data(gh_s + 39);
    const auto *gh_s_40 = buffer.data(gh_s + 40);
    const auto *gh_s_41 = buffer.data(gh_s + 41);
    const auto *gh_s_42 = buffer.data(gh_s + 42);
    const auto *gh_s_43 = buffer.data(gh_s + 43);
    const auto *gh_s_44 = buffer.data(gh_s + 44);
    const auto *gh_s_45 = buffer.data(gh_s + 45);
    const auto *gh_s_46 = buffer.data(gh_s + 46);
    const auto *gh_s_47 = buffer.data(gh_s + 47);
    const auto *gh_s_48 = buffer.data(gh_s + 48);
    const auto *gh_s_49 = buffer.data(gh_s + 49);
    const auto *gh_s_50 = buffer.data(gh_s + 50);
    const auto *gh_s_51 = buffer.data(gh_s + 51);
    const auto *gh_s_52 = buffer.data(gh_s + 52);
    const auto *gh_s_53 = buffer.data(gh_s + 53);
    const auto *gh_s_54 = buffer.data(gh_s + 54);
    const auto *gh_s_55 = buffer.data(gh_s + 55);
    const auto *gh_s_56 = buffer.data(gh_s + 56);
    const auto *gh_s_57 = buffer.data(gh_s + 57);
    const auto *gh_s_58 = buffer.data(gh_s + 58);
    const auto *gh_s_59 = buffer.data(gh_s + 59);
    const auto *gh_s_60 = buffer.data(gh_s + 60);
    const auto *gh_s_61 = buffer.data(gh_s + 61);
    const auto *gh_s_62 = buffer.data(gh_s + 62);
    const auto *gh_s_63 = buffer.data(gh_s + 63);
    const auto *gh_s_64 = buffer.data(gh_s + 64);
    const auto *gh_s_65 = buffer.data(gh_s + 65);
    const auto *gh_s_66 = buffer.data(gh_s + 66);
    const auto *gh_s_67 = buffer.data(gh_s + 67);
    const auto *gh_s_68 = buffer.data(gh_s + 68);
    const auto *gh_s_69 = buffer.data(gh_s + 69);
    const auto *gh_s_70 = buffer.data(gh_s + 70);
    const auto *gh_s_71 = buffer.data(gh_s + 71);
    const auto *gh_s_72 = buffer.data(gh_s + 72);
    const auto *gh_s_73 = buffer.data(gh_s + 73);
    const auto *gh_s_74 = buffer.data(gh_s + 74);
    const auto *gh_s_75 = buffer.data(gh_s + 75);
    const auto *gh_s_76 = buffer.data(gh_s + 76);
    const auto *gh_s_77 = buffer.data(gh_s + 77);
    const auto *gh_s_78 = buffer.data(gh_s + 78);
    const auto *gh_s_79 = buffer.data(gh_s + 79);
    const auto *gh_s_80 = buffer.data(gh_s + 80);
    const auto *gh_s_81 = buffer.data(gh_s + 81);
    const auto *gh_s_82 = buffer.data(gh_s + 82);
    const auto *gh_s_83 = buffer.data(gh_s + 83);
    const auto *gh_s_84 = buffer.data(gh_s + 84);
    const auto *gh_s_85 = buffer.data(gh_s + 85);
    const auto *gh_s_86 = buffer.data(gh_s + 86);
    const auto *gh_s_87 = buffer.data(gh_s + 87);
    const auto *gh_s_88 = buffer.data(gh_s + 88);
    const auto *gh_s_89 = buffer.data(gh_s + 89);
    const auto *gh_s_90 = buffer.data(gh_s + 90);
    const auto *gh_s_91 = buffer.data(gh_s + 91);
    const auto *gh_s_92 = buffer.data(gh_s + 92);
    const auto *gh_s_93 = buffer.data(gh_s + 93);
    const auto *gh_s_94 = buffer.data(gh_s + 94);
    const auto *gh_s_95 = buffer.data(gh_s + 95);
    const auto *gh_s_96 = buffer.data(gh_s + 96);
    const auto *gh_s_97 = buffer.data(gh_s + 97);
    const auto *gh_s_98 = buffer.data(gh_s + 98);
    const auto *gh_s_99 = buffer.data(gh_s + 99);
    const auto *gh_s_100 = buffer.data(gh_s + 100);
    const auto *gh_s_101 = buffer.data(gh_s + 101);
    const auto *gh_s_102 = buffer.data(gh_s + 102);
    const auto *gh_s_103 = buffer.data(gh_s + 103);
    const auto *gh_s_104 = buffer.data(gh_s + 104);
    const auto *gh_s_105 = buffer.data(gh_s + 105);
    const auto *gh_s_106 = buffer.data(gh_s + 106);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_1 = buffer.data(gf + 1);
    const auto *gf_2 = buffer.data(gf + 2);
    const auto *gf_6 = buffer.data(gf + 6);
    const auto *gf_8 = buffer.data(gf + 8);
    const auto *gf_9 = buffer.data(gf + 9);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_17 = buffer.data(gf + 17);
    const auto *gf_27 = buffer.data(gf + 27);
    const auto *gf_28 = buffer.data(gf + 28);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_32 = buffer.data(gf + 32);
    const auto *gf_33 = buffer.data(gf + 33);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_37 = buffer.data(gf + 37);
    const auto *gf_39 = buffer.data(gf + 39);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_6 = buffer.data(gg + 6);
    const auto *gg_9 = buffer.data(gg + 9);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_13 = buffer.data(gg + 13);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_16 = buffer.data(gg + 16);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_21 = buffer.data(gg + 21);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_26 = buffer.data(gg + 26);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_39 = buffer.data(gg + 39);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_43 = buffer.data(gg + 43);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_46 = buffer.data(gg + 46);
    const auto *gg_47 = buffer.data(gg + 47);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_56 = buffer.data(gg + 56);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, fg_0, gf_s_0, gh_s_0, gh_s_1, \
                         gh_s_2, gh_s_3, gf_0, gg_0, gg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * fg_0[k]
                 - f_1 * gf_s_0[k]
                 + f_2 * gh_s_0[k]
                 + f_0 * gf_0[k]
                 + pb_x[k] * gg_0[k];

        t_1[k] = f_2 * gh_s_1[k]
                 + pb_y[k] * gg_0[k];

        t_2[k] = f_2 * gh_s_2[k]
                 + pb_z[k] * gg_0[k];

        t_3[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_3[k]
                 + f_4 * gf_0[k]
                 + pb_y[k] * gg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, gf_s_0, gf_s_1, gh_s_4, gh_s_5, \
                         gh_s_6, gh_s_7, gf_0, gf_1, gg_2, gg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * gh_s_4[k]
                 + pb_y[k] * gg_2[k];

        t_5[k] = -f_3 * gf_s_0[k]
                 + f_2 * gh_s_5[k]
                 + f_4 * gf_0[k]
                 + pb_z[k] * gg_2[k];

        t_6[k] = -f_5 * gf_s_1[k]
                 + f_2 * gh_s_6[k]
                 + f_6 * gf_1[k]
                 + pb_y[k] * gg_3[k];

        t_7[k] = f_2 * gh_s_7[k]
                 + pb_z[k] * gg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, fg_10, gf_s_2, gh_s_8, gh_s_9, \
                         gh_s_10, gf_2, gg_5, gg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * gh_s_8[k]
                 + pb_y[k] * gg_5[k];

        t_9[k] = -f_5 * gf_s_2[k]
                 + f_2 * gh_s_9[k]
                 + f_6 * gf_2[k]
                 + pb_z[k] * gg_5[k];

        t_10[k] = f_0 * fg_10[k]
                  + f_2 * gh_s_10[k]
                  + pb_x[k] * gg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, fg_12, gh_s_11, gh_s_12, gh_s_13, \
                         gg_6, gg_9, gg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * gh_s_11[k]
                  + pb_z[k] * gg_6[k];

        t_12[k] = f_0 * fg_12[k]
                  + f_2 * gh_s_12[k]
                  + pb_x[k] * gg_12[k];

        t_13[k] = f_2 * gh_s_13[k]
                  + pb_y[k] * gg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, fg_14, gf_s_6, gh_s_14, gh_s_15, \
                         gh_s_16, gf_6, gg_10, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * fg_14[k]
                  + f_2 * gh_s_14[k]
                  + pb_x[k] * gg_14[k];

        t_15[k] = -f_1 * gf_s_6[k]
                  + f_2 * gh_s_15[k]
                  + f_0 * gf_6[k]
                  + pb_y[k] * gg_10[k];

        t_16[k] = f_2 * gh_s_16[k]
                  + pb_z[k] * gg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, gf_s_8, gf_s_9, gh_s_17, gh_s_18, gh_s_19, \
                         gf_8, gf_9, gg_12, gg_13, gg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_5 * gf_s_8[k]
                  + f_2 * gh_s_17[k]
                  + f_6 * gf_8[k]
                  + pb_y[k] * gg_12[k];

        t_18[k] = -f_3 * gf_s_9[k]
                  + f_2 * gh_s_18[k]
                  + f_4 * gf_9[k]
                  + pb_y[k] * gg_13[k];

        t_19[k] = f_2 * gh_s_19[k]
                  + pb_y[k] * gg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, fg_0, fh_0, gf_s_9, gh_s_20, \
                         gh_s_21, gh_s_22, gf_9, gg_14, gg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * gf_s_9[k]
                  + f_2 * gh_s_20[k]
                  + f_0 * gf_9[k]
                  + pb_z[k] * gg_14[k];

        t_21[k] = pa_y[k] * fh_0[k]
                  + f_2 * gh_s_21[k];

        t_22[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_22[k]
                  + pb_y[k] * gg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_z, fg_1, fh_3, fh_5, gh_s_23, \
                         gh_s_24, gh_s_25, gh_s_26, gg_15, gg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * gh_s_23[k]
                  + pb_z[k] * gg_15[k];

        t_24[k] = f_6 * fg_1[k]
                  + pa_y[k] * fh_3[k]
                  + f_2 * gh_s_24[k];

        t_25[k] = f_2 * gh_s_25[k]
                  + pb_z[k] * gg_16[k];

        t_26[k] = pa_y[k] * fh_5[k]
                  + f_2 * gh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, fg_3, fg_5, fh_6, gh_s_27, \
                         gh_s_28, gh_s_29, gg_18, gg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_7 * fg_3[k]
                  + pa_y[k] * fh_6[k]
                  + f_2 * gh_s_27[k];

        t_28[k] = f_2 * gh_s_28[k]
                  + pb_z[k] * gg_18[k];

        t_29[k] = f_4 * fg_5[k]
                  + f_2 * gh_s_29[k]
                  + pb_y[k] * gg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, fg_25, fh_9, gh_s_30, gh_s_31, \
                         gh_s_32, gg_21, gg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * fh_9[k]
                  + f_2 * gh_s_30[k];

        t_31[k] = f_7 * fg_25[k]
                  + f_2 * gh_s_31[k]
                  + pb_x[k] * gg_25[k];

        t_32[k] = f_2 * gh_s_32[k]
                  + pb_z[k] * gg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, fg_27, fg_28, fh_14, gh_s_33, gh_s_34, \
                         gh_s_35, gg_27, gg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * fg_27[k]
                  + f_2 * gh_s_33[k]
                  + pb_x[k] * gg_27[k];

        t_34[k] = f_7 * fg_28[k]
                  + f_2 * gh_s_34[k]
                  + pb_x[k] * gg_28[k];

        t_35[k] = pa_y[k] * fh_14[k]
                  + f_2 * gh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_z, dh_s_36, dh_36, fh_36, gf_s_16, \
                         gh_s_36, gh_s_37, gh_s_38, gf_16, gg_25, \
                         gg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_8 * dh_s_36[k]
                  + f_6 * dh_36[k]
                  + pa_x[k] * fh_36[k]
                  + f_2 * gh_s_36[k];

        t_37[k] = f_2 * gh_s_37[k]
                  + pb_z[k] * gg_25[k];

        t_38[k] = -f_3 * gf_s_16[k]
                  + f_2 * gh_s_38[k]
                  + f_4 * gf_16[k]
                  + pb_z[k] * gg_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, pb_z, fg_14, fh_20, gf_s_17, gh_s_39, \
                         gh_s_40, gh_s_41, gf_17, gg_27, gg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_5 * gf_s_17[k]
                  + f_2 * gh_s_39[k]
                  + f_6 * gf_17[k]
                  + pb_z[k] * gg_27[k];

        t_40[k] = f_4 * fg_14[k]
                  + f_2 * gh_s_40[k]
                  + pb_y[k] * gg_29[k];

        t_41[k] = pa_y[k] * fh_20[k]
                  + f_2 * gh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, fg_0, fh_0, fh_3, gh_s_42, \
                         gh_s_43, gh_s_44, gh_s_45, gg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * fh_0[k]
                  + f_2 * gh_s_42[k];

        t_43[k] = f_2 * gh_s_43[k]
                  + pb_y[k] * gg_30[k];

        t_44[k] = f_4 * fg_0[k]
                  + f_2 * gh_s_44[k]
                  + pb_z[k] * gg_30[k];

        t_45[k] = pa_z[k] * fh_3[k]
                  + f_2 * gh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, fg_2, fg_3, fh_5, fh_6, fh_7, \
                         gh_s_46, gh_s_47, gh_s_48, gh_s_49, gg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * gh_s_46[k]
                  + pb_y[k] * gg_32[k];

        t_47[k] = f_6 * fg_2[k]
                  + pa_z[k] * fh_5[k]
                  + f_2 * gh_s_47[k];

        t_48[k] = pa_z[k] * fh_6[k]
                  + f_2 * gh_s_48[k];

        t_49[k] = f_4 * fg_3[k]
                  + pa_z[k] * fh_7[k]
                  + f_2 * gh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_y, fg_5, fh_9, fh_10, gh_s_50, gh_s_51, \
                         gh_s_52, gg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * gh_s_50[k]
                  + pb_y[k] * gg_35[k];

        t_51[k] = f_7 * fg_5[k]
                  + pa_z[k] * fh_9[k]
                  + f_2 * gh_s_51[k];

        t_52[k] = pa_z[k] * fh_10[k]
                  + f_2 * gh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, fg_41, fg_42, gh_s_53, gh_s_54, \
                         gh_s_55, gg_39, gg_41, gg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_7 * fg_41[k]
                  + f_2 * gh_s_53[k]
                  + pb_x[k] * gg_41[k];

        t_54[k] = f_7 * fg_42[k]
                  + f_2 * gh_s_54[k]
                  + pb_x[k] * gg_42[k];

        t_55[k] = f_2 * gh_s_55[k]
                  + pb_y[k] * gg_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, pb_y, fg_44, fh_15, gf_s_27, gh_s_56, \
                         gh_s_57, gh_s_58, gf_27, gg_41, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_7 * fg_44[k]
                  + f_2 * gh_s_56[k]
                  + pb_x[k] * gg_44[k];

        t_57[k] = pa_z[k] * fh_15[k]
                  + f_2 * gh_s_57[k];

        t_58[k] = -f_9 * gf_s_27[k]
                  + f_2 * gh_s_58[k]
                  + f_7 * gf_27[k]
                  + pb_y[k] * gg_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, gf_s_28, gf_s_29, gh_s_59, gh_s_60, gh_s_61, \
                         gf_28, gf_29, gg_42, gg_43, gg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_5 * gf_s_28[k]
                  + f_2 * gh_s_59[k]
                  + f_6 * gf_28[k]
                  + pb_y[k] * gg_42[k];

        t_60[k] = -f_3 * gf_s_29[k]
                  + f_2 * gh_s_60[k]
                  + f_4 * gf_29[k]
                  + pb_y[k] * gg_43[k];

        t_61[k] = f_2 * gh_s_61[k]
                  + pb_y[k] * gg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, pa_x, pa_y, dh_s_0, dh_s_62, dh_0, dh_62, fh_21, fh_62, \
                         gh_s_62, gh_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_8 * dh_s_62[k]
                  + f_6 * dh_62[k]
                  + pa_x[k] * fh_62[k]
                  + f_2 * gh_s_62[k];

        t_63[k] = -f_10 * dh_s_0[k]
                  + f_4 * dh_0[k]
                  + pa_y[k] * fh_21[k]
                  + f_2 * gh_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pb_z, fg_15, fg_48, gf_s_33, gh_s_64, \
                         gh_s_65, gh_s_66, gf_33, gg_45, gg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_6 * fg_15[k]
                  + f_2 * gh_s_64[k]
                  + pb_y[k] * gg_45[k];

        t_65[k] = f_2 * gh_s_65[k]
                  + pb_z[k] * gg_45[k];

        t_66[k] = f_6 * fg_48[k]
                  - f_5 * gf_s_33[k]
                  + f_2 * gh_s_66[k]
                  + f_6 * gf_33[k]
                  + pb_x[k] * gg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_z, fg_51, gf_s_30, gf_s_36, gh_s_67, \
                         gh_s_68, gh_s_69, gf_30, gf_36, gg_46, gg_47, \
                         gg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * gh_s_67[k]
                  + pb_z[k] * gg_46[k];

        t_68[k] = -f_3 * gf_s_30[k]
                  + f_2 * gh_s_68[k]
                  + f_4 * gf_30[k]
                  + pb_z[k] * gg_47[k];

        t_69[k] = f_6 * fg_51[k]
                  - f_3 * gf_s_36[k]
                  + f_2 * gh_s_69[k]
                  + f_4 * gf_36[k]
                  + pb_x[k] * gg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, pb_z, fg_20, gf_s_32, gh_s_70, gh_s_71, \
                         gh_s_72, gf_32, gg_48, gg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * gh_s_70[k]
                  + pb_z[k] * gg_48[k];

        t_71[k] = f_6 * fg_20[k]
                  + f_2 * gh_s_71[k]
                  + pb_y[k] * gg_50[k];

        t_72[k] = -f_5 * gf_s_32[k]
                  + f_2 * gh_s_72[k]
                  + f_6 * gf_32[k]
                  + pb_z[k] * gg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, pb_z, fg_55, fg_57, gh_s_73, gh_s_74, \
                         gh_s_75, gg_51, gg_55, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_6 * fg_55[k]
                  + f_2 * gh_s_73[k]
                  + pb_x[k] * gg_55[k];

        t_74[k] = f_2 * gh_s_74[k]
                  + pb_z[k] * gg_51[k];

        t_75[k] = f_6 * fg_57[k]
                  + f_2 * gh_s_75[k]
                  + pb_x[k] * gg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, pb_x, dh_s_78, dh_78, fg_58, fg_59, fh_78, \
                         gh_s_76, gh_s_77, gh_s_78, gg_58, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_6 * fg_58[k]
                  + f_2 * gh_s_76[k]
                  + pb_x[k] * gg_58[k];

        t_77[k] = f_6 * fg_59[k]
                  + f_2 * gh_s_77[k]
                  + pb_x[k] * gg_59[k];

        t_78[k] = -f_10 * dh_s_78[k]
                  + f_4 * dh_78[k]
                  + pa_x[k] * fh_78[k]
                  + f_2 * gh_s_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_z, gf_s_36, gf_s_37, gh_s_79, gh_s_80, gh_s_81, \
                         gf_36, gf_37, gg_55, gg_56, gg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_2 * gh_s_79[k]
                  + pb_z[k] * gg_55[k];

        t_80[k] = -f_3 * gf_s_36[k]
                  + f_2 * gh_s_80[k]
                  + f_4 * gf_36[k]
                  + pb_z[k] * gg_56[k];

        t_81[k] = -f_5 * gf_s_37[k]
                  + f_2 * gh_s_81[k]
                  + f_6 * gf_37[k]
                  + pb_z[k] * gg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_y, pb_z, fg_29, fh_42, gf_s_39, gh_s_82, \
                         gh_s_83, gh_s_84, gf_39, gg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_6 * fg_29[k]
                  + f_2 * gh_s_82[k]
                  + pb_y[k] * gg_59[k];

        t_83[k] = -f_1 * gf_s_39[k]
                  + f_2 * gh_s_83[k]
                  + f_0 * gf_39[k]
                  + pb_z[k] * gg_59[k];

        t_84[k] = pa_y[k] * fh_42[k]
                  + f_2 * gh_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pb_y, fg_32, fh_22, fh_24, fh_44, \
                         gh_s_85, gh_s_86, gh_s_87, gh_s_88, gg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_z[k] * fh_22[k]
                  + f_2 * gh_s_85[k];

        t_86[k] = pa_y[k] * fh_44[k]
                  + f_2 * gh_s_86[k];

        t_87[k] = pa_z[k] * fh_24[k]
                  + f_2 * gh_s_87[k];

        t_88[k] = f_4 * fg_32[k]
                  + f_2 * gh_s_88[k]
                  + pb_y[k] * gg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_y, pa_z, pb_z, fg_18, fh_27, fh_47, gh_s_89, \
                         gh_s_90, gh_s_91, gg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * fh_47[k]
                  + f_2 * gh_s_89[k];

        t_90[k] = pa_z[k] * fh_27[k]
                  + f_2 * gh_s_90[k];

        t_91[k] = f_4 * fg_18[k]
                  + f_2 * gh_s_91[k]
                  + pb_z[k] * gg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pa_z, pb_y, fg_35, fh_31, fh_51, gh_s_92, \
                         gh_s_93, gh_s_94, gg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_4 * fg_35[k]
                  + f_2 * gh_s_92[k]
                  + pb_y[k] * gg_65[k];

        t_93[k] = pa_y[k] * fh_51[k]
                  + f_2 * gh_s_93[k];

        t_94[k] = pa_z[k] * fh_31[k]
                  + f_2 * gh_s_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, fg_71, fg_72, fg_73, gh_s_95, gh_s_96, \
                         gh_s_97, gg_71, gg_72, gg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_6 * fg_71[k]
                  + f_2 * gh_s_95[k]
                  + pb_x[k] * gg_71[k];

        t_96[k] = f_6 * fg_72[k]
                  + f_2 * gh_s_96[k]
                  + pb_x[k] * gg_72[k];

        t_97[k] = f_6 * fg_73[k]
                  + f_2 * gh_s_97[k]
                  + pb_x[k] * gg_73[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_y, pa_z, pb_z, fg_25, fh_36, fh_56, gh_s_98, \
                         gh_s_99, gh_s_100, gg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pa_y[k] * fh_56[k]
                  + f_2 * gh_s_98[k];

        t_99[k] = pa_z[k] * fh_36[k]
                  + f_2 * gh_s_99[k];

        t_100[k] = f_4 * fg_25[k]
                   + f_2 * gh_s_100[k]
                   + pb_z[k] * gg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_x, pb_y, dh_s_101, dh_s_102, dh_101, dh_102, \
                         fg_44, fh_101, fh_102, gh_s_101, gh_s_102, gh_s_103, \
                         gg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -f_10 * dh_s_101[k]
                   + f_4 * dh_101[k]
                   + pa_x[k] * fh_101[k]
                   + f_2 * gh_s_101[k];

        t_102[k] = -f_10 * dh_s_102[k]
                   + f_4 * dh_102[k]
                   + pa_x[k] * fh_102[k]
                   + f_2 * gh_s_102[k];

        t_103[k] = f_4 * fg_44[k]
                   + f_2 * gh_s_103[k]
                   + pb_y[k] * gg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pa_z, pb_y, dh_s_0, dh_0, fh_42, fh_62, \
                         gh_s_104, gh_s_105, gh_s_106, gg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * fh_62[k]
                   + f_2 * gh_s_104[k];

        t_105[k] = -f_10 * dh_s_0[k]
                   + f_4 * dh_0[k]
                   + pa_z[k] * fh_42[k]
                   + f_2 * gh_s_105[k];

        t_106[k] = f_2 * gh_s_106[k]
                   + pb_y[k] * gg_75[k];
    }
}

static auto
compute_prim_gh_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t dh_s,
                                        const size_t dh, const size_t fg, const size_t fh,
                                        const size_t gf_s, const size_t gh_s, const size_t gf,
                                        const size_t gg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.5 / p;

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
    auto *t_215 = buffer.data(target + 215);
    auto *t_216 = buffer.data(target + 216);

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_125 = buffer.data(dh_s + 125);

    const auto *dh_125 = buffer.data(dh + 125);

    const auto *fg_30 = buffer.data(fg + 30);
    const auto *fg_45 = buffer.data(fg + 45);
    const auto *fg_48 = buffer.data(fg + 48);
    const auto *fg_50 = buffer.data(fg + 50);
    const auto *fg_62 = buffer.data(fg + 62);
    const auto *fg_63 = buffer.data(fg + 63);
    const auto *fg_65 = buffer.data(fg + 65);
    const auto *fg_75 = buffer.data(fg + 75);
    const auto *fg_77 = buffer.data(fg + 77);
    const auto *fg_80 = buffer.data(fg + 80);
    const auto *fg_84 = buffer.data(fg + 84);
    const auto *fg_85 = buffer.data(fg + 85);
    const auto *fg_86 = buffer.data(fg + 86);
    const auto *fg_87 = buffer.data(fg + 87);
    const auto *fg_89 = buffer.data(fg + 89);
    const auto *fg_90 = buffer.data(fg + 90);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_95 = buffer.data(fg + 95);
    const auto *fg_96 = buffer.data(fg + 96);
    const auto *fg_99 = buffer.data(fg + 99);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_103 = buffer.data(fg + 103);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_110 = buffer.data(fg + 110);
    const auto *fg_114 = buffer.data(fg + 114);
    const auto *fg_116 = buffer.data(fg + 116);
    const auto *fg_117 = buffer.data(fg + 117);
    const auto *fg_118 = buffer.data(fg + 118);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_123 = buffer.data(fg + 123);
    const auto *fg_126 = buffer.data(fg + 126);
    const auto *fg_130 = buffer.data(fg + 130);
    const auto *fg_131 = buffer.data(fg + 131);
    const auto *fg_132 = buffer.data(fg + 132);
    const auto *fg_133 = buffer.data(fg + 133);
    const auto *fg_135 = buffer.data(fg + 135);
    const auto *fg_138 = buffer.data(fg + 138);
    const auto *fg_140 = buffer.data(fg + 140);
    const auto *fg_141 = buffer.data(fg + 141);
    const auto *fg_142 = buffer.data(fg + 142);
    const auto *fg_144 = buffer.data(fg + 144);
    const auto *fg_145 = buffer.data(fg + 145);
    const auto *fg_146 = buffer.data(fg + 146);
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_63 = buffer.data(fh + 63);
    const auto *fh_64 = buffer.data(fh + 64);
    const auto *fh_66 = buffer.data(fh + 66);
    const auto *fh_69 = buffer.data(fh + 69);
    const auto *fh_73 = buffer.data(fh + 73);
    const auto *fh_105 = buffer.data(fh + 105);
    const auto *fh_107 = buffer.data(fh + 107);
    const auto *fh_110 = buffer.data(fh + 110);
    const auto *fh_114 = buffer.data(fh + 114);
    const auto *fh_119 = buffer.data(fh + 119);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_131 = buffer.data(fh + 131);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_135 = buffer.data(fh + 135);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_145 = buffer.data(fh + 145);
    const auto *fh_146 = buffer.data(fh + 146);
    const auto *fh_152 = buffer.data(fh + 152);
    const auto *fh_156 = buffer.data(fh + 156);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_163 = buffer.data(fh + 163);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_166 = buffer.data(fh + 166);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_171 = buffer.data(fh + 171);
    const auto *fh_174 = buffer.data(fh + 174);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_184 = buffer.data(fh + 184);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);
    const auto *fh_187 = buffer.data(fh + 187);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_205 = buffer.data(fh + 205);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_51 = buffer.data(gf_s + 51);
    const auto *gf_s_52 = buffer.data(gf_s + 52);
    const auto *gf_s_55 = buffer.data(gf_s + 55);
    const auto *gf_s_56 = buffer.data(gf_s + 56);
    const auto *gf_s_57 = buffer.data(gf_s + 57);
    const auto *gf_s_58 = buffer.data(gf_s + 58);
    const auto *gf_s_59 = buffer.data(gf_s + 59);
    const auto *gf_s_100 = buffer.data(gf_s + 100);
    const auto *gf_s_101 = buffer.data(gf_s + 101);
    const auto *gf_s_103 = buffer.data(gf_s + 103);
    const auto *gf_s_105 = buffer.data(gf_s + 105);
    const auto *gf_s_106 = buffer.data(gf_s + 106);

    const auto *gh_s_107 = buffer.data(gh_s + 107);
    const auto *gh_s_108 = buffer.data(gh_s + 108);
    const auto *gh_s_109 = buffer.data(gh_s + 109);
    const auto *gh_s_110 = buffer.data(gh_s + 110);
    const auto *gh_s_111 = buffer.data(gh_s + 111);
    const auto *gh_s_112 = buffer.data(gh_s + 112);
    const auto *gh_s_113 = buffer.data(gh_s + 113);
    const auto *gh_s_114 = buffer.data(gh_s + 114);
    const auto *gh_s_115 = buffer.data(gh_s + 115);
    const auto *gh_s_116 = buffer.data(gh_s + 116);
    const auto *gh_s_117 = buffer.data(gh_s + 117);
    const auto *gh_s_118 = buffer.data(gh_s + 118);
    const auto *gh_s_119 = buffer.data(gh_s + 119);
    const auto *gh_s_120 = buffer.data(gh_s + 120);
    const auto *gh_s_121 = buffer.data(gh_s + 121);
    const auto *gh_s_122 = buffer.data(gh_s + 122);
    const auto *gh_s_123 = buffer.data(gh_s + 123);
    const auto *gh_s_124 = buffer.data(gh_s + 124);
    const auto *gh_s_125 = buffer.data(gh_s + 125);
    const auto *gh_s_126 = buffer.data(gh_s + 126);
    const auto *gh_s_127 = buffer.data(gh_s + 127);
    const auto *gh_s_128 = buffer.data(gh_s + 128);
    const auto *gh_s_129 = buffer.data(gh_s + 129);
    const auto *gh_s_130 = buffer.data(gh_s + 130);
    const auto *gh_s_131 = buffer.data(gh_s + 131);
    const auto *gh_s_132 = buffer.data(gh_s + 132);
    const auto *gh_s_133 = buffer.data(gh_s + 133);
    const auto *gh_s_134 = buffer.data(gh_s + 134);
    const auto *gh_s_135 = buffer.data(gh_s + 135);
    const auto *gh_s_136 = buffer.data(gh_s + 136);
    const auto *gh_s_137 = buffer.data(gh_s + 137);
    const auto *gh_s_138 = buffer.data(gh_s + 138);
    const auto *gh_s_139 = buffer.data(gh_s + 139);
    const auto *gh_s_140 = buffer.data(gh_s + 140);
    const auto *gh_s_141 = buffer.data(gh_s + 141);
    const auto *gh_s_142 = buffer.data(gh_s + 142);
    const auto *gh_s_143 = buffer.data(gh_s + 143);
    const auto *gh_s_144 = buffer.data(gh_s + 144);
    const auto *gh_s_145 = buffer.data(gh_s + 145);
    const auto *gh_s_146 = buffer.data(gh_s + 146);
    const auto *gh_s_147 = buffer.data(gh_s + 147);
    const auto *gh_s_148 = buffer.data(gh_s + 148);
    const auto *gh_s_149 = buffer.data(gh_s + 149);
    const auto *gh_s_150 = buffer.data(gh_s + 150);
    const auto *gh_s_151 = buffer.data(gh_s + 151);
    const auto *gh_s_152 = buffer.data(gh_s + 152);
    const auto *gh_s_153 = buffer.data(gh_s + 153);
    const auto *gh_s_154 = buffer.data(gh_s + 154);
    const auto *gh_s_155 = buffer.data(gh_s + 155);
    const auto *gh_s_156 = buffer.data(gh_s + 156);
    const auto *gh_s_157 = buffer.data(gh_s + 157);
    const auto *gh_s_158 = buffer.data(gh_s + 158);
    const auto *gh_s_159 = buffer.data(gh_s + 159);
    const auto *gh_s_160 = buffer.data(gh_s + 160);
    const auto *gh_s_161 = buffer.data(gh_s + 161);
    const auto *gh_s_162 = buffer.data(gh_s + 162);
    const auto *gh_s_163 = buffer.data(gh_s + 163);
    const auto *gh_s_164 = buffer.data(gh_s + 164);
    const auto *gh_s_165 = buffer.data(gh_s + 165);
    const auto *gh_s_166 = buffer.data(gh_s + 166);
    const auto *gh_s_167 = buffer.data(gh_s + 167);
    const auto *gh_s_168 = buffer.data(gh_s + 168);
    const auto *gh_s_169 = buffer.data(gh_s + 169);
    const auto *gh_s_170 = buffer.data(gh_s + 170);
    const auto *gh_s_171 = buffer.data(gh_s + 171);
    const auto *gh_s_172 = buffer.data(gh_s + 172);
    const auto *gh_s_173 = buffer.data(gh_s + 173);
    const auto *gh_s_174 = buffer.data(gh_s + 174);
    const auto *gh_s_175 = buffer.data(gh_s + 175);
    const auto *gh_s_176 = buffer.data(gh_s + 176);
    const auto *gh_s_177 = buffer.data(gh_s + 177);
    const auto *gh_s_178 = buffer.data(gh_s + 178);
    const auto *gh_s_179 = buffer.data(gh_s + 179);
    const auto *gh_s_180 = buffer.data(gh_s + 180);
    const auto *gh_s_181 = buffer.data(gh_s + 181);
    const auto *gh_s_182 = buffer.data(gh_s + 182);
    const auto *gh_s_183 = buffer.data(gh_s + 183);
    const auto *gh_s_184 = buffer.data(gh_s + 184);
    const auto *gh_s_185 = buffer.data(gh_s + 185);
    const auto *gh_s_186 = buffer.data(gh_s + 186);
    const auto *gh_s_187 = buffer.data(gh_s + 187);
    const auto *gh_s_188 = buffer.data(gh_s + 188);
    const auto *gh_s_189 = buffer.data(gh_s + 189);
    const auto *gh_s_190 = buffer.data(gh_s + 190);
    const auto *gh_s_191 = buffer.data(gh_s + 191);
    const auto *gh_s_192 = buffer.data(gh_s + 192);
    const auto *gh_s_193 = buffer.data(gh_s + 193);
    const auto *gh_s_194 = buffer.data(gh_s + 194);
    const auto *gh_s_195 = buffer.data(gh_s + 195);
    const auto *gh_s_196 = buffer.data(gh_s + 196);
    const auto *gh_s_197 = buffer.data(gh_s + 197);
    const auto *gh_s_198 = buffer.data(gh_s + 198);
    const auto *gh_s_199 = buffer.data(gh_s + 199);
    const auto *gh_s_200 = buffer.data(gh_s + 200);
    const auto *gh_s_201 = buffer.data(gh_s + 201);
    const auto *gh_s_202 = buffer.data(gh_s + 202);
    const auto *gh_s_203 = buffer.data(gh_s + 203);
    const auto *gh_s_204 = buffer.data(gh_s + 204);
    const auto *gh_s_205 = buffer.data(gh_s + 205);
    const auto *gh_s_206 = buffer.data(gh_s + 206);
    const auto *gh_s_207 = buffer.data(gh_s + 207);
    const auto *gh_s_208 = buffer.data(gh_s + 208);
    const auto *gh_s_209 = buffer.data(gh_s + 209);
    const auto *gh_s_210 = buffer.data(gh_s + 210);
    const auto *gh_s_211 = buffer.data(gh_s + 211);
    const auto *gh_s_212 = buffer.data(gh_s + 212);
    const auto *gh_s_213 = buffer.data(gh_s + 213);
    const auto *gh_s_214 = buffer.data(gh_s + 214);
    const auto *gh_s_215 = buffer.data(gh_s + 215);
    const auto *gh_s_216 = buffer.data(gh_s + 216);

    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_51 = buffer.data(gf + 51);
    const auto *gf_52 = buffer.data(gf + 52);
    const auto *gf_55 = buffer.data(gf + 55);
    const auto *gf_56 = buffer.data(gf + 56);
    const auto *gf_57 = buffer.data(gf + 57);
    const auto *gf_58 = buffer.data(gf + 58);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_100 = buffer.data(gf + 100);
    const auto *gf_101 = buffer.data(gf + 101);
    const auto *gf_103 = buffer.data(gf + 103);
    const auto *gf_105 = buffer.data(gf + 105);
    const auto *gf_106 = buffer.data(gf + 106);

    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_79 = buffer.data(gg + 79);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_88 = buffer.data(gg + 88);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_91 = buffer.data(gg + 91);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);

#pragma omp simd aligned(t_107, t_108, t_109, pb_y, pb_z, fg_30, gf_s_50, gh_s_107, gh_s_108, \
                         gh_s_109, gf_50, gg_75, gg_76, gg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_6 * fg_30[k]
                   + f_2 * gh_s_107[k]
                   + pb_z[k] * gg_75[k];

        t_108[k] = -f_3 * gf_s_50[k]
                   + f_2 * gh_s_108[k]
                   + f_4 * gf_50[k]
                   + pb_y[k] * gg_76[k];

        t_109[k] = f_2 * gh_s_109[k]
                   + pb_y[k] * gg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, pb_x, pb_y, fg_80, gf_s_51, gf_s_55, gh_s_110, \
                         gh_s_111, gf_51, gf_55, gg_78, gg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_6 * fg_80[k]
                   - f_5 * gf_s_55[k]
                   + f_2 * gh_s_110[k]
                   + f_6 * gf_55[k]
                   + pb_x[k] * gg_80[k];

        t_111[k] = -f_5 * gf_s_51[k]
                   + f_2 * gh_s_111[k]
                   + f_6 * gf_51[k]
                   + pb_y[k] * gg_78[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, pb_y, fg_84, gf_s_52, gf_s_59, gh_s_112, \
                         gh_s_113, gh_s_114, gf_52, gf_59, gg_79, gg_80, \
                         gg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -f_3 * gf_s_52[k]
                   + f_2 * gh_s_112[k]
                   + f_4 * gf_52[k]
                   + pb_y[k] * gg_79[k];

        t_113[k] = f_2 * gh_s_113[k]
                   + pb_y[k] * gg_80[k];

        t_114[k] = f_6 * fg_84[k]
                   - f_3 * gf_s_59[k]
                   + f_2 * gh_s_114[k]
                   + f_4 * gf_59[k]
                   + pb_x[k] * gg_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, fg_85, fg_86, fg_87, gh_s_115, gh_s_116, \
                         gh_s_117, gg_85, gg_86, gg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_6 * fg_85[k]
                   + f_2 * gh_s_115[k]
                   + pb_x[k] * gg_85[k];

        t_116[k] = f_6 * fg_86[k]
                   + f_2 * gh_s_116[k]
                   + pb_x[k] * gg_86[k];

        t_117[k] = f_6 * fg_87[k]
                   + f_2 * gh_s_117[k]
                   + pb_x[k] * gg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_y, fg_89, gf_s_56, gh_s_118, gh_s_119, \
                         gh_s_120, gf_56, gg_84, gg_85, gg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * gh_s_118[k]
                   + pb_y[k] * gg_84[k];

        t_119[k] = f_6 * fg_89[k]
                   + f_2 * gh_s_119[k]
                   + pb_x[k] * gg_89[k];

        t_120[k] = -f_1 * gf_s_56[k]
                   + f_2 * gh_s_120[k]
                   + f_0 * gf_56[k]
                   + pb_y[k] * gg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, gf_s_57, gf_s_58, gf_s_59, gh_s_121, \
                         gh_s_122, gh_s_123, gf_57, gf_58, gf_59, gg_86, gg_87, \
                         gg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_9 * gf_s_57[k]
                   + f_2 * gh_s_121[k]
                   + f_7 * gf_57[k]
                   + pb_y[k] * gg_86[k];

        t_122[k] = -f_5 * gf_s_58[k]
                   + f_2 * gh_s_122[k]
                   + f_6 * gf_58[k]
                   + pb_y[k] * gg_87[k];

        t_123[k] = -f_3 * gf_s_59[k]
                   + f_2 * gh_s_123[k]
                   + f_4 * gf_59[k]
                   + pb_y[k] * gg_88[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pb_y, dh_s_125, dh_125, fg_90, fh_125, \
                         fh_126, gh_s_124, gh_s_125, gh_s_126, gg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * gh_s_124[k]
                   + pb_y[k] * gg_89[k];

        t_125[k] = -f_10 * dh_s_125[k]
                   + f_4 * dh_125[k]
                   + pa_x[k] * fh_125[k]
                   + f_2 * gh_s_125[k];

        t_126[k] = f_11 * fg_90[k]
                   + pa_x[k] * fh_126[k]
                   + f_2 * gh_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, t_130, pa_x, pb_y, pb_z, fg_45, fg_93, fh_129, \
                         gh_s_127, gh_s_128, gh_s_129, gh_s_130, gg_90, \
                         gg_91 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_7 * fg_45[k]
                   + f_2 * gh_s_127[k]
                   + pb_y[k] * gg_90[k];

        t_128[k] = f_2 * gh_s_128[k]
                   + pb_z[k] * gg_90[k];

        t_129[k] = f_7 * fg_93[k]
                   + pa_x[k] * fh_129[k]
                   + f_2 * gh_s_129[k];

        t_130[k] = f_2 * gh_s_130[k]
                   + pb_z[k] * gg_91[k];
    }

#pragma omp simd aligned(t_131, t_132, t_133, pa_x, pb_z, fg_95, fg_96, fh_131, fh_132, \
                         gh_s_131, gh_s_132, gh_s_133, gg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_131[k] = f_7 * fg_95[k]
                   + pa_x[k] * fh_131[k]
                   + f_2 * gh_s_131[k];

        t_132[k] = f_6 * fg_96[k]
                   + pa_x[k] * fh_132[k]
                   + f_2 * gh_s_132[k];

        t_133[k] = f_2 * gh_s_133[k]
                   + pb_z[k] * gg_93[k];
    }

#pragma omp simd aligned(t_134, t_135, t_136, pa_x, pb_x, pb_y, fg_50, fg_99, fg_100, fh_135, \
                         gh_s_134, gh_s_135, gh_s_136, gg_95, gg_100 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_134[k] = f_7 * fg_50[k]
                   + f_2 * gh_s_134[k]
                   + pb_y[k] * gg_95[k];

        t_135[k] = f_6 * fg_99[k]
                   + pa_x[k] * fh_135[k]
                   + f_2 * gh_s_135[k];

        t_136[k] = f_4 * fg_100[k]
                   + f_2 * gh_s_136[k]
                   + pb_x[k] * gg_100[k];
    }

#pragma omp simd aligned(t_137, t_138, t_139, pb_x, pb_z, fg_102, fg_103, gh_s_137, gh_s_138, \
                         gh_s_139, gg_96, gg_102, gg_103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_137[k] = f_2 * gh_s_137[k]
                   + pb_z[k] * gg_96[k];

        t_138[k] = f_4 * fg_102[k]
                   + f_2 * gh_s_138[k]
                   + pb_x[k] * gg_102[k];

        t_139[k] = f_4 * fg_103[k]
                   + f_2 * gh_s_139[k]
                   + pb_x[k] * gg_103[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, t_143, pa_x, pb_x, pb_z, fg_104, fh_141, fh_143, \
                         gh_s_140, gh_s_141, gh_s_142, gh_s_143, gg_100, \
                         gg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_4 * fg_104[k]
                   + f_2 * gh_s_140[k]
                   + pb_x[k] * gg_104[k];

        t_141[k] = pa_x[k] * fh_141[k]
                   + f_2 * gh_s_141[k];

        t_142[k] = f_2 * gh_s_142[k]
                   + pb_z[k] * gg_100[k];

        t_143[k] = pa_x[k] * fh_143[k]
                   + f_2 * gh_s_143[k];
    }

#pragma omp simd aligned(t_144, t_145, t_146, t_147, pa_x, pa_z, fh_63, fh_144, fh_145, \
                         fh_146, gh_s_144, gh_s_145, gh_s_146, \
                         gh_s_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_144[k] = pa_x[k] * fh_144[k]
                   + f_2 * gh_s_144[k];

        t_145[k] = pa_x[k] * fh_145[k]
                   + f_2 * gh_s_145[k];

        t_146[k] = pa_x[k] * fh_146[k]
                   + f_2 * gh_s_146[k];

        t_147[k] = pa_z[k] * fh_63[k]
                   + f_2 * gh_s_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_z, pb_z, fg_45, fh_64, fh_66, gh_s_148, \
                         gh_s_149, gh_s_150, gg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pa_z[k] * fh_64[k]
                   + f_2 * gh_s_148[k];

        t_149[k] = f_4 * fg_45[k]
                   + f_2 * gh_s_149[k]
                   + pb_z[k] * gg_105[k];

        t_150[k] = pa_z[k] * fh_66[k]
                   + f_2 * gh_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_x, pa_z, pb_y, fg_62, fg_110, fh_69, fh_152, \
                         gh_s_151, gh_s_152, gh_s_153, gg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_6 * fg_62[k]
                   + f_2 * gh_s_151[k]
                   + pb_y[k] * gg_107[k];

        t_152[k] = f_7 * fg_110[k]
                   + pa_x[k] * fh_152[k]
                   + f_2 * gh_s_152[k];

        t_153[k] = pa_z[k] * fh_69[k]
                   + f_2 * gh_s_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_x, pb_y, pb_z, fg_48, fg_65, fg_114, fh_156, \
                         gh_s_154, gh_s_155, gh_s_156, gg_108, gg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_4 * fg_48[k]
                   + f_2 * gh_s_154[k]
                   + pb_z[k] * gg_108[k];

        t_155[k] = f_6 * fg_65[k]
                   + f_2 * gh_s_155[k]
                   + pb_y[k] * gg_110[k];

        t_156[k] = f_6 * fg_114[k]
                   + pa_x[k] * fh_156[k]
                   + f_2 * gh_s_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_z, pb_x, fg_116, fg_117, fh_73, gh_s_157, \
                         gh_s_158, gh_s_159, gg_116, gg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_z[k] * fh_73[k]
                   + f_2 * gh_s_157[k];

        t_158[k] = f_4 * fg_116[k]
                   + f_2 * gh_s_158[k]
                   + pb_x[k] * gg_116[k];

        t_159[k] = f_4 * fg_117[k]
                   + f_2 * gh_s_159[k]
                   + pb_x[k] * gg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, t_163, pa_x, pb_x, fg_118, fg_119, fh_162, \
                         fh_163, gh_s_160, gh_s_161, gh_s_162, gh_s_163, gg_118, \
                         gg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_4 * fg_118[k]
                   + f_2 * gh_s_160[k]
                   + pb_x[k] * gg_118[k];

        t_161[k] = f_4 * fg_119[k]
                   + f_2 * gh_s_161[k]
                   + pb_x[k] * gg_119[k];

        t_162[k] = pa_x[k] * fh_162[k]
                   + f_2 * gh_s_162[k];

        t_163[k] = pa_x[k] * fh_163[k]
                   + f_2 * gh_s_163[k];
    }

#pragma omp simd aligned(t_164, t_165, t_166, t_167, pa_x, fh_164, fh_165, fh_166, fh_167, \
                         gh_s_164, gh_s_165, gh_s_166, gh_s_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_164[k] = pa_x[k] * fh_164[k]
                   + f_2 * gh_s_164[k];

        t_165[k] = pa_x[k] * fh_165[k]
                   + f_2 * gh_s_165[k];

        t_166[k] = pa_x[k] * fh_166[k]
                   + f_2 * gh_s_166[k];

        t_167[k] = pa_x[k] * fh_167[k]
                   + f_2 * gh_s_167[k];
    }

#pragma omp simd aligned(t_168, t_169, t_170, pa_y, pb_y, fg_75, fh_105, fh_107, gh_s_168, \
                         gh_s_169, gh_s_170, gg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_168[k] = pa_y[k] * fh_105[k]
                   + f_2 * gh_s_168[k];

        t_169[k] = f_4 * fg_75[k]
                   + f_2 * gh_s_169[k]
                   + pb_y[k] * gg_120[k];

        t_170[k] = pa_y[k] * fh_107[k]
                   + f_2 * gh_s_170[k];
    }

#pragma omp simd aligned(t_171, t_172, t_173, pa_x, pa_y, pb_y, fg_77, fg_123, fh_110, fh_171, \
                         gh_s_171, gh_s_172, gh_s_173, gg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_171[k] = f_7 * fg_123[k]
                   + pa_x[k] * fh_171[k]
                   + f_2 * gh_s_171[k];

        t_172[k] = f_4 * fg_77[k]
                   + f_2 * gh_s_172[k]
                   + pb_y[k] * gg_122[k];

        t_173[k] = pa_y[k] * fh_110[k]
                   + f_2 * gh_s_173[k];
    }

#pragma omp simd aligned(t_174, t_175, t_176, pa_x, pb_y, pb_z, fg_63, fg_80, fg_126, fh_174, \
                         gh_s_174, gh_s_175, gh_s_176, gg_123, gg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_174[k] = f_6 * fg_126[k]
                   + pa_x[k] * fh_174[k]
                   + f_2 * gh_s_174[k];

        t_175[k] = f_6 * fg_63[k]
                   + f_2 * gh_s_175[k]
                   + pb_z[k] * gg_123[k];

        t_176[k] = f_4 * fg_80[k]
                   + f_2 * gh_s_176[k]
                   + pb_y[k] * gg_125[k];
    }

#pragma omp simd aligned(t_177, t_178, t_179, pa_y, pb_x, fg_130, fg_131, fh_114, gh_s_177, \
                         gh_s_178, gh_s_179, gg_130, gg_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_177[k] = pa_y[k] * fh_114[k]
                   + f_2 * gh_s_177[k];

        t_178[k] = f_4 * fg_130[k]
                   + f_2 * gh_s_178[k]
                   + pb_x[k] * gg_130[k];

        t_179[k] = f_4 * fg_131[k]
                   + f_2 * gh_s_179[k]
                   + pb_x[k] * gg_131[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pb_x, fg_132, fg_133, fh_119, gh_s_180, \
                         gh_s_181, gh_s_182, gg_132, gg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_4 * fg_132[k]
                   + f_2 * gh_s_180[k]
                   + pb_x[k] * gg_132[k];

        t_181[k] = f_4 * fg_133[k]
                   + f_2 * gh_s_181[k]
                   + pb_x[k] * gg_133[k];

        t_182[k] = pa_y[k] * fh_119[k]
                   + f_2 * gh_s_182[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, t_186, t_187, pa_x, fh_183, fh_184, fh_185, \
                         fh_186, fh_187, gh_s_183, gh_s_184, gh_s_185, gh_s_186, \
                         gh_s_187 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = pa_x[k] * fh_183[k]
                   + f_2 * gh_s_183[k];

        t_184[k] = pa_x[k] * fh_184[k]
                   + f_2 * gh_s_184[k];

        t_185[k] = pa_x[k] * fh_185[k]
                   + f_2 * gh_s_185[k];

        t_186[k] = pa_x[k] * fh_186[k]
                   + f_2 * gh_s_186[k];

        t_187[k] = pa_x[k] * fh_187[k]
                   + f_2 * gh_s_187[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, t_191, pa_x, pb_y, pb_z, fg_75, fg_135, fh_188, \
                         fh_189, gh_s_188, gh_s_189, gh_s_190, gh_s_191, \
                         gg_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = pa_x[k] * fh_188[k]
                   + f_2 * gh_s_188[k];

        t_189[k] = f_11 * fg_135[k]
                   + pa_x[k] * fh_189[k]
                   + f_2 * gh_s_189[k];

        t_190[k] = f_2 * gh_s_190[k]
                   + pb_y[k] * gg_135[k];

        t_191[k] = f_7 * fg_75[k]
                   + f_2 * gh_s_191[k]
                   + pb_z[k] * gg_135[k];
    }

#pragma omp simd aligned(t_192, t_193, t_194, pa_x, pb_y, fg_138, fg_140, fh_192, fh_194, \
                         gh_s_192, gh_s_193, gh_s_194, gg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_192[k] = f_7 * fg_138[k]
                   + pa_x[k] * fh_192[k]
                   + f_2 * gh_s_192[k];

        t_193[k] = f_2 * gh_s_193[k]
                   + pb_y[k] * gg_137[k];

        t_194[k] = f_7 * fg_140[k]
                   + pa_x[k] * fh_194[k]
                   + f_2 * gh_s_194[k];
    }

#pragma omp simd aligned(t_195, t_196, t_197, pa_x, pb_y, fg_141, fg_142, fh_195, fh_196, \
                         gh_s_195, gh_s_196, gh_s_197, gg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_195[k] = f_6 * fg_141[k]
                   + pa_x[k] * fh_195[k]
                   + f_2 * gh_s_195[k];

        t_196[k] = f_6 * fg_142[k]
                   + pa_x[k] * fh_196[k]
                   + f_2 * gh_s_196[k];

        t_197[k] = f_2 * gh_s_197[k]
                   + pb_y[k] * gg_140[k];
    }

#pragma omp simd aligned(t_198, t_199, t_200, pa_x, pb_x, fg_144, fg_145, fg_146, fh_198, \
                         gh_s_198, gh_s_199, gh_s_200, gg_145, gg_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_198[k] = f_6 * fg_144[k]
                   + pa_x[k] * fh_198[k]
                   + f_2 * gh_s_198[k];

        t_199[k] = f_4 * fg_145[k]
                   + f_2 * gh_s_199[k]
                   + pb_x[k] * gg_145[k];

        t_200[k] = f_4 * fg_146[k]
                   + f_2 * gh_s_200[k]
                   + pb_x[k] * gg_146[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, pb_y, fg_147, fg_149, gh_s_201, gh_s_202, \
                         gh_s_203, gg_144, gg_147, gg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_4 * fg_147[k]
                   + f_2 * gh_s_201[k]
                   + pb_x[k] * gg_147[k];

        t_202[k] = f_2 * gh_s_202[k]
                   + pb_y[k] * gg_144[k];

        t_203[k] = f_4 * fg_149[k]
                   + f_2 * gh_s_203[k]
                   + pb_x[k] * gg_149[k];
    }

#pragma omp simd aligned(t_204, t_205, t_206, t_207, pa_x, fh_204, fh_205, fh_206, fh_207, \
                         gh_s_204, gh_s_205, gh_s_206, gh_s_207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_x[k] * fh_204[k]
                   + f_2 * gh_s_204[k];

        t_205[k] = pa_x[k] * fh_205[k]
                   + f_2 * gh_s_205[k];

        t_206[k] = pa_x[k] * fh_206[k]
                   + f_2 * gh_s_206[k];

        t_207[k] = pa_x[k] * fh_207[k]
                   + f_2 * gh_s_207[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_x, pb_x, pb_y, fh_209, gf_s_100, gh_s_208, \
                         gh_s_209, gh_s_210, gf_100, gg_149, gg_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_2 * gh_s_208[k]
                   + pb_y[k] * gg_149[k];

        t_209[k] = pa_x[k] * fh_209[k]
                   + f_2 * gh_s_209[k];

        t_210[k] = -f_1 * gf_s_100[k]
                   + f_2 * gh_s_210[k]
                   + f_0 * gf_100[k]
                   + pb_x[k] * gg_150[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, pb_x, pb_z, gf_s_101, gf_s_103, gh_s_211, \
                         gh_s_212, gh_s_213, gf_101, gf_103, gg_150, gg_151, \
                         gg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = -f_9 * gf_s_101[k]
                   + f_2 * gh_s_211[k]
                   + f_7 * gf_101[k]
                   + pb_x[k] * gg_151[k];

        t_212[k] = f_2 * gh_s_212[k]
                   + pb_z[k] * gg_150[k];

        t_213[k] = -f_5 * gf_s_103[k]
                   + f_2 * gh_s_213[k]
                   + f_6 * gf_103[k]
                   + pb_x[k] * gg_153[k];
    }

#pragma omp simd aligned(t_214, t_215, t_216, pb_x, pb_z, gf_s_105, gf_s_106, gh_s_214, \
                         gh_s_215, gh_s_216, gf_105, gf_106, gg_151, gg_155, \
                         gg_156 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_214[k] = f_2 * gh_s_214[k]
                   + pb_z[k] * gg_151[k];

        t_215[k] = -f_5 * gf_s_105[k]
                   + f_2 * gh_s_215[k]
                   + f_6 * gf_105[k]
                   + pb_x[k] * gg_155[k];

        t_216[k] = -f_3 * gf_s_106[k]
                   + f_2 * gh_s_216[k]
                   + f_4 * gf_106[k]
                   + pb_x[k] * gg_156[k];
    }
}

static auto
compute_prim_gh_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t dh_s,
                                        const size_t dh, const size_t fg, const size_t fh,
                                        const size_t gf_s, const size_t gh_s, const size_t gf,
                                        const size_t gg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.0 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = alpha / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.0 * alpha / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 1.5 / p;
    const auto f_8 = 2.0 * beta / p;
    const auto f_9 = 3.0 * alpha / p;
    const auto f_10 = beta / p;
    const auto f_11 = 2.5 / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *dh_s_78 = buffer.data(dh_s + 78);
    const auto *dh_s_104 = buffer.data(dh_s + 104);
    const auto *dh_s_125 = buffer.data(dh_s + 125);

    const auto *dh_78 = buffer.data(dh + 78);
    const auto *dh_104 = buffer.data(dh + 104);
    const auto *dh_125 = buffer.data(dh + 125);

    const auto *fg_91 = buffer.data(fg + 91);
    const auto *fg_93 = buffer.data(fg + 93);
    const auto *fg_94 = buffer.data(fg + 94);
    const auto *fg_100 = buffer.data(fg + 100);
    const auto *fg_101 = buffer.data(fg + 101);
    const auto *fg_102 = buffer.data(fg + 102);
    const auto *fg_104 = buffer.data(fg + 104);
    const auto *fg_115 = buffer.data(fg + 115);
    const auto *fg_119 = buffer.data(fg + 119);
    const auto *fg_130 = buffer.data(fg + 130);
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
    const auto *fg_147 = buffer.data(fg + 147);
    const auto *fg_148 = buffer.data(fg + 148);
    const auto *fg_149 = buffer.data(fg + 149);

    const auto *fh_126 = buffer.data(fh + 126);
    const auto *fh_127 = buffer.data(fh + 127);
    const auto *fh_129 = buffer.data(fh + 129);
    const auto *fh_130 = buffer.data(fh + 130);
    const auto *fh_132 = buffer.data(fh + 132);
    const auto *fh_133 = buffer.data(fh + 133);
    const auto *fh_134 = buffer.data(fh + 134);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_143 = buffer.data(fh + 143);
    const auto *fh_144 = buffer.data(fh + 144);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_189 = buffer.data(fh + 189);
    const auto *fh_190 = buffer.data(fh + 190);
    const auto *fh_191 = buffer.data(fh + 191);
    const auto *fh_192 = buffer.data(fh + 192);
    const auto *fh_193 = buffer.data(fh + 193);
    const auto *fh_194 = buffer.data(fh + 194);
    const auto *fh_195 = buffer.data(fh + 195);
    const auto *fh_196 = buffer.data(fh + 196);
    const auto *fh_197 = buffer.data(fh + 197);
    const auto *fh_198 = buffer.data(fh + 198);
    const auto *fh_204 = buffer.data(fh + 204);
    const auto *fh_206 = buffer.data(fh + 206);
    const auto *fh_207 = buffer.data(fh + 207);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gf_s_106 = buffer.data(gf_s + 106);
    const auto *gf_s_107 = buffer.data(gf_s + 107);
    const auto *gf_s_108 = buffer.data(gf_s + 108);
    const auto *gf_s_109 = buffer.data(gf_s + 109);
    const auto *gf_s_112 = buffer.data(gf_s + 112);
    const auto *gf_s_115 = buffer.data(gf_s + 115);
    const auto *gf_s_119 = buffer.data(gf_s + 119);
    const auto *gf_s_120 = buffer.data(gf_s + 120);
    const auto *gf_s_121 = buffer.data(gf_s + 121);
    const auto *gf_s_122 = buffer.data(gf_s + 122);
    const auto *gf_s_123 = buffer.data(gf_s + 123);
    const auto *gf_s_124 = buffer.data(gf_s + 124);
    const auto *gf_s_125 = buffer.data(gf_s + 125);
    const auto *gf_s_126 = buffer.data(gf_s + 126);
    const auto *gf_s_127 = buffer.data(gf_s + 127);
    const auto *gf_s_128 = buffer.data(gf_s + 128);
    const auto *gf_s_129 = buffer.data(gf_s + 129);
    const auto *gf_s_140 = buffer.data(gf_s + 140);
    const auto *gf_s_142 = buffer.data(gf_s + 142);
    const auto *gf_s_143 = buffer.data(gf_s + 143);
    const auto *gf_s_145 = buffer.data(gf_s + 145);
    const auto *gf_s_146 = buffer.data(gf_s + 146);
    const auto *gf_s_147 = buffer.data(gf_s + 147);
    const auto *gf_s_148 = buffer.data(gf_s + 148);
    const auto *gf_s_149 = buffer.data(gf_s + 149);

    const auto *gh_s_217 = buffer.data(gh_s + 217);
    const auto *gh_s_218 = buffer.data(gh_s + 218);
    const auto *gh_s_219 = buffer.data(gh_s + 219);
    const auto *gh_s_220 = buffer.data(gh_s + 220);
    const auto *gh_s_221 = buffer.data(gh_s + 221);
    const auto *gh_s_222 = buffer.data(gh_s + 222);
    const auto *gh_s_223 = buffer.data(gh_s + 223);
    const auto *gh_s_224 = buffer.data(gh_s + 224);
    const auto *gh_s_225 = buffer.data(gh_s + 225);
    const auto *gh_s_226 = buffer.data(gh_s + 226);
    const auto *gh_s_227 = buffer.data(gh_s + 227);
    const auto *gh_s_228 = buffer.data(gh_s + 228);
    const auto *gh_s_229 = buffer.data(gh_s + 229);
    const auto *gh_s_230 = buffer.data(gh_s + 230);
    const auto *gh_s_231 = buffer.data(gh_s + 231);
    const auto *gh_s_232 = buffer.data(gh_s + 232);
    const auto *gh_s_233 = buffer.data(gh_s + 233);
    const auto *gh_s_234 = buffer.data(gh_s + 234);
    const auto *gh_s_235 = buffer.data(gh_s + 235);
    const auto *gh_s_236 = buffer.data(gh_s + 236);
    const auto *gh_s_237 = buffer.data(gh_s + 237);
    const auto *gh_s_238 = buffer.data(gh_s + 238);
    const auto *gh_s_239 = buffer.data(gh_s + 239);
    const auto *gh_s_240 = buffer.data(gh_s + 240);
    const auto *gh_s_241 = buffer.data(gh_s + 241);
    const auto *gh_s_242 = buffer.data(gh_s + 242);
    const auto *gh_s_243 = buffer.data(gh_s + 243);
    const auto *gh_s_244 = buffer.data(gh_s + 244);
    const auto *gh_s_245 = buffer.data(gh_s + 245);
    const auto *gh_s_246 = buffer.data(gh_s + 246);
    const auto *gh_s_247 = buffer.data(gh_s + 247);
    const auto *gh_s_248 = buffer.data(gh_s + 248);
    const auto *gh_s_249 = buffer.data(gh_s + 249);
    const auto *gh_s_250 = buffer.data(gh_s + 250);
    const auto *gh_s_251 = buffer.data(gh_s + 251);
    const auto *gh_s_252 = buffer.data(gh_s + 252);
    const auto *gh_s_253 = buffer.data(gh_s + 253);
    const auto *gh_s_254 = buffer.data(gh_s + 254);
    const auto *gh_s_255 = buffer.data(gh_s + 255);
    const auto *gh_s_256 = buffer.data(gh_s + 256);
    const auto *gh_s_257 = buffer.data(gh_s + 257);
    const auto *gh_s_258 = buffer.data(gh_s + 258);
    const auto *gh_s_259 = buffer.data(gh_s + 259);
    const auto *gh_s_260 = buffer.data(gh_s + 260);
    const auto *gh_s_261 = buffer.data(gh_s + 261);
    const auto *gh_s_262 = buffer.data(gh_s + 262);
    const auto *gh_s_263 = buffer.data(gh_s + 263);
    const auto *gh_s_264 = buffer.data(gh_s + 264);
    const auto *gh_s_265 = buffer.data(gh_s + 265);
    const auto *gh_s_266 = buffer.data(gh_s + 266);
    const auto *gh_s_267 = buffer.data(gh_s + 267);
    const auto *gh_s_268 = buffer.data(gh_s + 268);
    const auto *gh_s_269 = buffer.data(gh_s + 269);
    const auto *gh_s_270 = buffer.data(gh_s + 270);
    const auto *gh_s_271 = buffer.data(gh_s + 271);
    const auto *gh_s_272 = buffer.data(gh_s + 272);
    const auto *gh_s_273 = buffer.data(gh_s + 273);
    const auto *gh_s_274 = buffer.data(gh_s + 274);
    const auto *gh_s_275 = buffer.data(gh_s + 275);
    const auto *gh_s_276 = buffer.data(gh_s + 276);
    const auto *gh_s_277 = buffer.data(gh_s + 277);
    const auto *gh_s_278 = buffer.data(gh_s + 278);
    const auto *gh_s_279 = buffer.data(gh_s + 279);
    const auto *gh_s_280 = buffer.data(gh_s + 280);
    const auto *gh_s_281 = buffer.data(gh_s + 281);
    const auto *gh_s_282 = buffer.data(gh_s + 282);
    const auto *gh_s_283 = buffer.data(gh_s + 283);
    const auto *gh_s_284 = buffer.data(gh_s + 284);
    const auto *gh_s_285 = buffer.data(gh_s + 285);
    const auto *gh_s_286 = buffer.data(gh_s + 286);
    const auto *gh_s_287 = buffer.data(gh_s + 287);
    const auto *gh_s_288 = buffer.data(gh_s + 288);
    const auto *gh_s_289 = buffer.data(gh_s + 289);
    const auto *gh_s_290 = buffer.data(gh_s + 290);
    const auto *gh_s_291 = buffer.data(gh_s + 291);
    const auto *gh_s_292 = buffer.data(gh_s + 292);
    const auto *gh_s_293 = buffer.data(gh_s + 293);
    const auto *gh_s_294 = buffer.data(gh_s + 294);
    const auto *gh_s_295 = buffer.data(gh_s + 295);
    const auto *gh_s_296 = buffer.data(gh_s + 296);
    const auto *gh_s_297 = buffer.data(gh_s + 297);
    const auto *gh_s_298 = buffer.data(gh_s + 298);
    const auto *gh_s_299 = buffer.data(gh_s + 299);
    const auto *gh_s_300 = buffer.data(gh_s + 300);
    const auto *gh_s_301 = buffer.data(gh_s + 301);
    const auto *gh_s_302 = buffer.data(gh_s + 302);
    const auto *gh_s_303 = buffer.data(gh_s + 303);
    const auto *gh_s_304 = buffer.data(gh_s + 304);
    const auto *gh_s_305 = buffer.data(gh_s + 305);
    const auto *gh_s_306 = buffer.data(gh_s + 306);
    const auto *gh_s_307 = buffer.data(gh_s + 307);
    const auto *gh_s_308 = buffer.data(gh_s + 308);
    const auto *gh_s_309 = buffer.data(gh_s + 309);
    const auto *gh_s_310 = buffer.data(gh_s + 310);
    const auto *gh_s_311 = buffer.data(gh_s + 311);
    const auto *gh_s_312 = buffer.data(gh_s + 312);
    const auto *gh_s_313 = buffer.data(gh_s + 313);
    const auto *gh_s_314 = buffer.data(gh_s + 314);

    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_107 = buffer.data(gf + 107);
    const auto *gf_108 = buffer.data(gf + 108);
    const auto *gf_109 = buffer.data(gf + 109);
    const auto *gf_112 = buffer.data(gf + 112);
    const auto *gf_115 = buffer.data(gf + 115);
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
    const auto *gf_140 = buffer.data(gf + 140);
    const auto *gf_142 = buffer.data(gf + 142);
    const auto *gf_143 = buffer.data(gf + 143);
    const auto *gf_145 = buffer.data(gf + 145);
    const auto *gf_146 = buffer.data(gf + 146);
    const auto *gf_147 = buffer.data(gf + 147);
    const auto *gf_148 = buffer.data(gf + 148);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_158 = buffer.data(gg + 158);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_167 = buffer.data(gg + 167);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_181 = buffer.data(gg + 181);
    const auto *gg_182 = buffer.data(gg + 182);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_184 = buffer.data(gg + 184);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_187 = buffer.data(gg + 187);
    const auto *gg_188 = buffer.data(gg + 188);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

#pragma omp simd aligned(t_217, t_218, t_219, pb_x, pb_z, gf_s_108, gf_s_109, gh_s_217, \
                         gh_s_218, gh_s_219, gf_108, gf_109, gg_153, gg_158, \
                         gg_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_2 * gh_s_217[k]
                   + pb_z[k] * gg_153[k];

        t_218[k] = -f_3 * gf_s_108[k]
                   + f_2 * gh_s_218[k]
                   + f_4 * gf_108[k]
                   + pb_x[k] * gg_158[k];

        t_219[k] = -f_3 * gf_s_109[k]
                   + f_2 * gh_s_219[k]
                   + f_4 * gf_109[k]
                   + pb_x[k] * gg_159[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, t_224, pb_x, gh_s_220, gh_s_221, \
                         gh_s_222, gh_s_223, gh_s_224, gg_160, gg_161, gg_162, gg_163, \
                         gg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = f_2 * gh_s_220[k]
                   + pb_x[k] * gg_160[k];

        t_221[k] = f_2 * gh_s_221[k]
                   + pb_x[k] * gg_161[k];

        t_222[k] = f_2 * gh_s_222[k]
                   + pb_x[k] * gg_162[k];

        t_223[k] = f_2 * gh_s_223[k]
                   + pb_x[k] * gg_163[k];

        t_224[k] = f_2 * gh_s_224[k]
                   + pb_x[k] * gg_164[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pb_y, pb_z, fg_100, gf_s_106, gh_s_225, \
                         gh_s_226, gh_s_227, gf_106, gg_160, gg_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_0 * fg_100[k]
                   - f_1 * gf_s_106[k]
                   + f_2 * gh_s_225[k]
                   + f_0 * gf_106[k]
                   + pb_y[k] * gg_160[k];

        t_226[k] = f_2 * gh_s_226[k]
                   + pb_z[k] * gg_160[k];

        t_227[k] = -f_3 * gf_s_106[k]
                   + f_2 * gh_s_227[k]
                   + f_4 * gf_106[k]
                   + pb_z[k] * gg_161[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pb_y, pb_z, fg_104, gf_s_107, gf_s_109, \
                         gh_s_228, gh_s_229, gh_s_230, gf_107, gf_109, gg_162, \
                         gg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = -f_5 * gf_s_107[k]
                   + f_2 * gh_s_228[k]
                   + f_6 * gf_107[k]
                   + pb_z[k] * gg_162[k];

        t_229[k] = f_0 * fg_104[k]
                   + f_2 * gh_s_229[k]
                   + pb_y[k] * gg_164[k];

        t_230[k] = -f_1 * gf_s_109[k]
                   + f_2 * gh_s_230[k]
                   + f_0 * gf_109[k]
                   + pb_z[k] * gg_164[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, t_234, pa_z, pb_x, fh_126, fh_127, fh_129, \
                         gf_s_112, gh_s_231, gh_s_232, gh_s_233, gh_s_234, gf_112, \
                         gg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = pa_z[k] * fh_126[k]
                   + f_2 * gh_s_231[k];

        t_232[k] = pa_z[k] * fh_127[k]
                   + f_2 * gh_s_232[k];

        t_233[k] = -f_9 * gf_s_112[k]
                   + f_2 * gh_s_233[k]
                   + f_7 * gf_112[k]
                   + pb_x[k] * gg_167[k];

        t_234[k] = pa_z[k] * fh_129[k]
                   + f_2 * gh_s_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_z, pb_x, fg_91, fh_130, fh_132, gf_s_115, \
                         gh_s_235, gh_s_236, gh_s_237, gf_115, gg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_4 * fg_91[k]
                   + pa_z[k] * fh_130[k]
                   + f_2 * gh_s_235[k];

        t_236[k] = -f_5 * gf_s_115[k]
                   + f_2 * gh_s_236[k]
                   + f_6 * gf_115[k]
                   + pb_x[k] * gg_170[k];

        t_237[k] = pa_z[k] * fh_132[k]
                   + f_2 * gh_s_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_z, pb_x, fg_93, fg_94, fh_133, fh_134, \
                         gf_s_119, gh_s_238, gh_s_239, gh_s_240, gf_119, \
                         gg_174 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_4 * fg_93[k]
                   + pa_z[k] * fh_133[k]
                   + f_2 * gh_s_238[k];

        t_239[k] = f_6 * fg_94[k]
                   + pa_z[k] * fh_134[k]
                   + f_2 * gh_s_239[k];

        t_240[k] = -f_3 * gf_s_119[k]
                   + f_2 * gh_s_240[k]
                   + f_4 * gf_119[k]
                   + pb_x[k] * gg_174[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, t_244, t_245, pb_x, gh_s_241, gh_s_242, \
                         gh_s_243, gh_s_244, gh_s_245, gg_175, gg_176, gg_177, gg_178, \
                         gg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = f_2 * gh_s_241[k]
                   + pb_x[k] * gg_175[k];

        t_242[k] = f_2 * gh_s_242[k]
                   + pb_x[k] * gg_176[k];

        t_243[k] = f_2 * gh_s_243[k]
                   + pb_x[k] * gg_177[k];

        t_244[k] = f_2 * gh_s_244[k]
                   + pb_x[k] * gg_178[k];

        t_245[k] = f_2 * gh_s_245[k]
                   + pb_x[k] * gg_179[k];
    }

#pragma omp simd aligned(t_246, t_247, t_248, pa_z, pb_z, fg_100, fg_101, fh_141, fh_143, \
                         gh_s_246, gh_s_247, gh_s_248, gg_175 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_246[k] = pa_z[k] * fh_141[k]
                   + f_2 * gh_s_246[k];

        t_247[k] = f_4 * fg_100[k]
                   + f_2 * gh_s_247[k]
                   + pb_z[k] * gg_175[k];

        t_248[k] = f_6 * fg_101[k]
                   + pa_z[k] * fh_143[k]
                   + f_2 * gh_s_248[k];
    }

#pragma omp simd aligned(t_249, t_250, t_251, pa_y, pa_z, pb_y, dh_s_104, dh_104, fg_102, \
                         fg_119, fh_144, fh_167, gh_s_249, gh_s_250, gh_s_251, \
                         gg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_249[k] = f_7 * fg_102[k]
                   + pa_z[k] * fh_144[k]
                   + f_2 * gh_s_249[k];

        t_250[k] = f_7 * fg_119[k]
                   + f_2 * gh_s_250[k]
                   + pb_y[k] * gg_179[k];

        t_251[k] = -f_8 * dh_s_104[k]
                   + f_6 * dh_104[k]
                   + pa_y[k] * fh_167[k]
                   + f_2 * gh_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pb_x, gf_s_120, gf_s_121, gf_s_122, gh_s_252, \
                         gh_s_253, gh_s_254, gf_120, gf_121, gf_122, gg_180, gg_181, \
                         gg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = -f_1 * gf_s_120[k]
                   + f_2 * gh_s_252[k]
                   + f_0 * gf_120[k]
                   + pb_x[k] * gg_180[k];

        t_253[k] = -f_9 * gf_s_121[k]
                   + f_2 * gh_s_253[k]
                   + f_7 * gf_121[k]
                   + pb_x[k] * gg_181[k];

        t_254[k] = -f_9 * gf_s_122[k]
                   + f_2 * gh_s_254[k]
                   + f_7 * gf_122[k]
                   + pb_x[k] * gg_182[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pb_x, gf_s_123, gf_s_124, gf_s_125, gh_s_255, \
                         gh_s_256, gh_s_257, gf_123, gf_124, gf_125, gg_183, gg_184, \
                         gg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = -f_5 * gf_s_123[k]
                   + f_2 * gh_s_255[k]
                   + f_6 * gf_123[k]
                   + pb_x[k] * gg_183[k];

        t_256[k] = -f_5 * gf_s_124[k]
                   + f_2 * gh_s_256[k]
                   + f_6 * gf_124[k]
                   + pb_x[k] * gg_184[k];

        t_257[k] = -f_5 * gf_s_125[k]
                   + f_2 * gh_s_257[k]
                   + f_6 * gf_125[k]
                   + pb_x[k] * gg_185[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pb_x, gf_s_126, gf_s_127, gf_s_128, gh_s_258, \
                         gh_s_259, gh_s_260, gf_126, gf_127, gf_128, gg_186, gg_187, \
                         gg_188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = -f_3 * gf_s_126[k]
                   + f_2 * gh_s_258[k]
                   + f_4 * gf_126[k]
                   + pb_x[k] * gg_186[k];

        t_259[k] = -f_3 * gf_s_127[k]
                   + f_2 * gh_s_259[k]
                   + f_4 * gf_127[k]
                   + pb_x[k] * gg_187[k];

        t_260[k] = -f_3 * gf_s_128[k]
                   + f_2 * gh_s_260[k]
                   + f_4 * gf_128[k]
                   + pb_x[k] * gg_188[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, t_264, pb_x, gf_s_129, gh_s_261, gh_s_262, \
                         gh_s_263, gh_s_264, gf_129, gg_189, gg_190, gg_191, \
                         gg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = -f_3 * gf_s_129[k]
                   + f_2 * gh_s_261[k]
                   + f_4 * gf_129[k]
                   + pb_x[k] * gg_189[k];

        t_262[k] = f_2 * gh_s_262[k]
                   + pb_x[k] * gg_190[k];

        t_263[k] = f_2 * gh_s_263[k]
                   + pb_x[k] * gg_191[k];

        t_264[k] = f_2 * gh_s_264[k]
                   + pb_x[k] * gg_192[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, pa_z, pb_x, dh_s_78, dh_78, fh_162, gh_s_265, \
                         gh_s_266, gh_s_267, gg_193, gg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = f_2 * gh_s_265[k]
                   + pb_x[k] * gg_193[k];

        t_266[k] = f_2 * gh_s_266[k]
                   + pb_x[k] * gg_194[k];

        t_267[k] = -f_10 * dh_s_78[k]
                   + f_4 * dh_78[k]
                   + pa_z[k] * fh_162[k]
                   + f_2 * gh_s_267[k];
    }

#pragma omp simd aligned(t_268, t_269, pb_y, pb_z, fg_115, fg_132, gf_s_128, gh_s_268, \
                         gh_s_269, gf_128, gg_190, gg_192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_268[k] = f_6 * fg_115[k]
                   + f_2 * gh_s_268[k]
                   + pb_z[k] * gg_190[k];

        t_269[k] = f_6 * fg_132[k]
                   - f_5 * gf_s_128[k]
                   + f_2 * gh_s_269[k]
                   + f_6 * gf_128[k]
                   + pb_y[k] * gg_192[k];
    }

#pragma omp simd aligned(t_270, t_271, pb_y, fg_133, fg_134, gf_s_129, gh_s_270, gh_s_271, \
                         gf_129, gg_193, gg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_6 * fg_133[k]
                   - f_3 * gf_s_129[k]
                   + f_2 * gh_s_270[k]
                   + f_4 * gf_129[k]
                   + pb_y[k] * gg_193[k];

        t_271[k] = f_6 * fg_134[k]
                   + f_2 * gh_s_271[k]
                   + pb_y[k] * gg_194[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_y, dh_s_125, dh_125, fg_135, fh_188, \
                         fh_189, fh_190, fh_191, gh_s_272, gh_s_273, gh_s_274, \
                         gh_s_275 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = -f_10 * dh_s_125[k]
                   + f_4 * dh_125[k]
                   + pa_y[k] * fh_188[k]
                   + f_2 * gh_s_272[k];

        t_273[k] = pa_y[k] * fh_189[k]
                   + f_2 * gh_s_273[k];

        t_274[k] = f_4 * fg_135[k]
                   + pa_y[k] * fh_190[k]
                   + f_2 * gh_s_274[k];

        t_275[k] = pa_y[k] * fh_191[k]
                   + f_2 * gh_s_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, t_279, pa_y, fg_136, fg_137, fg_138, fh_192, \
                         fh_193, fh_194, fh_195, gh_s_276, gh_s_277, gh_s_278, \
                         gh_s_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_6 * fg_136[k]
                   + pa_y[k] * fh_192[k]
                   + f_2 * gh_s_276[k];

        t_277[k] = f_4 * fg_137[k]
                   + pa_y[k] * fh_193[k]
                   + f_2 * gh_s_277[k];

        t_278[k] = pa_y[k] * fh_194[k]
                   + f_2 * gh_s_278[k];

        t_279[k] = f_7 * fg_138[k]
                   + pa_y[k] * fh_195[k]
                   + f_2 * gh_s_279[k];
    }

#pragma omp simd aligned(t_280, t_281, t_282, t_283, pa_y, pb_x, fg_139, fg_140, fh_196, \
                         fh_197, fh_198, gh_s_280, gh_s_281, gh_s_282, gh_s_283, \
                         gg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_280[k] = f_6 * fg_139[k]
                   + pa_y[k] * fh_196[k]
                   + f_2 * gh_s_280[k];

        t_281[k] = f_4 * fg_140[k]
                   + pa_y[k] * fh_197[k]
                   + f_2 * gh_s_281[k];

        t_282[k] = pa_y[k] * fh_198[k]
                   + f_2 * gh_s_282[k];

        t_283[k] = f_2 * gh_s_283[k]
                   + pb_x[k] * gg_205[k];
    }

#pragma omp simd aligned(t_284, t_285, t_286, t_287, pb_x, gh_s_284, gh_s_285, gh_s_286, \
                         gh_s_287, gg_206, gg_207, gg_208, gg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_284[k] = f_2 * gh_s_284[k]
                   + pb_x[k] * gg_206[k];

        t_285[k] = f_2 * gh_s_285[k]
                   + pb_x[k] * gg_207[k];

        t_286[k] = f_2 * gh_s_286[k]
                   + pb_x[k] * gg_208[k];

        t_287[k] = f_2 * gh_s_287[k]
                   + pb_x[k] * gg_209[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, pa_y, pb_z, fg_130, fg_145, fg_147, fh_204, \
                         fh_206, gh_s_288, gh_s_289, gh_s_290, gg_205 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = f_11 * fg_145[k]
                   + pa_y[k] * fh_204[k]
                   + f_2 * gh_s_288[k];

        t_289[k] = f_7 * fg_130[k]
                   + f_2 * gh_s_289[k]
                   + pb_z[k] * gg_205[k];

        t_290[k] = f_7 * fg_147[k]
                   + pa_y[k] * fh_206[k]
                   + f_2 * gh_s_290[k];
    }

#pragma omp simd aligned(t_291, t_292, t_293, pa_y, pb_y, fg_148, fg_149, fh_207, fh_209, \
                         gh_s_291, gh_s_292, gh_s_293, gg_209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_291[k] = f_6 * fg_148[k]
                   + pa_y[k] * fh_207[k]
                   + f_2 * gh_s_291[k];

        t_292[k] = f_4 * fg_149[k]
                   + f_2 * gh_s_292[k]
                   + pb_y[k] * gg_209[k];

        t_293[k] = pa_y[k] * fh_209[k]
                   + f_2 * gh_s_293[k];
    }

#pragma omp simd aligned(t_294, t_295, t_296, pb_x, pb_y, gf_s_140, gf_s_142, gh_s_294, \
                         gh_s_295, gh_s_296, gf_140, gf_142, gg_210, \
                         gg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_294[k] = -f_1 * gf_s_140[k]
                   + f_2 * gh_s_294[k]
                   + f_0 * gf_140[k]
                   + pb_x[k] * gg_210[k];

        t_295[k] = f_2 * gh_s_295[k]
                   + pb_y[k] * gg_210[k];

        t_296[k] = -f_9 * gf_s_142[k]
                   + f_2 * gh_s_296[k]
                   + f_7 * gf_142[k]
                   + pb_x[k] * gg_212[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pb_x, pb_y, gf_s_143, gf_s_145, gh_s_297, \
                         gh_s_298, gh_s_299, gf_143, gf_145, gg_212, gg_213, \
                         gg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = -f_5 * gf_s_143[k]
                   + f_2 * gh_s_297[k]
                   + f_6 * gf_143[k]
                   + pb_x[k] * gg_213[k];

        t_298[k] = f_2 * gh_s_298[k]
                   + pb_y[k] * gg_212[k];

        t_299[k] = -f_5 * gf_s_145[k]
                   + f_2 * gh_s_299[k]
                   + f_6 * gf_145[k]
                   + pb_x[k] * gg_215[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pb_x, pb_y, gf_s_146, gf_s_147, gh_s_300, \
                         gh_s_301, gh_s_302, gf_146, gf_147, gg_215, gg_216, \
                         gg_217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = -f_3 * gf_s_146[k]
                   + f_2 * gh_s_300[k]
                   + f_4 * gf_146[k]
                   + pb_x[k] * gg_216[k];

        t_301[k] = -f_3 * gf_s_147[k]
                   + f_2 * gh_s_301[k]
                   + f_4 * gf_147[k]
                   + pb_x[k] * gg_217[k];

        t_302[k] = f_2 * gh_s_302[k]
                   + pb_y[k] * gg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, t_306, pb_x, gf_s_149, gh_s_303, gh_s_304, \
                         gh_s_305, gh_s_306, gf_149, gg_219, gg_220, gg_221, \
                         gg_222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = -f_3 * gf_s_149[k]
                   + f_2 * gh_s_303[k]
                   + f_4 * gf_149[k]
                   + pb_x[k] * gg_219[k];

        t_304[k] = f_2 * gh_s_304[k]
                   + pb_x[k] * gg_220[k];

        t_305[k] = f_2 * gh_s_305[k]
                   + pb_x[k] * gg_221[k];

        t_306[k] = f_2 * gh_s_306[k]
                   + pb_x[k] * gg_222[k];
    }

#pragma omp simd aligned(t_307, t_308, t_309, pb_x, pb_y, gf_s_146, gh_s_307, gh_s_308, \
                         gh_s_309, gf_146, gg_220, gg_223, gg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_307[k] = f_2 * gh_s_307[k]
                   + pb_x[k] * gg_223[k];

        t_308[k] = f_2 * gh_s_308[k]
                   + pb_x[k] * gg_224[k];

        t_309[k] = -f_1 * gf_s_146[k]
                   + f_2 * gh_s_309[k]
                   + f_0 * gf_146[k]
                   + pb_y[k] * gg_220[k];
    }

#pragma omp simd aligned(t_310, t_311, t_312, pb_y, gf_s_147, gf_s_148, gf_s_149, gh_s_310, \
                         gh_s_311, gh_s_312, gf_147, gf_148, gf_149, gg_221, gg_222, \
                         gg_223 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_310[k] = -f_9 * gf_s_147[k]
                   + f_2 * gh_s_310[k]
                   + f_7 * gf_147[k]
                   + pb_y[k] * gg_221[k];

        t_311[k] = -f_5 * gf_s_148[k]
                   + f_2 * gh_s_311[k]
                   + f_6 * gf_148[k]
                   + pb_y[k] * gg_222[k];

        t_312[k] = -f_3 * gf_s_149[k]
                   + f_2 * gh_s_312[k]
                   + f_4 * gf_149[k]
                   + pb_y[k] * gg_223[k];
    }

#pragma omp simd aligned(t_313, t_314, pb_y, pb_z, fg_149, gf_s_149, gh_s_313, gh_s_314, \
                         gf_149, gg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_2 * gh_s_313[k]
                   + pb_y[k] * gg_224[k];

        t_314[k] = f_0 * fg_149[k]
                   - f_1 * gf_s_149[k]
                   + f_2 * gh_s_314[k]
                   + f_0 * gf_149[k]
                   + pb_z[k] * gg_224[k];
    }
}

auto
compute_prim_gh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t dh_s, const size_t dh,
                                 const size_t fg, const size_t fh, const size_t gf_s,
                                 const size_t gh_s, const size_t gf, const size_t gg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_gh_kinetic_energy_0_piece0(buffer, target, pa, pb, dh_s, dh, fg, fh, gf_s,
                                            gh_s, gf, gg, ncols, alpha, beta, p);

    compute_prim_gh_kinetic_energy_0_piece1(buffer, target, pa, pb, dh_s, dh, fg, fh, gf_s,
                                            gh_s, gf, gg, ncols, alpha, beta, p);

    compute_prim_gh_kinetic_energy_0_piece2(buffer, target, pa, pb, dh_s, dh, fg, fh, gf_s,
                                            gh_s, gf, gg, ncols, alpha, beta, p);
}

}  // namespace simdkin
