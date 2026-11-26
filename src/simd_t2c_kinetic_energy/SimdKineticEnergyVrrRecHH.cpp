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


#include "SimdKineticEnergyVrrRecHH.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_hh_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fh_s,
                                        const size_t fh, const size_t gg, const size_t gh,
                                        const size_t hf_s, const size_t hh_s, const size_t hf,
                                        const size_t hg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 3.0 * beta / p;
    const auto f_10 = 3.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 2.0 * beta / p;

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

    const auto *fh_s_0 = buffer.data(fh_s + 0);
    const auto *fh_s_36 = buffer.data(fh_s + 36);
    const auto *fh_s_62 = buffer.data(fh_s + 62);
    const auto *fh_s_78 = buffer.data(fh_s + 78);
    const auto *fh_s_101 = buffer.data(fh_s + 101);
    const auto *fh_s_102 = buffer.data(fh_s + 102);

    const auto *fh_0 = buffer.data(fh + 0);
    const auto *fh_36 = buffer.data(fh + 36);
    const auto *fh_62 = buffer.data(fh + 62);
    const auto *fh_78 = buffer.data(fh + 78);
    const auto *fh_101 = buffer.data(fh + 101);
    const auto *fh_102 = buffer.data(fh + 102);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_1 = buffer.data(gg + 1);
    const auto *gg_2 = buffer.data(gg + 2);
    const auto *gg_3 = buffer.data(gg + 3);
    const auto *gg_5 = buffer.data(gg + 5);
    const auto *gg_10 = buffer.data(gg + 10);
    const auto *gg_12 = buffer.data(gg + 12);
    const auto *gg_14 = buffer.data(gg + 14);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_18 = buffer.data(gg + 18);
    const auto *gg_20 = buffer.data(gg + 20);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_27 = buffer.data(gg + 27);
    const auto *gg_28 = buffer.data(gg + 28);
    const auto *gg_29 = buffer.data(gg + 29);
    const auto *gg_32 = buffer.data(gg + 32);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_41 = buffer.data(gg + 41);
    const auto *gg_42 = buffer.data(gg + 42);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_51 = buffer.data(gg + 51);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_57 = buffer.data(gg + 57);
    const auto *gg_58 = buffer.data(gg + 58);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_71 = buffer.data(gg + 71);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_73 = buffer.data(gg + 73);

    const auto *gh_0 = buffer.data(gh + 0);
    const auto *gh_3 = buffer.data(gh + 3);
    const auto *gh_5 = buffer.data(gh + 5);
    const auto *gh_6 = buffer.data(gh + 6);
    const auto *gh_7 = buffer.data(gh + 7);
    const auto *gh_9 = buffer.data(gh + 9);
    const auto *gh_10 = buffer.data(gh + 10);
    const auto *gh_14 = buffer.data(gh + 14);
    const auto *gh_15 = buffer.data(gh + 15);
    const auto *gh_20 = buffer.data(gh + 20);
    const auto *gh_21 = buffer.data(gh + 21);
    const auto *gh_22 = buffer.data(gh + 22);
    const auto *gh_24 = buffer.data(gh + 24);
    const auto *gh_27 = buffer.data(gh + 27);
    const auto *gh_31 = buffer.data(gh + 31);
    const auto *gh_36 = buffer.data(gh + 36);
    const auto *gh_42 = buffer.data(gh + 42);
    const auto *gh_44 = buffer.data(gh + 44);
    const auto *gh_47 = buffer.data(gh + 47);
    const auto *gh_51 = buffer.data(gh + 51);
    const auto *gh_56 = buffer.data(gh + 56);
    const auto *gh_62 = buffer.data(gh + 62);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_101 = buffer.data(gh + 101);
    const auto *gh_102 = buffer.data(gh + 102);

    const auto *hf_s_0 = buffer.data(hf_s + 0);
    const auto *hf_s_1 = buffer.data(hf_s + 1);
    const auto *hf_s_2 = buffer.data(hf_s + 2);
    const auto *hf_s_6 = buffer.data(hf_s + 6);
    const auto *hf_s_8 = buffer.data(hf_s + 8);
    const auto *hf_s_9 = buffer.data(hf_s + 9);
    const auto *hf_s_16 = buffer.data(hf_s + 16);
    const auto *hf_s_17 = buffer.data(hf_s + 17);
    const auto *hf_s_27 = buffer.data(hf_s + 27);
    const auto *hf_s_28 = buffer.data(hf_s + 28);
    const auto *hf_s_29 = buffer.data(hf_s + 29);
    const auto *hf_s_30 = buffer.data(hf_s + 30);
    const auto *hf_s_32 = buffer.data(hf_s + 32);
    const auto *hf_s_33 = buffer.data(hf_s + 33);
    const auto *hf_s_36 = buffer.data(hf_s + 36);
    const auto *hf_s_37 = buffer.data(hf_s + 37);
    const auto *hf_s_39 = buffer.data(hf_s + 39);

    const auto *hh_s_0 = buffer.data(hh_s + 0);
    const auto *hh_s_1 = buffer.data(hh_s + 1);
    const auto *hh_s_2 = buffer.data(hh_s + 2);
    const auto *hh_s_3 = buffer.data(hh_s + 3);
    const auto *hh_s_4 = buffer.data(hh_s + 4);
    const auto *hh_s_5 = buffer.data(hh_s + 5);
    const auto *hh_s_6 = buffer.data(hh_s + 6);
    const auto *hh_s_7 = buffer.data(hh_s + 7);
    const auto *hh_s_8 = buffer.data(hh_s + 8);
    const auto *hh_s_9 = buffer.data(hh_s + 9);
    const auto *hh_s_10 = buffer.data(hh_s + 10);
    const auto *hh_s_11 = buffer.data(hh_s + 11);
    const auto *hh_s_12 = buffer.data(hh_s + 12);
    const auto *hh_s_13 = buffer.data(hh_s + 13);
    const auto *hh_s_14 = buffer.data(hh_s + 14);
    const auto *hh_s_15 = buffer.data(hh_s + 15);
    const auto *hh_s_16 = buffer.data(hh_s + 16);
    const auto *hh_s_17 = buffer.data(hh_s + 17);
    const auto *hh_s_18 = buffer.data(hh_s + 18);
    const auto *hh_s_19 = buffer.data(hh_s + 19);
    const auto *hh_s_20 = buffer.data(hh_s + 20);
    const auto *hh_s_21 = buffer.data(hh_s + 21);
    const auto *hh_s_22 = buffer.data(hh_s + 22);
    const auto *hh_s_23 = buffer.data(hh_s + 23);
    const auto *hh_s_24 = buffer.data(hh_s + 24);
    const auto *hh_s_25 = buffer.data(hh_s + 25);
    const auto *hh_s_26 = buffer.data(hh_s + 26);
    const auto *hh_s_27 = buffer.data(hh_s + 27);
    const auto *hh_s_28 = buffer.data(hh_s + 28);
    const auto *hh_s_29 = buffer.data(hh_s + 29);
    const auto *hh_s_30 = buffer.data(hh_s + 30);
    const auto *hh_s_31 = buffer.data(hh_s + 31);
    const auto *hh_s_32 = buffer.data(hh_s + 32);
    const auto *hh_s_33 = buffer.data(hh_s + 33);
    const auto *hh_s_34 = buffer.data(hh_s + 34);
    const auto *hh_s_35 = buffer.data(hh_s + 35);
    const auto *hh_s_36 = buffer.data(hh_s + 36);
    const auto *hh_s_37 = buffer.data(hh_s + 37);
    const auto *hh_s_38 = buffer.data(hh_s + 38);
    const auto *hh_s_39 = buffer.data(hh_s + 39);
    const auto *hh_s_40 = buffer.data(hh_s + 40);
    const auto *hh_s_41 = buffer.data(hh_s + 41);
    const auto *hh_s_42 = buffer.data(hh_s + 42);
    const auto *hh_s_43 = buffer.data(hh_s + 43);
    const auto *hh_s_44 = buffer.data(hh_s + 44);
    const auto *hh_s_45 = buffer.data(hh_s + 45);
    const auto *hh_s_46 = buffer.data(hh_s + 46);
    const auto *hh_s_47 = buffer.data(hh_s + 47);
    const auto *hh_s_48 = buffer.data(hh_s + 48);
    const auto *hh_s_49 = buffer.data(hh_s + 49);
    const auto *hh_s_50 = buffer.data(hh_s + 50);
    const auto *hh_s_51 = buffer.data(hh_s + 51);
    const auto *hh_s_52 = buffer.data(hh_s + 52);
    const auto *hh_s_53 = buffer.data(hh_s + 53);
    const auto *hh_s_54 = buffer.data(hh_s + 54);
    const auto *hh_s_55 = buffer.data(hh_s + 55);
    const auto *hh_s_56 = buffer.data(hh_s + 56);
    const auto *hh_s_57 = buffer.data(hh_s + 57);
    const auto *hh_s_58 = buffer.data(hh_s + 58);
    const auto *hh_s_59 = buffer.data(hh_s + 59);
    const auto *hh_s_60 = buffer.data(hh_s + 60);
    const auto *hh_s_61 = buffer.data(hh_s + 61);
    const auto *hh_s_62 = buffer.data(hh_s + 62);
    const auto *hh_s_63 = buffer.data(hh_s + 63);
    const auto *hh_s_64 = buffer.data(hh_s + 64);
    const auto *hh_s_65 = buffer.data(hh_s + 65);
    const auto *hh_s_66 = buffer.data(hh_s + 66);
    const auto *hh_s_67 = buffer.data(hh_s + 67);
    const auto *hh_s_68 = buffer.data(hh_s + 68);
    const auto *hh_s_69 = buffer.data(hh_s + 69);
    const auto *hh_s_70 = buffer.data(hh_s + 70);
    const auto *hh_s_71 = buffer.data(hh_s + 71);
    const auto *hh_s_72 = buffer.data(hh_s + 72);
    const auto *hh_s_73 = buffer.data(hh_s + 73);
    const auto *hh_s_74 = buffer.data(hh_s + 74);
    const auto *hh_s_75 = buffer.data(hh_s + 75);
    const auto *hh_s_76 = buffer.data(hh_s + 76);
    const auto *hh_s_77 = buffer.data(hh_s + 77);
    const auto *hh_s_78 = buffer.data(hh_s + 78);
    const auto *hh_s_79 = buffer.data(hh_s + 79);
    const auto *hh_s_80 = buffer.data(hh_s + 80);
    const auto *hh_s_81 = buffer.data(hh_s + 81);
    const auto *hh_s_82 = buffer.data(hh_s + 82);
    const auto *hh_s_83 = buffer.data(hh_s + 83);
    const auto *hh_s_84 = buffer.data(hh_s + 84);
    const auto *hh_s_85 = buffer.data(hh_s + 85);
    const auto *hh_s_86 = buffer.data(hh_s + 86);
    const auto *hh_s_87 = buffer.data(hh_s + 87);
    const auto *hh_s_88 = buffer.data(hh_s + 88);
    const auto *hh_s_89 = buffer.data(hh_s + 89);
    const auto *hh_s_90 = buffer.data(hh_s + 90);
    const auto *hh_s_91 = buffer.data(hh_s + 91);
    const auto *hh_s_92 = buffer.data(hh_s + 92);
    const auto *hh_s_93 = buffer.data(hh_s + 93);
    const auto *hh_s_94 = buffer.data(hh_s + 94);
    const auto *hh_s_95 = buffer.data(hh_s + 95);
    const auto *hh_s_96 = buffer.data(hh_s + 96);
    const auto *hh_s_97 = buffer.data(hh_s + 97);
    const auto *hh_s_98 = buffer.data(hh_s + 98);
    const auto *hh_s_99 = buffer.data(hh_s + 99);
    const auto *hh_s_100 = buffer.data(hh_s + 100);
    const auto *hh_s_101 = buffer.data(hh_s + 101);
    const auto *hh_s_102 = buffer.data(hh_s + 102);
    const auto *hh_s_103 = buffer.data(hh_s + 103);
    const auto *hh_s_104 = buffer.data(hh_s + 104);
    const auto *hh_s_105 = buffer.data(hh_s + 105);
    const auto *hh_s_106 = buffer.data(hh_s + 106);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_8 = buffer.data(hf + 8);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_17 = buffer.data(hf + 17);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_28 = buffer.data(hf + 28);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_32 = buffer.data(hf + 32);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_37 = buffer.data(hf + 37);
    const auto *hf_39 = buffer.data(hf + 39);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_1 = buffer.data(hg + 1);
    const auto *hg_2 = buffer.data(hg + 2);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_12 = buffer.data(hg + 12);
    const auto *hg_13 = buffer.data(hg + 13);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_20 = buffer.data(hg + 20);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_26 = buffer.data(hg + 26);
    const auto *hg_27 = buffer.data(hg + 27);
    const auto *hg_28 = buffer.data(hg + 28);
    const auto *hg_29 = buffer.data(hg + 29);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_41 = buffer.data(hg + 41);
    const auto *hg_42 = buffer.data(hg + 42);
    const auto *hg_43 = buffer.data(hg + 43);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_47 = buffer.data(hg + 47);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_50 = buffer.data(hg + 50);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_56 = buffer.data(hg + 56);
    const auto *hg_57 = buffer.data(hg + 57);
    const auto *hg_58 = buffer.data(hg + 58);
    const auto *hg_59 = buffer.data(hg + 59);
    const auto *hg_62 = buffer.data(hg + 62);
    const auto *hg_63 = buffer.data(hg + 63);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_70 = buffer.data(hg + 70);
    const auto *hg_71 = buffer.data(hg + 71);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_73 = buffer.data(hg + 73);
    const auto *hg_74 = buffer.data(hg + 74);
    const auto *hg_75 = buffer.data(hg + 75);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, gg_0, hf_s_0, hh_s_0, hh_s_1, \
                         hh_s_2, hh_s_3, hf_0, hg_0, hg_1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * gg_0[k]
                 - f_1 * hf_s_0[k]
                 + f_2 * hh_s_0[k]
                 + f_3 * hf_0[k]
                 + pb_x[k] * hg_0[k];

        t_1[k] = f_2 * hh_s_1[k]
                 + pb_y[k] * hg_0[k];

        t_2[k] = f_2 * hh_s_2[k]
                 + pb_z[k] * hg_0[k];

        t_3[k] = -f_4 * hf_s_0[k]
                 + f_2 * hh_s_3[k]
                 + f_5 * hf_0[k]
                 + pb_y[k] * hg_1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, t_7, pb_y, pb_z, hf_s_0, hf_s_1, hh_s_4, hh_s_5, \
                         hh_s_6, hh_s_7, hf_0, hf_1, hg_2, hg_3 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * hh_s_4[k]
                 + pb_y[k] * hg_2[k];

        t_5[k] = -f_4 * hf_s_0[k]
                 + f_2 * hh_s_5[k]
                 + f_5 * hf_0[k]
                 + pb_z[k] * hg_2[k];

        t_6[k] = -f_6 * hf_s_1[k]
                 + f_2 * hh_s_6[k]
                 + f_7 * hf_1[k]
                 + pb_y[k] * hg_3[k];

        t_7[k] = f_2 * hh_s_7[k]
                 + pb_z[k] * hg_3[k];
    }

#pragma omp simd aligned(t_8, t_9, t_10, pb_x, pb_y, pb_z, gg_10, hf_s_2, hh_s_8, hh_s_9, \
                         hh_s_10, hf_2, hg_5, hg_10 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_8[k] = f_2 * hh_s_8[k]
                 + pb_y[k] * hg_5[k];

        t_9[k] = -f_6 * hf_s_2[k]
                 + f_2 * hh_s_9[k]
                 + f_7 * hf_2[k]
                 + pb_z[k] * hg_5[k];

        t_10[k] = f_0 * gg_10[k]
                  + f_2 * hh_s_10[k]
                  + pb_x[k] * hg_10[k];
    }

#pragma omp simd aligned(t_11, t_12, t_13, pb_x, pb_y, pb_z, gg_12, hh_s_11, hh_s_12, hh_s_13, \
                         hg_6, hg_9, hg_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_11[k] = f_2 * hh_s_11[k]
                  + pb_z[k] * hg_6[k];

        t_12[k] = f_0 * gg_12[k]
                  + f_2 * hh_s_12[k]
                  + pb_x[k] * hg_12[k];

        t_13[k] = f_2 * hh_s_13[k]
                  + pb_y[k] * hg_9[k];
    }

#pragma omp simd aligned(t_14, t_15, t_16, pb_x, pb_y, pb_z, gg_14, hf_s_6, hh_s_14, hh_s_15, \
                         hh_s_16, hf_6, hg_10, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_14[k] = f_0 * gg_14[k]
                  + f_2 * hh_s_14[k]
                  + pb_x[k] * hg_14[k];

        t_15[k] = -f_1 * hf_s_6[k]
                  + f_2 * hh_s_15[k]
                  + f_3 * hf_6[k]
                  + pb_y[k] * hg_10[k];

        t_16[k] = f_2 * hh_s_16[k]
                  + pb_z[k] * hg_10[k];
    }

#pragma omp simd aligned(t_17, t_18, t_19, pb_y, hf_s_8, hf_s_9, hh_s_17, hh_s_18, hh_s_19, \
                         hf_8, hf_9, hg_12, hg_13, hg_14 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_17[k] = -f_6 * hf_s_8[k]
                  + f_2 * hh_s_17[k]
                  + f_7 * hf_8[k]
                  + pb_y[k] * hg_12[k];

        t_18[k] = -f_4 * hf_s_9[k]
                  + f_2 * hh_s_18[k]
                  + f_5 * hf_9[k]
                  + pb_y[k] * hg_13[k];

        t_19[k] = f_2 * hh_s_19[k]
                  + pb_y[k] * hg_14[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_y, pb_z, gg_0, gh_0, hf_s_9, hh_s_20, \
                         hh_s_21, hh_s_22, hf_9, hg_14, hg_15 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = -f_1 * hf_s_9[k]
                  + f_2 * hh_s_20[k]
                  + f_3 * hf_9[k]
                  + pb_z[k] * hg_14[k];

        t_21[k] = pa_y[k] * gh_0[k]
                  + f_2 * hh_s_21[k];

        t_22[k] = f_5 * gg_0[k]
                  + f_2 * hh_s_22[k]
                  + pb_y[k] * hg_15[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_y, pb_z, gg_1, gh_3, gh_5, hh_s_23, \
                         hh_s_24, hh_s_25, hh_s_26, hg_15, hg_16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_2 * hh_s_23[k]
                  + pb_z[k] * hg_15[k];

        t_24[k] = f_7 * gg_1[k]
                  + pa_y[k] * gh_3[k]
                  + f_2 * hh_s_24[k];

        t_25[k] = f_2 * hh_s_25[k]
                  + pb_z[k] * hg_16[k];

        t_26[k] = pa_y[k] * gh_5[k]
                  + f_2 * hh_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_y, pb_y, pb_z, gg_3, gg_5, gh_6, hh_s_27, \
                         hh_s_28, hh_s_29, hg_18, hg_20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = f_8 * gg_3[k]
                  + pa_y[k] * gh_6[k]
                  + f_2 * hh_s_27[k];

        t_28[k] = f_2 * hh_s_28[k]
                  + pb_z[k] * hg_18[k];

        t_29[k] = f_5 * gg_5[k]
                  + f_2 * hh_s_29[k]
                  + pb_y[k] * hg_20[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_x, pb_z, gg_25, gh_9, hh_s_30, hh_s_31, \
                         hh_s_32, hg_21, hg_25 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = pa_y[k] * gh_9[k]
                  + f_2 * hh_s_30[k];

        t_31[k] = f_3 * gg_25[k]
                  + f_2 * hh_s_31[k]
                  + pb_x[k] * hg_25[k];

        t_32[k] = f_2 * hh_s_32[k]
                  + pb_z[k] * hg_21[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pa_y, pb_x, gg_27, gg_28, gh_14, hh_s_33, hh_s_34, \
                         hh_s_35, hg_27, hg_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_3 * gg_27[k]
                  + f_2 * hh_s_33[k]
                  + pb_x[k] * hg_27[k];

        t_34[k] = f_3 * gg_28[k]
                  + f_2 * hh_s_34[k]
                  + pb_x[k] * hg_28[k];

        t_35[k] = pa_y[k] * gh_14[k]
                  + f_2 * hh_s_35[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_z, fh_s_36, fh_36, gh_36, hf_s_16, \
                         hh_s_36, hh_s_37, hh_s_38, hf_16, hg_25, \
                         hg_26 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_9 * fh_s_36[k]
                  + f_8 * fh_36[k]
                  + pa_x[k] * gh_36[k]
                  + f_2 * hh_s_36[k];

        t_37[k] = f_2 * hh_s_37[k]
                  + pb_z[k] * hg_25[k];

        t_38[k] = -f_4 * hf_s_16[k]
                  + f_2 * hh_s_38[k]
                  + f_5 * hf_16[k]
                  + pb_z[k] * hg_26[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pb_y, pb_z, gg_14, gh_20, hf_s_17, hh_s_39, \
                         hh_s_40, hh_s_41, hf_17, hg_27, hg_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_6 * hf_s_17[k]
                  + f_2 * hh_s_39[k]
                  + f_7 * hf_17[k]
                  + pb_z[k] * hg_27[k];

        t_40[k] = f_5 * gg_14[k]
                  + f_2 * hh_s_40[k]
                  + pb_y[k] * hg_29[k];

        t_41[k] = pa_y[k] * gh_20[k]
                  + f_2 * hh_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_z, pb_y, pb_z, gg_0, gh_0, gh_3, hh_s_42, \
                         hh_s_43, hh_s_44, hh_s_45, hg_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_z[k] * gh_0[k]
                  + f_2 * hh_s_42[k];

        t_43[k] = f_2 * hh_s_43[k]
                  + pb_y[k] * hg_30[k];

        t_44[k] = f_5 * gg_0[k]
                  + f_2 * hh_s_44[k]
                  + pb_z[k] * hg_30[k];

        t_45[k] = pa_z[k] * gh_3[k]
                  + f_2 * hh_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, t_49, pa_z, pb_y, gg_2, gg_3, gh_5, gh_6, gh_7, \
                         hh_s_46, hh_s_47, hh_s_48, hh_s_49, hg_32 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_2 * hh_s_46[k]
                  + pb_y[k] * hg_32[k];

        t_47[k] = f_7 * gg_2[k]
                  + pa_z[k] * gh_5[k]
                  + f_2 * hh_s_47[k];

        t_48[k] = pa_z[k] * gh_6[k]
                  + f_2 * hh_s_48[k];

        t_49[k] = f_5 * gg_3[k]
                  + pa_z[k] * gh_7[k]
                  + f_2 * hh_s_49[k];
    }

#pragma omp simd aligned(t_50, t_51, t_52, pa_z, pb_y, gg_5, gh_9, gh_10, hh_s_50, hh_s_51, \
                         hh_s_52, hg_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_50[k] = f_2 * hh_s_50[k]
                  + pb_y[k] * hg_35[k];

        t_51[k] = f_8 * gg_5[k]
                  + pa_z[k] * gh_9[k]
                  + f_2 * hh_s_51[k];

        t_52[k] = pa_z[k] * gh_10[k]
                  + f_2 * hh_s_52[k];
    }

#pragma omp simd aligned(t_53, t_54, t_55, pb_x, pb_y, gg_41, gg_42, hh_s_53, hh_s_54, \
                         hh_s_55, hg_39, hg_41, hg_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_53[k] = f_3 * gg_41[k]
                  + f_2 * hh_s_53[k]
                  + pb_x[k] * hg_41[k];

        t_54[k] = f_3 * gg_42[k]
                  + f_2 * hh_s_54[k]
                  + pb_x[k] * hg_42[k];

        t_55[k] = f_2 * hh_s_55[k]
                  + pb_y[k] * hg_39[k];
    }

#pragma omp simd aligned(t_56, t_57, t_58, pa_z, pb_x, pb_y, gg_44, gh_15, hf_s_27, hh_s_56, \
                         hh_s_57, hh_s_58, hf_27, hg_41, hg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_56[k] = f_3 * gg_44[k]
                  + f_2 * hh_s_56[k]
                  + pb_x[k] * hg_44[k];

        t_57[k] = pa_z[k] * gh_15[k]
                  + f_2 * hh_s_57[k];

        t_58[k] = -f_10 * hf_s_27[k]
                  + f_2 * hh_s_58[k]
                  + f_8 * hf_27[k]
                  + pb_y[k] * hg_41[k];
    }

#pragma omp simd aligned(t_59, t_60, t_61, pb_y, hf_s_28, hf_s_29, hh_s_59, hh_s_60, hh_s_61, \
                         hf_28, hf_29, hg_42, hg_43, hg_44 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_59[k] = -f_6 * hf_s_28[k]
                  + f_2 * hh_s_59[k]
                  + f_7 * hf_28[k]
                  + pb_y[k] * hg_42[k];

        t_60[k] = -f_4 * hf_s_29[k]
                  + f_2 * hh_s_60[k]
                  + f_5 * hf_29[k]
                  + pb_y[k] * hg_43[k];

        t_61[k] = f_2 * hh_s_61[k]
                  + pb_y[k] * hg_44[k];
    }

#pragma omp simd aligned(t_62, t_63, pa_x, pa_y, fh_s_0, fh_s_62, fh_0, fh_62, gh_21, gh_62, \
                         hh_s_62, hh_s_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_62[k] = -f_9 * fh_s_62[k]
                  + f_8 * fh_62[k]
                  + pa_x[k] * gh_62[k]
                  + f_2 * hh_s_62[k];

        t_63[k] = -f_11 * fh_s_0[k]
                  + f_5 * fh_0[k]
                  + pa_y[k] * gh_21[k]
                  + f_2 * hh_s_63[k];
    }

#pragma omp simd aligned(t_64, t_65, t_66, pb_x, pb_y, pb_z, gg_15, gg_48, hf_s_33, hh_s_64, \
                         hh_s_65, hh_s_66, hf_33, hg_45, hg_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_64[k] = f_7 * gg_15[k]
                  + f_2 * hh_s_64[k]
                  + pb_y[k] * hg_45[k];

        t_65[k] = f_2 * hh_s_65[k]
                  + pb_z[k] * hg_45[k];

        t_66[k] = f_8 * gg_48[k]
                  - f_6 * hf_s_33[k]
                  + f_2 * hh_s_66[k]
                  + f_7 * hf_33[k]
                  + pb_x[k] * hg_48[k];
    }

#pragma omp simd aligned(t_67, t_68, t_69, pb_x, pb_z, gg_51, hf_s_30, hf_s_36, hh_s_67, \
                         hh_s_68, hh_s_69, hf_30, hf_36, hg_46, hg_47, \
                         hg_51 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_67[k] = f_2 * hh_s_67[k]
                  + pb_z[k] * hg_46[k];

        t_68[k] = -f_4 * hf_s_30[k]
                  + f_2 * hh_s_68[k]
                  + f_5 * hf_30[k]
                  + pb_z[k] * hg_47[k];

        t_69[k] = f_8 * gg_51[k]
                  - f_4 * hf_s_36[k]
                  + f_2 * hh_s_69[k]
                  + f_5 * hf_36[k]
                  + pb_x[k] * hg_51[k];
    }

#pragma omp simd aligned(t_70, t_71, t_72, pb_y, pb_z, gg_20, hf_s_32, hh_s_70, hh_s_71, \
                         hh_s_72, hf_32, hg_48, hg_50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_70[k] = f_2 * hh_s_70[k]
                  + pb_z[k] * hg_48[k];

        t_71[k] = f_7 * gg_20[k]
                  + f_2 * hh_s_71[k]
                  + pb_y[k] * hg_50[k];

        t_72[k] = -f_6 * hf_s_32[k]
                  + f_2 * hh_s_72[k]
                  + f_7 * hf_32[k]
                  + pb_z[k] * hg_50[k];
    }

#pragma omp simd aligned(t_73, t_74, t_75, pb_x, pb_z, gg_55, gg_57, hh_s_73, hh_s_74, \
                         hh_s_75, hg_51, hg_55, hg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_73[k] = f_8 * gg_55[k]
                  + f_2 * hh_s_73[k]
                  + pb_x[k] * hg_55[k];

        t_74[k] = f_2 * hh_s_74[k]
                  + pb_z[k] * hg_51[k];

        t_75[k] = f_8 * gg_57[k]
                  + f_2 * hh_s_75[k]
                  + pb_x[k] * hg_57[k];
    }

#pragma omp simd aligned(t_76, t_77, t_78, pa_x, pb_x, fh_s_78, fh_78, gg_58, gg_59, gh_78, \
                         hh_s_76, hh_s_77, hh_s_78, hg_58, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_76[k] = f_8 * gg_58[k]
                  + f_2 * hh_s_76[k]
                  + pb_x[k] * hg_58[k];

        t_77[k] = f_8 * gg_59[k]
                  + f_2 * hh_s_77[k]
                  + pb_x[k] * hg_59[k];

        t_78[k] = -f_12 * fh_s_78[k]
                  + f_7 * fh_78[k]
                  + pa_x[k] * gh_78[k]
                  + f_2 * hh_s_78[k];
    }

#pragma omp simd aligned(t_79, t_80, t_81, pb_z, hf_s_36, hf_s_37, hh_s_79, hh_s_80, hh_s_81, \
                         hf_36, hf_37, hg_55, hg_56, hg_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_79[k] = f_2 * hh_s_79[k]
                  + pb_z[k] * hg_55[k];

        t_80[k] = -f_4 * hf_s_36[k]
                  + f_2 * hh_s_80[k]
                  + f_5 * hf_36[k]
                  + pb_z[k] * hg_56[k];

        t_81[k] = -f_6 * hf_s_37[k]
                  + f_2 * hh_s_81[k]
                  + f_7 * hf_37[k]
                  + pb_z[k] * hg_57[k];
    }

#pragma omp simd aligned(t_82, t_83, t_84, pa_y, pb_y, pb_z, gg_29, gh_42, hf_s_39, hh_s_82, \
                         hh_s_83, hh_s_84, hf_39, hg_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_82[k] = f_7 * gg_29[k]
                  + f_2 * hh_s_82[k]
                  + pb_y[k] * hg_59[k];

        t_83[k] = -f_1 * hf_s_39[k]
                  + f_2 * hh_s_83[k]
                  + f_3 * hf_39[k]
                  + pb_z[k] * hg_59[k];

        t_84[k] = pa_y[k] * gh_42[k]
                  + f_2 * hh_s_84[k];
    }

#pragma omp simd aligned(t_85, t_86, t_87, t_88, pa_y, pa_z, pb_y, gg_32, gh_22, gh_24, gh_44, \
                         hh_s_85, hh_s_86, hh_s_87, hh_s_88, hg_62 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_85[k] = pa_z[k] * gh_22[k]
                  + f_2 * hh_s_85[k];

        t_86[k] = pa_y[k] * gh_44[k]
                  + f_2 * hh_s_86[k];

        t_87[k] = pa_z[k] * gh_24[k]
                  + f_2 * hh_s_87[k];

        t_88[k] = f_5 * gg_32[k]
                  + f_2 * hh_s_88[k]
                  + pb_y[k] * hg_62[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_y, pa_z, pb_z, gg_18, gh_27, gh_47, hh_s_89, \
                         hh_s_90, hh_s_91, hg_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * gh_47[k]
                  + f_2 * hh_s_89[k];

        t_90[k] = pa_z[k] * gh_27[k]
                  + f_2 * hh_s_90[k];

        t_91[k] = f_5 * gg_18[k]
                  + f_2 * hh_s_91[k]
                  + pb_z[k] * hg_63[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pa_y, pa_z, pb_y, gg_35, gh_31, gh_51, hh_s_92, \
                         hh_s_93, hh_s_94, hg_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_5 * gg_35[k]
                  + f_2 * hh_s_92[k]
                  + pb_y[k] * hg_65[k];

        t_93[k] = pa_y[k] * gh_51[k]
                  + f_2 * hh_s_93[k];

        t_94[k] = pa_z[k] * gh_31[k]
                  + f_2 * hh_s_94[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, gg_71, gg_72, gg_73, hh_s_95, hh_s_96, \
                         hh_s_97, hg_71, hg_72, hg_73 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_8 * gg_71[k]
                  + f_2 * hh_s_95[k]
                  + pb_x[k] * hg_71[k];

        t_96[k] = f_8 * gg_72[k]
                  + f_2 * hh_s_96[k]
                  + pb_x[k] * hg_72[k];

        t_97[k] = f_8 * gg_73[k]
                  + f_2 * hh_s_97[k]
                  + pb_x[k] * hg_73[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_y, pa_z, pb_z, gg_25, gh_36, gh_56, hh_s_98, \
                         hh_s_99, hh_s_100, hg_70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = pa_y[k] * gh_56[k]
                  + f_2 * hh_s_98[k];

        t_99[k] = pa_z[k] * gh_36[k]
                  + f_2 * hh_s_99[k];

        t_100[k] = f_5 * gg_25[k]
                   + f_2 * hh_s_100[k]
                   + pb_z[k] * hg_70[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, pa_x, pb_y, fh_s_101, fh_s_102, fh_101, fh_102, \
                         gg_44, gh_101, gh_102, hh_s_101, hh_s_102, hh_s_103, \
                         hg_74 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = -f_12 * fh_s_101[k]
                   + f_7 * fh_101[k]
                   + pa_x[k] * gh_101[k]
                   + f_2 * hh_s_101[k];

        t_102[k] = -f_12 * fh_s_102[k]
                   + f_7 * fh_102[k]
                   + pa_x[k] * gh_102[k]
                   + f_2 * hh_s_102[k];

        t_103[k] = f_5 * gg_44[k]
                   + f_2 * hh_s_103[k]
                   + pb_y[k] * hg_74[k];
    }

#pragma omp simd aligned(t_104, t_105, t_106, pa_y, pa_z, pb_y, fh_s_0, fh_0, gh_42, gh_62, \
                         hh_s_104, hh_s_105, hh_s_106, hg_75 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_104[k] = pa_y[k] * gh_62[k]
                   + f_2 * hh_s_104[k];

        t_105[k] = -f_11 * fh_s_0[k]
                   + f_5 * fh_0[k]
                   + pa_z[k] * gh_42[k]
                   + f_2 * hh_s_105[k];

        t_106[k] = f_2 * hh_s_106[k]
                   + pb_y[k] * hg_75[k];
    }
}

static auto
compute_prim_hh_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fh_s,
                                        const size_t fh, const size_t gg, const size_t gh,
                                        const size_t hf_s, const size_t hh_s, const size_t hf,
                                        const size_t hg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_10 = 3.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_s_21 = buffer.data(fh_s + 21);
    const auto *fh_s_42 = buffer.data(fh_s + 42);
    const auto *fh_s_47 = buffer.data(fh_s + 47);
    const auto *fh_s_51 = buffer.data(fh_s + 51);
    const auto *fh_s_125 = buffer.data(fh_s + 125);
    const auto *fh_s_141 = buffer.data(fh_s + 141);
    const auto *fh_s_164 = buffer.data(fh_s + 164);
    const auto *fh_s_165 = buffer.data(fh_s + 165);
    const auto *fh_s_167 = buffer.data(fh_s + 167);
    const auto *fh_s_183 = buffer.data(fh_s + 183);
    const auto *fh_s_185 = buffer.data(fh_s + 185);
    const auto *fh_s_186 = buffer.data(fh_s + 186);

    const auto *fh_21 = buffer.data(fh + 21);
    const auto *fh_42 = buffer.data(fh + 42);
    const auto *fh_47 = buffer.data(fh + 47);
    const auto *fh_51 = buffer.data(fh + 51);
    const auto *fh_125 = buffer.data(fh + 125);
    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_164 = buffer.data(fh + 164);
    const auto *fh_165 = buffer.data(fh + 165);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_183 = buffer.data(fh + 183);
    const auto *fh_185 = buffer.data(fh + 185);
    const auto *fh_186 = buffer.data(fh + 186);

    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_50 = buffer.data(gg + 50);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_59 = buffer.data(gg + 59);
    const auto *gg_62 = buffer.data(gg + 62);
    const auto *gg_63 = buffer.data(gg + 63);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_70 = buffer.data(gg + 70);
    const auto *gg_74 = buffer.data(gg + 74);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_76 = buffer.data(gg + 76);
    const auto *gg_77 = buffer.data(gg + 77);
    const auto *gg_78 = buffer.data(gg + 78);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_84 = buffer.data(gg + 84);
    const auto *gg_85 = buffer.data(gg + 85);
    const auto *gg_86 = buffer.data(gg + 86);
    const auto *gg_87 = buffer.data(gg + 87);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_96 = buffer.data(gg + 96);
    const auto *gg_100 = buffer.data(gg + 100);
    const auto *gg_102 = buffer.data(gg + 102);
    const auto *gg_103 = buffer.data(gg + 103);
    const auto *gg_104 = buffer.data(gg + 104);
    const auto *gg_116 = buffer.data(gg + 116);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_118 = buffer.data(gg + 118);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_131 = buffer.data(gg + 131);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_133 = buffer.data(gg + 133);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_144 = buffer.data(gg + 144);
    const auto *gg_145 = buffer.data(gg + 145);
    const auto *gg_146 = buffer.data(gg + 146);
    const auto *gg_147 = buffer.data(gg + 147);
    const auto *gg_149 = buffer.data(gg + 149);

    const auto *gh_63 = buffer.data(gh + 63);
    const auto *gh_64 = buffer.data(gh + 64);
    const auto *gh_66 = buffer.data(gh + 66);
    const auto *gh_69 = buffer.data(gh + 69);
    const auto *gh_73 = buffer.data(gh + 73);
    const auto *gh_78 = buffer.data(gh + 78);
    const auto *gh_89 = buffer.data(gh + 89);
    const auto *gh_93 = buffer.data(gh + 93);
    const auto *gh_105 = buffer.data(gh + 105);
    const auto *gh_107 = buffer.data(gh + 107);
    const auto *gh_108 = buffer.data(gh + 108);
    const auto *gh_110 = buffer.data(gh + 110);
    const auto *gh_111 = buffer.data(gh + 111);
    const auto *gh_114 = buffer.data(gh + 114);
    const auto *gh_119 = buffer.data(gh + 119);
    const auto *gh_125 = buffer.data(gh + 125);
    const auto *gh_141 = buffer.data(gh + 141);
    const auto *gh_164 = buffer.data(gh + 164);
    const auto *gh_165 = buffer.data(gh + 165);
    const auto *gh_167 = buffer.data(gh + 167);
    const auto *gh_183 = buffer.data(gh + 183);
    const auto *gh_185 = buffer.data(gh + 185);
    const auto *gh_186 = buffer.data(gh + 186);

    const auto *hf_s_50 = buffer.data(hf_s + 50);
    const auto *hf_s_51 = buffer.data(hf_s + 51);
    const auto *hf_s_52 = buffer.data(hf_s + 52);
    const auto *hf_s_55 = buffer.data(hf_s + 55);
    const auto *hf_s_56 = buffer.data(hf_s + 56);
    const auto *hf_s_57 = buffer.data(hf_s + 57);
    const auto *hf_s_58 = buffer.data(hf_s + 58);
    const auto *hf_s_59 = buffer.data(hf_s + 59);
    const auto *hf_s_60 = buffer.data(hf_s + 60);
    const auto *hf_s_62 = buffer.data(hf_s + 62);
    const auto *hf_s_63 = buffer.data(hf_s + 63);
    const auto *hf_s_66 = buffer.data(hf_s + 66);
    const auto *hf_s_67 = buffer.data(hf_s + 67);
    const auto *hf_s_69 = buffer.data(hf_s + 69);
    const auto *hf_s_90 = buffer.data(hf_s + 90);
    const auto *hf_s_91 = buffer.data(hf_s + 91);
    const auto *hf_s_92 = buffer.data(hf_s + 92);
    const auto *hf_s_95 = buffer.data(hf_s + 95);
    const auto *hf_s_96 = buffer.data(hf_s + 96);
    const auto *hf_s_99 = buffer.data(hf_s + 99);

    const auto *hh_s_107 = buffer.data(hh_s + 107);
    const auto *hh_s_108 = buffer.data(hh_s + 108);
    const auto *hh_s_109 = buffer.data(hh_s + 109);
    const auto *hh_s_110 = buffer.data(hh_s + 110);
    const auto *hh_s_111 = buffer.data(hh_s + 111);
    const auto *hh_s_112 = buffer.data(hh_s + 112);
    const auto *hh_s_113 = buffer.data(hh_s + 113);
    const auto *hh_s_114 = buffer.data(hh_s + 114);
    const auto *hh_s_115 = buffer.data(hh_s + 115);
    const auto *hh_s_116 = buffer.data(hh_s + 116);
    const auto *hh_s_117 = buffer.data(hh_s + 117);
    const auto *hh_s_118 = buffer.data(hh_s + 118);
    const auto *hh_s_119 = buffer.data(hh_s + 119);
    const auto *hh_s_120 = buffer.data(hh_s + 120);
    const auto *hh_s_121 = buffer.data(hh_s + 121);
    const auto *hh_s_122 = buffer.data(hh_s + 122);
    const auto *hh_s_123 = buffer.data(hh_s + 123);
    const auto *hh_s_124 = buffer.data(hh_s + 124);
    const auto *hh_s_125 = buffer.data(hh_s + 125);
    const auto *hh_s_126 = buffer.data(hh_s + 126);
    const auto *hh_s_127 = buffer.data(hh_s + 127);
    const auto *hh_s_128 = buffer.data(hh_s + 128);
    const auto *hh_s_129 = buffer.data(hh_s + 129);
    const auto *hh_s_130 = buffer.data(hh_s + 130);
    const auto *hh_s_131 = buffer.data(hh_s + 131);
    const auto *hh_s_132 = buffer.data(hh_s + 132);
    const auto *hh_s_133 = buffer.data(hh_s + 133);
    const auto *hh_s_134 = buffer.data(hh_s + 134);
    const auto *hh_s_135 = buffer.data(hh_s + 135);
    const auto *hh_s_136 = buffer.data(hh_s + 136);
    const auto *hh_s_137 = buffer.data(hh_s + 137);
    const auto *hh_s_138 = buffer.data(hh_s + 138);
    const auto *hh_s_139 = buffer.data(hh_s + 139);
    const auto *hh_s_140 = buffer.data(hh_s + 140);
    const auto *hh_s_141 = buffer.data(hh_s + 141);
    const auto *hh_s_142 = buffer.data(hh_s + 142);
    const auto *hh_s_143 = buffer.data(hh_s + 143);
    const auto *hh_s_144 = buffer.data(hh_s + 144);
    const auto *hh_s_145 = buffer.data(hh_s + 145);
    const auto *hh_s_146 = buffer.data(hh_s + 146);
    const auto *hh_s_147 = buffer.data(hh_s + 147);
    const auto *hh_s_148 = buffer.data(hh_s + 148);
    const auto *hh_s_149 = buffer.data(hh_s + 149);
    const auto *hh_s_150 = buffer.data(hh_s + 150);
    const auto *hh_s_151 = buffer.data(hh_s + 151);
    const auto *hh_s_152 = buffer.data(hh_s + 152);
    const auto *hh_s_153 = buffer.data(hh_s + 153);
    const auto *hh_s_154 = buffer.data(hh_s + 154);
    const auto *hh_s_155 = buffer.data(hh_s + 155);
    const auto *hh_s_156 = buffer.data(hh_s + 156);
    const auto *hh_s_157 = buffer.data(hh_s + 157);
    const auto *hh_s_158 = buffer.data(hh_s + 158);
    const auto *hh_s_159 = buffer.data(hh_s + 159);
    const auto *hh_s_160 = buffer.data(hh_s + 160);
    const auto *hh_s_161 = buffer.data(hh_s + 161);
    const auto *hh_s_162 = buffer.data(hh_s + 162);
    const auto *hh_s_163 = buffer.data(hh_s + 163);
    const auto *hh_s_164 = buffer.data(hh_s + 164);
    const auto *hh_s_165 = buffer.data(hh_s + 165);
    const auto *hh_s_166 = buffer.data(hh_s + 166);
    const auto *hh_s_167 = buffer.data(hh_s + 167);
    const auto *hh_s_168 = buffer.data(hh_s + 168);
    const auto *hh_s_169 = buffer.data(hh_s + 169);
    const auto *hh_s_170 = buffer.data(hh_s + 170);
    const auto *hh_s_171 = buffer.data(hh_s + 171);
    const auto *hh_s_172 = buffer.data(hh_s + 172);
    const auto *hh_s_173 = buffer.data(hh_s + 173);
    const auto *hh_s_174 = buffer.data(hh_s + 174);
    const auto *hh_s_175 = buffer.data(hh_s + 175);
    const auto *hh_s_176 = buffer.data(hh_s + 176);
    const auto *hh_s_177 = buffer.data(hh_s + 177);
    const auto *hh_s_178 = buffer.data(hh_s + 178);
    const auto *hh_s_179 = buffer.data(hh_s + 179);
    const auto *hh_s_180 = buffer.data(hh_s + 180);
    const auto *hh_s_181 = buffer.data(hh_s + 181);
    const auto *hh_s_182 = buffer.data(hh_s + 182);
    const auto *hh_s_183 = buffer.data(hh_s + 183);
    const auto *hh_s_184 = buffer.data(hh_s + 184);
    const auto *hh_s_185 = buffer.data(hh_s + 185);
    const auto *hh_s_186 = buffer.data(hh_s + 186);
    const auto *hh_s_187 = buffer.data(hh_s + 187);
    const auto *hh_s_188 = buffer.data(hh_s + 188);
    const auto *hh_s_189 = buffer.data(hh_s + 189);
    const auto *hh_s_190 = buffer.data(hh_s + 190);
    const auto *hh_s_191 = buffer.data(hh_s + 191);
    const auto *hh_s_192 = buffer.data(hh_s + 192);
    const auto *hh_s_193 = buffer.data(hh_s + 193);
    const auto *hh_s_194 = buffer.data(hh_s + 194);
    const auto *hh_s_195 = buffer.data(hh_s + 195);
    const auto *hh_s_196 = buffer.data(hh_s + 196);
    const auto *hh_s_197 = buffer.data(hh_s + 197);
    const auto *hh_s_198 = buffer.data(hh_s + 198);
    const auto *hh_s_199 = buffer.data(hh_s + 199);
    const auto *hh_s_200 = buffer.data(hh_s + 200);
    const auto *hh_s_201 = buffer.data(hh_s + 201);
    const auto *hh_s_202 = buffer.data(hh_s + 202);
    const auto *hh_s_203 = buffer.data(hh_s + 203);
    const auto *hh_s_204 = buffer.data(hh_s + 204);

    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_58 = buffer.data(hf + 58);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_62 = buffer.data(hf + 62);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_67 = buffer.data(hf + 67);
    const auto *hf_69 = buffer.data(hf + 69);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_91 = buffer.data(hf + 91);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_96 = buffer.data(hf + 96);
    const auto *hf_99 = buffer.data(hf + 99);

    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_76 = buffer.data(hg + 76);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_79 = buffer.data(hg + 79);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_85 = buffer.data(hg + 85);
    const auto *hg_86 = buffer.data(hg + 86);
    const auto *hg_87 = buffer.data(hg + 87);
    const auto *hg_88 = buffer.data(hg + 88);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_92 = buffer.data(hg + 92);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_95 = buffer.data(hg + 95);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_101 = buffer.data(hg + 101);
    const auto *hg_102 = buffer.data(hg + 102);
    const auto *hg_103 = buffer.data(hg + 103);
    const auto *hg_104 = buffer.data(hg + 104);
    const auto *hg_105 = buffer.data(hg + 105);
    const auto *hg_107 = buffer.data(hg + 107);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_115 = buffer.data(hg + 115);
    const auto *hg_116 = buffer.data(hg + 116);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_118 = buffer.data(hg + 118);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_122 = buffer.data(hg + 122);
    const auto *hg_123 = buffer.data(hg + 123);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_131 = buffer.data(hg + 131);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_133 = buffer.data(hg + 133);
    const auto *hg_134 = buffer.data(hg + 134);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_136 = buffer.data(hg + 136);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_139 = buffer.data(hg + 139);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_145 = buffer.data(hg + 145);
    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_149 = buffer.data(hg + 149);

#pragma omp simd aligned(t_107, t_108, t_109, pb_y, pb_z, gg_30, hf_s_50, hh_s_107, hh_s_108, \
                         hh_s_109, hf_50, hg_75, hg_76, hg_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_107[k] = f_7 * gg_30[k]
                   + f_2 * hh_s_107[k]
                   + pb_z[k] * hg_75[k];

        t_108[k] = -f_4 * hf_s_50[k]
                   + f_2 * hh_s_108[k]
                   + f_5 * hf_50[k]
                   + pb_y[k] * hg_76[k];

        t_109[k] = f_2 * hh_s_109[k]
                   + pb_y[k] * hg_77[k];
    }

#pragma omp simd aligned(t_110, t_111, pb_x, pb_y, gg_80, hf_s_51, hf_s_55, hh_s_110, \
                         hh_s_111, hf_51, hf_55, hg_78, hg_80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_110[k] = f_8 * gg_80[k]
                   - f_6 * hf_s_55[k]
                   + f_2 * hh_s_110[k]
                   + f_7 * hf_55[k]
                   + pb_x[k] * hg_80[k];

        t_111[k] = -f_6 * hf_s_51[k]
                   + f_2 * hh_s_111[k]
                   + f_7 * hf_51[k]
                   + pb_y[k] * hg_78[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, pb_y, gg_84, hf_s_52, hf_s_59, hh_s_112, \
                         hh_s_113, hh_s_114, hf_52, hf_59, hg_79, hg_80, \
                         hg_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = -f_4 * hf_s_52[k]
                   + f_2 * hh_s_112[k]
                   + f_5 * hf_52[k]
                   + pb_y[k] * hg_79[k];

        t_113[k] = f_2 * hh_s_113[k]
                   + pb_y[k] * hg_80[k];

        t_114[k] = f_8 * gg_84[k]
                   - f_4 * hf_s_59[k]
                   + f_2 * hh_s_114[k]
                   + f_5 * hf_59[k]
                   + pb_x[k] * hg_84[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pb_x, gg_85, gg_86, gg_87, hh_s_115, hh_s_116, \
                         hh_s_117, hg_85, hg_86, hg_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = f_8 * gg_85[k]
                   + f_2 * hh_s_115[k]
                   + pb_x[k] * hg_85[k];

        t_116[k] = f_8 * gg_86[k]
                   + f_2 * hh_s_116[k]
                   + pb_x[k] * hg_86[k];

        t_117[k] = f_8 * gg_87[k]
                   + f_2 * hh_s_117[k]
                   + pb_x[k] * hg_87[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pb_x, pb_y, gg_89, hf_s_56, hh_s_118, hh_s_119, \
                         hh_s_120, hf_56, hg_84, hg_85, hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_2 * hh_s_118[k]
                   + pb_y[k] * hg_84[k];

        t_119[k] = f_8 * gg_89[k]
                   + f_2 * hh_s_119[k]
                   + pb_x[k] * hg_89[k];

        t_120[k] = -f_1 * hf_s_56[k]
                   + f_2 * hh_s_120[k]
                   + f_3 * hf_56[k]
                   + pb_y[k] * hg_85[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pb_y, hf_s_57, hf_s_58, hf_s_59, hh_s_121, \
                         hh_s_122, hh_s_123, hf_57, hf_58, hf_59, hg_86, hg_87, \
                         hg_88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = -f_10 * hf_s_57[k]
                   + f_2 * hh_s_121[k]
                   + f_8 * hf_57[k]
                   + pb_y[k] * hg_86[k];

        t_122[k] = -f_6 * hf_s_58[k]
                   + f_2 * hh_s_122[k]
                   + f_7 * hf_58[k]
                   + pb_y[k] * hg_87[k];

        t_123[k] = -f_4 * hf_s_59[k]
                   + f_2 * hh_s_123[k]
                   + f_5 * hf_59[k]
                   + pb_y[k] * hg_88[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_x, pa_y, pb_y, fh_s_21, fh_s_125, fh_21, \
                         fh_125, gh_63, gh_125, hh_s_124, hh_s_125, hh_s_126, \
                         hg_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_2 * hh_s_124[k]
                   + pb_y[k] * hg_89[k];

        t_125[k] = -f_12 * fh_s_125[k]
                   + f_7 * fh_125[k]
                   + pa_x[k] * gh_125[k]
                   + f_2 * hh_s_125[k];

        t_126[k] = -f_12 * fh_s_21[k]
                   + f_7 * fh_21[k]
                   + pa_y[k] * gh_63[k]
                   + f_2 * hh_s_126[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pb_x, pb_y, pb_z, gg_45, gg_93, hf_s_63, \
                         hh_s_127, hh_s_128, hh_s_129, hf_63, hg_90, \
                         hg_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_8 * gg_45[k]
                   + f_2 * hh_s_127[k]
                   + pb_y[k] * hg_90[k];

        t_128[k] = f_2 * hh_s_128[k]
                   + pb_z[k] * hg_90[k];

        t_129[k] = f_7 * gg_93[k]
                   - f_6 * hf_s_63[k]
                   + f_2 * hh_s_129[k]
                   + f_7 * hf_63[k]
                   + pb_x[k] * hg_93[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pb_x, pb_z, gg_96, hf_s_60, hf_s_66, hh_s_130, \
                         hh_s_131, hh_s_132, hf_60, hf_66, hg_91, hg_92, \
                         hg_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = f_2 * hh_s_130[k]
                   + pb_z[k] * hg_91[k];

        t_131[k] = -f_4 * hf_s_60[k]
                   + f_2 * hh_s_131[k]
                   + f_5 * hf_60[k]
                   + pb_z[k] * hg_92[k];

        t_132[k] = f_7 * gg_96[k]
                   - f_4 * hf_s_66[k]
                   + f_2 * hh_s_132[k]
                   + f_5 * hf_66[k]
                   + pb_x[k] * hg_96[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pb_y, pb_z, gg_50, hf_s_62, hh_s_133, hh_s_134, \
                         hh_s_135, hf_62, hg_93, hg_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_2 * hh_s_133[k]
                   + pb_z[k] * hg_93[k];

        t_134[k] = f_8 * gg_50[k]
                   + f_2 * hh_s_134[k]
                   + pb_y[k] * hg_95[k];

        t_135[k] = -f_6 * hf_s_62[k]
                   + f_2 * hh_s_135[k]
                   + f_7 * hf_62[k]
                   + pb_z[k] * hg_95[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pb_x, pb_z, gg_100, gg_102, hh_s_136, hh_s_137, \
                         hh_s_138, hg_96, hg_100, hg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_7 * gg_100[k]
                   + f_2 * hh_s_136[k]
                   + pb_x[k] * hg_100[k];

        t_137[k] = f_2 * hh_s_137[k]
                   + pb_z[k] * hg_96[k];

        t_138[k] = f_7 * gg_102[k]
                   + f_2 * hh_s_138[k]
                   + pb_x[k] * hg_102[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_x, pb_x, fh_s_141, fh_141, gg_103, gg_104, \
                         gh_141, hh_s_139, hh_s_140, hh_s_141, hg_103, \
                         hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = f_7 * gg_103[k]
                   + f_2 * hh_s_139[k]
                   + pb_x[k] * hg_103[k];

        t_140[k] = f_7 * gg_104[k]
                   + f_2 * hh_s_140[k]
                   + pb_x[k] * hg_104[k];

        t_141[k] = -f_11 * fh_s_141[k]
                   + f_5 * fh_141[k]
                   + pa_x[k] * gh_141[k]
                   + f_2 * hh_s_141[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_z, hf_s_66, hf_s_67, hh_s_142, hh_s_143, \
                         hh_s_144, hf_66, hf_67, hg_100, hg_101, \
                         hg_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_2 * hh_s_142[k]
                   + pb_z[k] * hg_100[k];

        t_143[k] = -f_4 * hf_s_66[k]
                   + f_2 * hh_s_143[k]
                   + f_5 * hf_66[k]
                   + pb_z[k] * hg_101[k];

        t_144[k] = -f_6 * hf_s_67[k]
                   + f_2 * hh_s_144[k]
                   + f_7 * hf_67[k]
                   + pb_z[k] * hg_102[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pa_z, pb_y, pb_z, gg_59, gh_63, hf_s_69, \
                         hh_s_145, hh_s_146, hh_s_147, hf_69, hg_104 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_8 * gg_59[k]
                   + f_2 * hh_s_145[k]
                   + pb_y[k] * hg_104[k];

        t_146[k] = -f_1 * hf_s_69[k]
                   + f_2 * hh_s_146[k]
                   + f_3 * hf_69[k]
                   + pb_z[k] * hg_104[k];

        t_147[k] = pa_z[k] * gh_63[k]
                   + f_2 * hh_s_147[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_z, pb_z, gg_45, gh_64, gh_66, hh_s_148, \
                         hh_s_149, hh_s_150, hg_105 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = pa_z[k] * gh_64[k]
                   + f_2 * hh_s_148[k];

        t_149[k] = f_5 * gg_45[k]
                   + f_2 * hh_s_149[k]
                   + pb_z[k] * hg_105[k];

        t_150[k] = pa_z[k] * gh_66[k]
                   + f_2 * hh_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pa_y, pa_z, pb_y, fh_s_47, fh_47, gg_62, gh_69, \
                         gh_89, hh_s_151, hh_s_152, hh_s_153, hg_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_7 * gg_62[k]
                   + f_2 * hh_s_151[k]
                   + pb_y[k] * hg_107[k];

        t_152[k] = -f_11 * fh_s_47[k]
                   + f_5 * fh_47[k]
                   + pa_y[k] * gh_89[k]
                   + f_2 * hh_s_152[k];

        t_153[k] = pa_z[k] * gh_69[k]
                   + f_2 * hh_s_153[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pa_y, pb_y, pb_z, fh_s_51, fh_51, gg_48, gg_65, \
                         gh_93, hh_s_154, hh_s_155, hh_s_156, hg_108, \
                         hg_110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_5 * gg_48[k]
                   + f_2 * hh_s_154[k]
                   + pb_z[k] * hg_108[k];

        t_155[k] = f_7 * gg_65[k]
                   + f_2 * hh_s_155[k]
                   + pb_y[k] * hg_110[k];

        t_156[k] = -f_11 * fh_s_51[k]
                   + f_5 * fh_51[k]
                   + pa_y[k] * gh_93[k]
                   + f_2 * hh_s_156[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pa_z, pb_x, gg_116, gg_117, gh_73, hh_s_157, \
                         hh_s_158, hh_s_159, hg_116, hg_117 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = pa_z[k] * gh_73[k]
                   + f_2 * hh_s_157[k];

        t_158[k] = f_7 * gg_116[k]
                   + f_2 * hh_s_158[k]
                   + pb_x[k] * hg_116[k];

        t_159[k] = f_7 * gg_117[k]
                   + f_2 * hh_s_159[k]
                   + pb_x[k] * hg_117[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_z, pb_x, gg_118, gg_119, gh_78, hh_s_160, \
                         hh_s_161, hh_s_162, hg_118, hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = f_7 * gg_118[k]
                   + f_2 * hh_s_160[k]
                   + pb_x[k] * hg_118[k];

        t_161[k] = f_7 * gg_119[k]
                   + f_2 * hh_s_161[k]
                   + pb_x[k] * hg_119[k];

        t_162[k] = pa_z[k] * gh_78[k]
                   + f_2 * hh_s_162[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_x, pb_z, fh_s_164, fh_s_165, fh_164, fh_165, \
                         gg_55, gh_164, gh_165, hh_s_163, hh_s_164, hh_s_165, \
                         hg_115 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_5 * gg_55[k]
                   + f_2 * hh_s_163[k]
                   + pb_z[k] * hg_115[k];

        t_164[k] = -f_11 * fh_s_164[k]
                   + f_5 * fh_164[k]
                   + pa_x[k] * gh_164[k]
                   + f_2 * hh_s_164[k];

        t_165[k] = -f_11 * fh_s_165[k]
                   + f_5 * fh_165[k]
                   + pa_x[k] * gh_165[k]
                   + f_2 * hh_s_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_x, pa_y, pb_y, fh_s_167, fh_167, gg_74, \
                         gh_105, gh_167, hh_s_166, hh_s_167, hh_s_168, \
                         hg_119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = f_7 * gg_74[k]
                   + f_2 * hh_s_166[k]
                   + pb_y[k] * hg_119[k];

        t_167[k] = -f_11 * fh_s_167[k]
                   + f_5 * fh_167[k]
                   + pa_x[k] * gh_167[k]
                   + f_2 * hh_s_167[k];

        t_168[k] = pa_y[k] * gh_105[k]
                   + f_2 * hh_s_168[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pa_y, pb_y, gg_75, gg_76, gh_107, gh_108, \
                         hh_s_169, hh_s_170, hh_s_171, hg_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_5 * gg_75[k]
                   + f_2 * hh_s_169[k]
                   + pb_y[k] * hg_120[k];

        t_170[k] = pa_y[k] * gh_107[k]
                   + f_2 * hh_s_170[k];

        t_171[k] = f_7 * gg_76[k]
                   + pa_y[k] * gh_108[k]
                   + f_2 * hh_s_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pa_y, pb_y, gg_77, gg_78, gh_110, gh_111, \
                         hh_s_172, hh_s_173, hh_s_174, hg_122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_5 * gg_77[k]
                   + f_2 * hh_s_172[k]
                   + pb_y[k] * hg_122[k];

        t_173[k] = pa_y[k] * gh_110[k]
                   + f_2 * hh_s_173[k];

        t_174[k] = f_8 * gg_78[k]
                   + pa_y[k] * gh_111[k]
                   + f_2 * hh_s_174[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_y, pb_y, pb_z, gg_63, gg_80, gh_114, \
                         hh_s_175, hh_s_176, hh_s_177, hg_123, hg_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = f_7 * gg_63[k]
                   + f_2 * hh_s_175[k]
                   + pb_z[k] * hg_123[k];

        t_176[k] = f_5 * gg_80[k]
                   + f_2 * hh_s_176[k]
                   + pb_y[k] * hg_125[k];

        t_177[k] = pa_y[k] * gh_114[k]
                   + f_2 * hh_s_177[k];
    }

#pragma omp simd aligned(t_178, t_179, t_180, pb_x, gg_130, gg_131, gg_132, hh_s_178, \
                         hh_s_179, hh_s_180, hg_130, hg_131, hg_132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_7 * gg_130[k]
                   + f_2 * hh_s_178[k]
                   + pb_x[k] * hg_130[k];

        t_179[k] = f_7 * gg_131[k]
                   + f_2 * hh_s_179[k]
                   + pb_x[k] * hg_131[k];

        t_180[k] = f_7 * gg_132[k]
                   + f_2 * hh_s_180[k]
                   + pb_x[k] * hg_132[k];
    }

#pragma omp simd aligned(t_181, t_182, t_183, pa_x, pa_y, pb_x, fh_s_183, fh_183, gg_133, \
                         gh_119, gh_183, hh_s_181, hh_s_182, hh_s_183, \
                         hg_133 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_181[k] = f_7 * gg_133[k]
                   + f_2 * hh_s_181[k]
                   + pb_x[k] * hg_133[k];

        t_182[k] = pa_y[k] * gh_119[k]
                   + f_2 * hh_s_182[k];

        t_183[k] = -f_11 * fh_s_183[k]
                   + f_5 * fh_183[k]
                   + pa_x[k] * gh_183[k]
                   + f_2 * hh_s_183[k];
    }

#pragma omp simd aligned(t_184, t_185, t_186, pa_x, pb_z, fh_s_185, fh_s_186, fh_185, fh_186, \
                         gg_70, gh_185, gh_186, hh_s_184, hh_s_185, hh_s_186, \
                         hg_130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_184[k] = f_7 * gg_70[k]
                   + f_2 * hh_s_184[k]
                   + pb_z[k] * hg_130[k];

        t_185[k] = -f_11 * fh_s_185[k]
                   + f_5 * fh_185[k]
                   + pa_x[k] * gh_185[k]
                   + f_2 * hh_s_185[k];

        t_186[k] = -f_11 * fh_s_186[k]
                   + f_5 * fh_186[k]
                   + pa_x[k] * gh_186[k]
                   + f_2 * hh_s_186[k];
    }

#pragma omp simd aligned(t_187, t_188, t_189, pa_y, pa_z, pb_y, fh_s_42, fh_42, gg_89, gh_105, \
                         gh_125, hh_s_187, hh_s_188, hh_s_189, hg_134 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_187[k] = f_5 * gg_89[k]
                   + f_2 * hh_s_187[k]
                   + pb_y[k] * hg_134[k];

        t_188[k] = pa_y[k] * gh_125[k]
                   + f_2 * hh_s_188[k];

        t_189[k] = -f_12 * fh_s_42[k]
                   + f_7 * fh_42[k]
                   + pa_z[k] * gh_105[k]
                   + f_2 * hh_s_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, t_193, pb_y, pb_z, gg_75, hf_s_90, hh_s_190, \
                         hh_s_191, hh_s_192, hh_s_193, hf_90, hg_135, hg_136, \
                         hg_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = f_2 * hh_s_190[k]
                   + pb_y[k] * hg_135[k];

        t_191[k] = f_8 * gg_75[k]
                   + f_2 * hh_s_191[k]
                   + pb_z[k] * hg_135[k];

        t_192[k] = -f_4 * hf_s_90[k]
                   + f_2 * hh_s_192[k]
                   + f_5 * hf_90[k]
                   + pb_y[k] * hg_136[k];

        t_193[k] = f_2 * hh_s_193[k]
                   + pb_y[k] * hg_137[k];
    }

#pragma omp simd aligned(t_194, t_195, pb_x, pb_y, gg_140, hf_s_91, hf_s_95, hh_s_194, \
                         hh_s_195, hf_91, hf_95, hg_138, hg_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = f_7 * gg_140[k]
                   - f_6 * hf_s_95[k]
                   + f_2 * hh_s_194[k]
                   + f_7 * hf_95[k]
                   + pb_x[k] * hg_140[k];

        t_195[k] = -f_6 * hf_s_91[k]
                   + f_2 * hh_s_195[k]
                   + f_7 * hf_91[k]
                   + pb_y[k] * hg_138[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, pb_x, pb_y, gg_144, hf_s_92, hf_s_99, hh_s_196, \
                         hh_s_197, hh_s_198, hf_92, hf_99, hg_139, hg_140, \
                         hg_144 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = -f_4 * hf_s_92[k]
                   + f_2 * hh_s_196[k]
                   + f_5 * hf_92[k]
                   + pb_y[k] * hg_139[k];

        t_197[k] = f_2 * hh_s_197[k]
                   + pb_y[k] * hg_140[k];

        t_198[k] = f_7 * gg_144[k]
                   - f_4 * hf_s_99[k]
                   + f_2 * hh_s_198[k]
                   + f_5 * hf_99[k]
                   + pb_x[k] * hg_144[k];
    }

#pragma omp simd aligned(t_199, t_200, t_201, pb_x, gg_145, gg_146, gg_147, hh_s_199, \
                         hh_s_200, hh_s_201, hg_145, hg_146, hg_147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_199[k] = f_7 * gg_145[k]
                   + f_2 * hh_s_199[k]
                   + pb_x[k] * hg_145[k];

        t_200[k] = f_7 * gg_146[k]
                   + f_2 * hh_s_200[k]
                   + pb_x[k] * hg_146[k];

        t_201[k] = f_7 * gg_147[k]
                   + f_2 * hh_s_201[k]
                   + pb_x[k] * hg_147[k];
    }

#pragma omp simd aligned(t_202, t_203, t_204, pb_x, pb_y, gg_149, hf_s_96, hh_s_202, hh_s_203, \
                         hh_s_204, hf_96, hg_144, hg_145, hg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_202[k] = f_2 * hh_s_202[k]
                   + pb_y[k] * hg_144[k];

        t_203[k] = f_7 * gg_149[k]
                   + f_2 * hh_s_203[k]
                   + pb_x[k] * hg_149[k];

        t_204[k] = -f_1 * hf_s_96[k]
                   + f_2 * hh_s_204[k]
                   + f_3 * hf_96[k]
                   + pb_y[k] * hg_145[k];
    }
}

static auto
compute_prim_hh_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fh_s,
                                        const size_t fh, const size_t gg, const size_t gh,
                                        const size_t hf_s, const size_t hh_s, const size_t hf,
                                        const size_t hg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_10 = 3.0 * alpha / p;
    const auto f_11 = beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_s_209 = buffer.data(fh_s + 209);

    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gg_90 = buffer.data(gg + 90);
    const auto *gg_93 = buffer.data(gg + 93);
    const auto *gg_95 = buffer.data(gg + 95);
    const auto *gg_105 = buffer.data(gg + 105);
    const auto *gg_107 = buffer.data(gg + 107);
    const auto *gg_108 = buffer.data(gg + 108);
    const auto *gg_110 = buffer.data(gg + 110);
    const auto *gg_120 = buffer.data(gg + 120);
    const auto *gg_122 = buffer.data(gg + 122);
    const auto *gg_123 = buffer.data(gg + 123);
    const auto *gg_125 = buffer.data(gg + 125);
    const auto *gg_135 = buffer.data(gg + 135);
    const auto *gg_137 = buffer.data(gg + 137);
    const auto *gg_140 = buffer.data(gg + 140);
    const auto *gg_150 = buffer.data(gg + 150);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_155 = buffer.data(gg + 155);
    const auto *gg_156 = buffer.data(gg + 156);
    const auto *gg_159 = buffer.data(gg + 159);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_163 = buffer.data(gg + 163);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_170 = buffer.data(gg + 170);
    const auto *gg_174 = buffer.data(gg + 174);
    const auto *gg_176 = buffer.data(gg + 176);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_178 = buffer.data(gg + 178);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_180 = buffer.data(gg + 180);
    const auto *gg_183 = buffer.data(gg + 183);
    const auto *gg_185 = buffer.data(gg + 185);
    const auto *gg_186 = buffer.data(gg + 186);
    const auto *gg_189 = buffer.data(gg + 189);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_191 = buffer.data(gg + 191);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_198 = buffer.data(gg + 198);
    const auto *gg_201 = buffer.data(gg + 201);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_206 = buffer.data(gg + 206);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_216 = buffer.data(gg + 216);
    const auto *gg_217 = buffer.data(gg + 217);
    const auto *gg_219 = buffer.data(gg + 219);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_221 = buffer.data(gg + 221);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *gh_126 = buffer.data(gh + 126);
    const auto *gh_127 = buffer.data(gh + 127);
    const auto *gh_129 = buffer.data(gh + 129);
    const auto *gh_132 = buffer.data(gh + 132);
    const auto *gh_136 = buffer.data(gh + 136);
    const auto *gh_189 = buffer.data(gh + 189);
    const auto *gh_191 = buffer.data(gh + 191);
    const auto *gh_194 = buffer.data(gh + 194);
    const auto *gh_198 = buffer.data(gh + 198);
    const auto *gh_203 = buffer.data(gh + 203);
    const auto *gh_209 = buffer.data(gh + 209);
    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_215 = buffer.data(gh + 215);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_219 = buffer.data(gh + 219);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_229 = buffer.data(gh + 229);
    const auto *gh_230 = buffer.data(gh + 230);
    const auto *gh_236 = buffer.data(gh + 236);
    const auto *gh_240 = buffer.data(gh + 240);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_247 = buffer.data(gh + 247);
    const auto *gh_248 = buffer.data(gh + 248);
    const auto *gh_249 = buffer.data(gh + 249);
    const auto *gh_250 = buffer.data(gh + 250);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_252 = buffer.data(gh + 252);
    const auto *gh_255 = buffer.data(gh + 255);
    const auto *gh_257 = buffer.data(gh + 257);
    const auto *gh_258 = buffer.data(gh + 258);
    const auto *gh_261 = buffer.data(gh + 261);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_268 = buffer.data(gh + 268);
    const auto *gh_269 = buffer.data(gh + 269);
    const auto *gh_270 = buffer.data(gh + 270);
    const auto *gh_271 = buffer.data(gh + 271);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_276 = buffer.data(gh + 276);
    const auto *gh_279 = buffer.data(gh + 279);
    const auto *gh_288 = buffer.data(gh + 288);
    const auto *gh_289 = buffer.data(gh + 289);
    const auto *gh_290 = buffer.data(gh + 290);
    const auto *gh_291 = buffer.data(gh + 291);
    const auto *gh_292 = buffer.data(gh + 292);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_310 = buffer.data(gh + 310);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *hf_s_97 = buffer.data(hf_s + 97);
    const auto *hf_s_98 = buffer.data(hf_s + 98);
    const auto *hf_s_99 = buffer.data(hf_s + 99);
    const auto *hf_s_150 = buffer.data(hf_s + 150);
    const auto *hf_s_151 = buffer.data(hf_s + 151);
    const auto *hf_s_153 = buffer.data(hf_s + 153);

    const auto *hh_s_205 = buffer.data(hh_s + 205);
    const auto *hh_s_206 = buffer.data(hh_s + 206);
    const auto *hh_s_207 = buffer.data(hh_s + 207);
    const auto *hh_s_208 = buffer.data(hh_s + 208);
    const auto *hh_s_209 = buffer.data(hh_s + 209);
    const auto *hh_s_210 = buffer.data(hh_s + 210);
    const auto *hh_s_211 = buffer.data(hh_s + 211);
    const auto *hh_s_212 = buffer.data(hh_s + 212);
    const auto *hh_s_213 = buffer.data(hh_s + 213);
    const auto *hh_s_214 = buffer.data(hh_s + 214);
    const auto *hh_s_215 = buffer.data(hh_s + 215);
    const auto *hh_s_216 = buffer.data(hh_s + 216);
    const auto *hh_s_217 = buffer.data(hh_s + 217);
    const auto *hh_s_218 = buffer.data(hh_s + 218);
    const auto *hh_s_219 = buffer.data(hh_s + 219);
    const auto *hh_s_220 = buffer.data(hh_s + 220);
    const auto *hh_s_221 = buffer.data(hh_s + 221);
    const auto *hh_s_222 = buffer.data(hh_s + 222);
    const auto *hh_s_223 = buffer.data(hh_s + 223);
    const auto *hh_s_224 = buffer.data(hh_s + 224);
    const auto *hh_s_225 = buffer.data(hh_s + 225);
    const auto *hh_s_226 = buffer.data(hh_s + 226);
    const auto *hh_s_227 = buffer.data(hh_s + 227);
    const auto *hh_s_228 = buffer.data(hh_s + 228);
    const auto *hh_s_229 = buffer.data(hh_s + 229);
    const auto *hh_s_230 = buffer.data(hh_s + 230);
    const auto *hh_s_231 = buffer.data(hh_s + 231);
    const auto *hh_s_232 = buffer.data(hh_s + 232);
    const auto *hh_s_233 = buffer.data(hh_s + 233);
    const auto *hh_s_234 = buffer.data(hh_s + 234);
    const auto *hh_s_235 = buffer.data(hh_s + 235);
    const auto *hh_s_236 = buffer.data(hh_s + 236);
    const auto *hh_s_237 = buffer.data(hh_s + 237);
    const auto *hh_s_238 = buffer.data(hh_s + 238);
    const auto *hh_s_239 = buffer.data(hh_s + 239);
    const auto *hh_s_240 = buffer.data(hh_s + 240);
    const auto *hh_s_241 = buffer.data(hh_s + 241);
    const auto *hh_s_242 = buffer.data(hh_s + 242);
    const auto *hh_s_243 = buffer.data(hh_s + 243);
    const auto *hh_s_244 = buffer.data(hh_s + 244);
    const auto *hh_s_245 = buffer.data(hh_s + 245);
    const auto *hh_s_246 = buffer.data(hh_s + 246);
    const auto *hh_s_247 = buffer.data(hh_s + 247);
    const auto *hh_s_248 = buffer.data(hh_s + 248);
    const auto *hh_s_249 = buffer.data(hh_s + 249);
    const auto *hh_s_250 = buffer.data(hh_s + 250);
    const auto *hh_s_251 = buffer.data(hh_s + 251);
    const auto *hh_s_252 = buffer.data(hh_s + 252);
    const auto *hh_s_253 = buffer.data(hh_s + 253);
    const auto *hh_s_254 = buffer.data(hh_s + 254);
    const auto *hh_s_255 = buffer.data(hh_s + 255);
    const auto *hh_s_256 = buffer.data(hh_s + 256);
    const auto *hh_s_257 = buffer.data(hh_s + 257);
    const auto *hh_s_258 = buffer.data(hh_s + 258);
    const auto *hh_s_259 = buffer.data(hh_s + 259);
    const auto *hh_s_260 = buffer.data(hh_s + 260);
    const auto *hh_s_261 = buffer.data(hh_s + 261);
    const auto *hh_s_262 = buffer.data(hh_s + 262);
    const auto *hh_s_263 = buffer.data(hh_s + 263);
    const auto *hh_s_264 = buffer.data(hh_s + 264);
    const auto *hh_s_265 = buffer.data(hh_s + 265);
    const auto *hh_s_266 = buffer.data(hh_s + 266);
    const auto *hh_s_267 = buffer.data(hh_s + 267);
    const auto *hh_s_268 = buffer.data(hh_s + 268);
    const auto *hh_s_269 = buffer.data(hh_s + 269);
    const auto *hh_s_270 = buffer.data(hh_s + 270);
    const auto *hh_s_271 = buffer.data(hh_s + 271);
    const auto *hh_s_272 = buffer.data(hh_s + 272);
    const auto *hh_s_273 = buffer.data(hh_s + 273);
    const auto *hh_s_274 = buffer.data(hh_s + 274);
    const auto *hh_s_275 = buffer.data(hh_s + 275);
    const auto *hh_s_276 = buffer.data(hh_s + 276);
    const auto *hh_s_277 = buffer.data(hh_s + 277);
    const auto *hh_s_278 = buffer.data(hh_s + 278);
    const auto *hh_s_279 = buffer.data(hh_s + 279);
    const auto *hh_s_280 = buffer.data(hh_s + 280);
    const auto *hh_s_281 = buffer.data(hh_s + 281);
    const auto *hh_s_282 = buffer.data(hh_s + 282);
    const auto *hh_s_283 = buffer.data(hh_s + 283);
    const auto *hh_s_284 = buffer.data(hh_s + 284);
    const auto *hh_s_285 = buffer.data(hh_s + 285);
    const auto *hh_s_286 = buffer.data(hh_s + 286);
    const auto *hh_s_287 = buffer.data(hh_s + 287);
    const auto *hh_s_288 = buffer.data(hh_s + 288);
    const auto *hh_s_289 = buffer.data(hh_s + 289);
    const auto *hh_s_290 = buffer.data(hh_s + 290);
    const auto *hh_s_291 = buffer.data(hh_s + 291);
    const auto *hh_s_292 = buffer.data(hh_s + 292);
    const auto *hh_s_293 = buffer.data(hh_s + 293);
    const auto *hh_s_294 = buffer.data(hh_s + 294);
    const auto *hh_s_295 = buffer.data(hh_s + 295);
    const auto *hh_s_296 = buffer.data(hh_s + 296);
    const auto *hh_s_297 = buffer.data(hh_s + 297);
    const auto *hh_s_298 = buffer.data(hh_s + 298);
    const auto *hh_s_299 = buffer.data(hh_s + 299);
    const auto *hh_s_300 = buffer.data(hh_s + 300);
    const auto *hh_s_301 = buffer.data(hh_s + 301);
    const auto *hh_s_302 = buffer.data(hh_s + 302);
    const auto *hh_s_303 = buffer.data(hh_s + 303);
    const auto *hh_s_304 = buffer.data(hh_s + 304);
    const auto *hh_s_305 = buffer.data(hh_s + 305);
    const auto *hh_s_306 = buffer.data(hh_s + 306);
    const auto *hh_s_307 = buffer.data(hh_s + 307);
    const auto *hh_s_308 = buffer.data(hh_s + 308);
    const auto *hh_s_309 = buffer.data(hh_s + 309);
    const auto *hh_s_310 = buffer.data(hh_s + 310);
    const auto *hh_s_311 = buffer.data(hh_s + 311);
    const auto *hh_s_312 = buffer.data(hh_s + 312);
    const auto *hh_s_313 = buffer.data(hh_s + 313);
    const auto *hh_s_314 = buffer.data(hh_s + 314);
    const auto *hh_s_315 = buffer.data(hh_s + 315);
    const auto *hh_s_316 = buffer.data(hh_s + 316);
    const auto *hh_s_317 = buffer.data(hh_s + 317);
    const auto *hh_s_318 = buffer.data(hh_s + 318);

    const auto *hf_97 = buffer.data(hf + 97);
    const auto *hf_98 = buffer.data(hf + 98);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_153 = buffer.data(hf + 153);

    const auto *hg_146 = buffer.data(hg + 146);
    const auto *hg_147 = buffer.data(hg + 147);
    const auto *hg_148 = buffer.data(hg + 148);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_155 = buffer.data(hg + 155);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_162 = buffer.data(hg + 162);
    const auto *hg_163 = buffer.data(hg + 163);
    const auto *hg_164 = buffer.data(hg + 164);
    const auto *hg_165 = buffer.data(hg + 165);
    const auto *hg_167 = buffer.data(hg + 167);
    const auto *hg_168 = buffer.data(hg + 168);
    const auto *hg_170 = buffer.data(hg + 170);
    const auto *hg_176 = buffer.data(hg + 176);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_178 = buffer.data(hg + 178);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_180 = buffer.data(hg + 180);
    const auto *hg_182 = buffer.data(hg + 182);
    const auto *hg_183 = buffer.data(hg + 183);
    const auto *hg_185 = buffer.data(hg + 185);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_191 = buffer.data(hg + 191);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_193 = buffer.data(hg + 193);
    const auto *hg_194 = buffer.data(hg + 194);
    const auto *hg_195 = buffer.data(hg + 195);
    const auto *hg_197 = buffer.data(hg + 197);
    const auto *hg_198 = buffer.data(hg + 198);
    const auto *hg_200 = buffer.data(hg + 200);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_206 = buffer.data(hg + 206);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_208 = buffer.data(hg + 208);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_220 = buffer.data(hg + 220);
    const auto *hg_221 = buffer.data(hg + 221);
    const auto *hg_222 = buffer.data(hg + 222);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);

#pragma omp simd aligned(t_205, t_206, t_207, pb_y, hf_s_97, hf_s_98, hf_s_99, hh_s_205, \
                         hh_s_206, hh_s_207, hf_97, hf_98, hf_99, hg_146, hg_147, \
                         hg_148 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_205[k] = -f_10 * hf_s_97[k]
                   + f_2 * hh_s_205[k]
                   + f_8 * hf_97[k]
                   + pb_y[k] * hg_146[k];

        t_206[k] = -f_6 * hf_s_98[k]
                   + f_2 * hh_s_206[k]
                   + f_7 * hf_98[k]
                   + pb_y[k] * hg_147[k];

        t_207[k] = -f_4 * hf_s_99[k]
                   + f_2 * hh_s_207[k]
                   + f_5 * hf_99[k]
                   + pb_y[k] * hg_148[k];
    }

#pragma omp simd aligned(t_208, t_209, t_210, pa_x, pb_y, fh_s_209, fh_209, gg_150, gh_209, \
                         gh_210, hh_s_208, hh_s_209, hh_s_210, hg_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_208[k] = f_2 * hh_s_208[k]
                   + pb_y[k] * hg_149[k];

        t_209[k] = -f_11 * fh_s_209[k]
                   + f_5 * fh_209[k]
                   + pa_x[k] * gh_209[k]
                   + f_2 * hh_s_209[k];

        t_210[k] = f_0 * gg_150[k]
                   + pa_x[k] * gh_210[k]
                   + f_2 * hh_s_210[k];
    }

#pragma omp simd aligned(t_211, t_212, t_213, t_214, pa_x, pb_y, pb_z, gg_90, gg_153, gh_213, \
                         hh_s_211, hh_s_212, hh_s_213, hh_s_214, hg_150, \
                         hg_151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_211[k] = f_3 * gg_90[k]
                   + f_2 * hh_s_211[k]
                   + pb_y[k] * hg_150[k];

        t_212[k] = f_2 * hh_s_212[k]
                   + pb_z[k] * hg_150[k];

        t_213[k] = f_8 * gg_153[k]
                   + pa_x[k] * gh_213[k]
                   + f_2 * hh_s_213[k];

        t_214[k] = f_2 * hh_s_214[k]
                   + pb_z[k] * hg_151[k];
    }

#pragma omp simd aligned(t_215, t_216, t_217, pa_x, pb_z, gg_155, gg_156, gh_215, gh_216, \
                         hh_s_215, hh_s_216, hh_s_217, hg_153 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_215[k] = f_8 * gg_155[k]
                   + pa_x[k] * gh_215[k]
                   + f_2 * hh_s_215[k];

        t_216[k] = f_7 * gg_156[k]
                   + pa_x[k] * gh_216[k]
                   + f_2 * hh_s_216[k];

        t_217[k] = f_2 * hh_s_217[k]
                   + pb_z[k] * hg_153[k];
    }

#pragma omp simd aligned(t_218, t_219, t_220, pa_x, pb_x, pb_y, gg_95, gg_159, gg_160, gh_219, \
                         hh_s_218, hh_s_219, hh_s_220, hg_155, hg_160 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_218[k] = f_3 * gg_95[k]
                   + f_2 * hh_s_218[k]
                   + pb_y[k] * hg_155[k];

        t_219[k] = f_7 * gg_159[k]
                   + pa_x[k] * gh_219[k]
                   + f_2 * hh_s_219[k];

        t_220[k] = f_5 * gg_160[k]
                   + f_2 * hh_s_220[k]
                   + pb_x[k] * hg_160[k];
    }

#pragma omp simd aligned(t_221, t_222, t_223, pb_x, pb_z, gg_162, gg_163, hh_s_221, hh_s_222, \
                         hh_s_223, hg_156, hg_162, hg_163 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_221[k] = f_2 * hh_s_221[k]
                   + pb_z[k] * hg_156[k];

        t_222[k] = f_5 * gg_162[k]
                   + f_2 * hh_s_222[k]
                   + pb_x[k] * hg_162[k];

        t_223[k] = f_5 * gg_163[k]
                   + f_2 * hh_s_223[k]
                   + pb_x[k] * hg_163[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_x, pb_x, pb_z, gg_164, gh_225, gh_227, \
                         hh_s_224, hh_s_225, hh_s_226, hh_s_227, hg_160, \
                         hg_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_5 * gg_164[k]
                   + f_2 * hh_s_224[k]
                   + pb_x[k] * hg_164[k];

        t_225[k] = pa_x[k] * gh_225[k]
                   + f_2 * hh_s_225[k];

        t_226[k] = f_2 * hh_s_226[k]
                   + pb_z[k] * hg_160[k];

        t_227[k] = pa_x[k] * gh_227[k]
                   + f_2 * hh_s_227[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, t_231, pa_x, pa_z, gh_126, gh_228, gh_229, \
                         gh_230, hh_s_228, hh_s_229, hh_s_230, \
                         hh_s_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = pa_x[k] * gh_228[k]
                   + f_2 * hh_s_228[k];

        t_229[k] = pa_x[k] * gh_229[k]
                   + f_2 * hh_s_229[k];

        t_230[k] = pa_x[k] * gh_230[k]
                   + f_2 * hh_s_230[k];

        t_231[k] = pa_z[k] * gh_126[k]
                   + f_2 * hh_s_231[k];
    }

#pragma omp simd aligned(t_232, t_233, t_234, pa_z, pb_z, gg_90, gh_127, gh_129, hh_s_232, \
                         hh_s_233, hh_s_234, hg_165 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_232[k] = pa_z[k] * gh_127[k]
                   + f_2 * hh_s_232[k];

        t_233[k] = f_5 * gg_90[k]
                   + f_2 * hh_s_233[k]
                   + pb_z[k] * hg_165[k];

        t_234[k] = pa_z[k] * gh_129[k]
                   + f_2 * hh_s_234[k];
    }

#pragma omp simd aligned(t_235, t_236, t_237, pa_x, pa_z, pb_y, gg_107, gg_170, gh_132, \
                         gh_236, hh_s_235, hh_s_236, hh_s_237, hg_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_235[k] = f_8 * gg_107[k]
                   + f_2 * hh_s_235[k]
                   + pb_y[k] * hg_167[k];

        t_236[k] = f_8 * gg_170[k]
                   + pa_x[k] * gh_236[k]
                   + f_2 * hh_s_236[k];

        t_237[k] = pa_z[k] * gh_132[k]
                   + f_2 * hh_s_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, pa_x, pb_y, pb_z, gg_93, gg_110, gg_174, gh_240, \
                         hh_s_238, hh_s_239, hh_s_240, hg_168, hg_170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = f_5 * gg_93[k]
                   + f_2 * hh_s_238[k]
                   + pb_z[k] * hg_168[k];

        t_239[k] = f_8 * gg_110[k]
                   + f_2 * hh_s_239[k]
                   + pb_y[k] * hg_170[k];

        t_240[k] = f_7 * gg_174[k]
                   + pa_x[k] * gh_240[k]
                   + f_2 * hh_s_240[k];
    }

#pragma omp simd aligned(t_241, t_242, t_243, pa_z, pb_x, gg_176, gg_177, gh_136, hh_s_241, \
                         hh_s_242, hh_s_243, hg_176, hg_177 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_241[k] = pa_z[k] * gh_136[k]
                   + f_2 * hh_s_241[k];

        t_242[k] = f_5 * gg_176[k]
                   + f_2 * hh_s_242[k]
                   + pb_x[k] * hg_176[k];

        t_243[k] = f_5 * gg_177[k]
                   + f_2 * hh_s_243[k]
                   + pb_x[k] * hg_177[k];
    }

#pragma omp simd aligned(t_244, t_245, t_246, t_247, pa_x, pb_x, gg_178, gg_179, gh_246, \
                         gh_247, hh_s_244, hh_s_245, hh_s_246, hh_s_247, hg_178, \
                         hg_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_244[k] = f_5 * gg_178[k]
                   + f_2 * hh_s_244[k]
                   + pb_x[k] * hg_178[k];

        t_245[k] = f_5 * gg_179[k]
                   + f_2 * hh_s_245[k]
                   + pb_x[k] * hg_179[k];

        t_246[k] = pa_x[k] * gh_246[k]
                   + f_2 * hh_s_246[k];

        t_247[k] = pa_x[k] * gh_247[k]
                   + f_2 * hh_s_247[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_x, gh_248, gh_249, gh_250, gh_251, \
                         hh_s_248, hh_s_249, hh_s_250, hh_s_251 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = pa_x[k] * gh_248[k]
                   + f_2 * hh_s_248[k];

        t_249[k] = pa_x[k] * gh_249[k]
                   + f_2 * hh_s_249[k];

        t_250[k] = pa_x[k] * gh_250[k]
                   + f_2 * hh_s_250[k];

        t_251[k] = pa_x[k] * gh_251[k]
                   + f_2 * hh_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, pa_x, pb_y, pb_z, gg_105, gg_120, gg_180, \
                         gh_252, hh_s_252, hh_s_253, hh_s_254, hg_180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = f_0 * gg_180[k]
                   + pa_x[k] * gh_252[k]
                   + f_2 * hh_s_252[k];

        t_253[k] = f_7 * gg_120[k]
                   + f_2 * hh_s_253[k]
                   + pb_y[k] * hg_180[k];

        t_254[k] = f_7 * gg_105[k]
                   + f_2 * hh_s_254[k]
                   + pb_z[k] * hg_180[k];
    }

#pragma omp simd aligned(t_255, t_256, t_257, pa_x, pb_y, gg_122, gg_183, gg_185, gh_255, \
                         gh_257, hh_s_255, hh_s_256, hh_s_257, hg_182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_255[k] = f_8 * gg_183[k]
                   + pa_x[k] * gh_255[k]
                   + f_2 * hh_s_255[k];

        t_256[k] = f_7 * gg_122[k]
                   + f_2 * hh_s_256[k]
                   + pb_y[k] * hg_182[k];

        t_257[k] = f_8 * gg_185[k]
                   + pa_x[k] * gh_257[k]
                   + f_2 * hh_s_257[k];
    }

#pragma omp simd aligned(t_258, t_259, t_260, pa_x, pb_y, pb_z, gg_108, gg_125, gg_186, \
                         gh_258, hh_s_258, hh_s_259, hh_s_260, hg_183, \
                         hg_185 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_258[k] = f_7 * gg_186[k]
                   + pa_x[k] * gh_258[k]
                   + f_2 * hh_s_258[k];

        t_259[k] = f_7 * gg_108[k]
                   + f_2 * hh_s_259[k]
                   + pb_z[k] * hg_183[k];

        t_260[k] = f_7 * gg_125[k]
                   + f_2 * hh_s_260[k]
                   + pb_y[k] * hg_185[k];
    }

#pragma omp simd aligned(t_261, t_262, t_263, pa_x, pb_x, gg_189, gg_190, gg_191, gh_261, \
                         hh_s_261, hh_s_262, hh_s_263, hg_190, hg_191 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_261[k] = f_7 * gg_189[k]
                   + pa_x[k] * gh_261[k]
                   + f_2 * hh_s_261[k];

        t_262[k] = f_5 * gg_190[k]
                   + f_2 * hh_s_262[k]
                   + pb_x[k] * hg_190[k];

        t_263[k] = f_5 * gg_191[k]
                   + f_2 * hh_s_263[k]
                   + pb_x[k] * hg_191[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pb_x, gg_192, gg_193, gg_194, hh_s_264, \
                         hh_s_265, hh_s_266, hg_192, hg_193, hg_194 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_5 * gg_192[k]
                   + f_2 * hh_s_264[k]
                   + pb_x[k] * hg_192[k];

        t_265[k] = f_5 * gg_193[k]
                   + f_2 * hh_s_265[k]
                   + pb_x[k] * hg_193[k];

        t_266[k] = f_5 * gg_194[k]
                   + f_2 * hh_s_266[k]
                   + pb_x[k] * hg_194[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, t_270, t_271, pa_x, gh_267, gh_268, gh_269, \
                         gh_270, gh_271, hh_s_267, hh_s_268, hh_s_269, hh_s_270, \
                         hh_s_271 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = pa_x[k] * gh_267[k]
                   + f_2 * hh_s_267[k];

        t_268[k] = pa_x[k] * gh_268[k]
                   + f_2 * hh_s_268[k];

        t_269[k] = pa_x[k] * gh_269[k]
                   + f_2 * hh_s_269[k];

        t_270[k] = pa_x[k] * gh_270[k]
                   + f_2 * hh_s_270[k];

        t_271[k] = pa_x[k] * gh_271[k]
                   + f_2 * hh_s_271[k];
    }

#pragma omp simd aligned(t_272, t_273, t_274, t_275, pa_x, pa_y, pb_y, gg_135, gh_189, gh_191, \
                         gh_272, hh_s_272, hh_s_273, hh_s_274, hh_s_275, \
                         hg_195 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_272[k] = pa_x[k] * gh_272[k]
                   + f_2 * hh_s_272[k];

        t_273[k] = pa_y[k] * gh_189[k]
                   + f_2 * hh_s_273[k];

        t_274[k] = f_5 * gg_135[k]
                   + f_2 * hh_s_274[k]
                   + pb_y[k] * hg_195[k];

        t_275[k] = pa_y[k] * gh_191[k]
                   + f_2 * hh_s_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pa_x, pa_y, pb_y, gg_137, gg_198, gh_194, \
                         gh_276, hh_s_276, hh_s_277, hh_s_278, hg_197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_8 * gg_198[k]
                   + pa_x[k] * gh_276[k]
                   + f_2 * hh_s_276[k];

        t_277[k] = f_5 * gg_137[k]
                   + f_2 * hh_s_277[k]
                   + pb_y[k] * hg_197[k];

        t_278[k] = pa_y[k] * gh_194[k]
                   + f_2 * hh_s_278[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, pa_x, pb_y, pb_z, gg_123, gg_140, gg_201, \
                         gh_279, hh_s_279, hh_s_280, hh_s_281, hg_198, \
                         hg_200 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_7 * gg_201[k]
                   + pa_x[k] * gh_279[k]
                   + f_2 * hh_s_279[k];

        t_280[k] = f_8 * gg_123[k]
                   + f_2 * hh_s_280[k]
                   + pb_z[k] * hg_198[k];

        t_281[k] = f_5 * gg_140[k]
                   + f_2 * hh_s_281[k]
                   + pb_y[k] * hg_200[k];
    }

#pragma omp simd aligned(t_282, t_283, t_284, pa_y, pb_x, gg_205, gg_206, gh_198, hh_s_282, \
                         hh_s_283, hh_s_284, hg_205, hg_206 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_282[k] = pa_y[k] * gh_198[k]
                   + f_2 * hh_s_282[k];

        t_283[k] = f_5 * gg_205[k]
                   + f_2 * hh_s_283[k]
                   + pb_x[k] * hg_205[k];

        t_284[k] = f_5 * gg_206[k]
                   + f_2 * hh_s_284[k]
                   + pb_x[k] * hg_206[k];
    }

#pragma omp simd aligned(t_285, t_286, t_287, pa_y, pb_x, gg_207, gg_208, gh_203, hh_s_285, \
                         hh_s_286, hh_s_287, hg_207, hg_208 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_285[k] = f_5 * gg_207[k]
                   + f_2 * hh_s_285[k]
                   + pb_x[k] * hg_207[k];

        t_286[k] = f_5 * gg_208[k]
                   + f_2 * hh_s_286[k]
                   + pb_x[k] * hg_208[k];

        t_287[k] = pa_y[k] * gh_203[k]
                   + f_2 * hh_s_287[k];
    }

#pragma omp simd aligned(t_288, t_289, t_290, t_291, t_292, pa_x, gh_288, gh_289, gh_290, \
                         gh_291, gh_292, hh_s_288, hh_s_289, hh_s_290, hh_s_291, \
                         hh_s_292 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_288[k] = pa_x[k] * gh_288[k]
                   + f_2 * hh_s_288[k];

        t_289[k] = pa_x[k] * gh_289[k]
                   + f_2 * hh_s_289[k];

        t_290[k] = pa_x[k] * gh_290[k]
                   + f_2 * hh_s_290[k];

        t_291[k] = pa_x[k] * gh_291[k]
                   + f_2 * hh_s_291[k];

        t_292[k] = pa_x[k] * gh_292[k]
                   + f_2 * hh_s_292[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_x, pb_y, pb_z, gg_135, gg_210, gh_293, \
                         gh_294, hh_s_293, hh_s_294, hh_s_295, hh_s_296, \
                         hg_210 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = pa_x[k] * gh_293[k]
                   + f_2 * hh_s_293[k];

        t_294[k] = f_0 * gg_210[k]
                   + pa_x[k] * gh_294[k]
                   + f_2 * hh_s_294[k];

        t_295[k] = f_2 * hh_s_295[k]
                   + pb_y[k] * hg_210[k];

        t_296[k] = f_3 * gg_135[k]
                   + f_2 * hh_s_296[k]
                   + pb_z[k] * hg_210[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, pa_x, pb_y, gg_213, gg_215, gh_297, gh_299, \
                         hh_s_297, hh_s_298, hh_s_299, hg_212 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = f_8 * gg_213[k]
                   + pa_x[k] * gh_297[k]
                   + f_2 * hh_s_297[k];

        t_298[k] = f_2 * hh_s_298[k]
                   + pb_y[k] * hg_212[k];

        t_299[k] = f_8 * gg_215[k]
                   + pa_x[k] * gh_299[k]
                   + f_2 * hh_s_299[k];
    }

#pragma omp simd aligned(t_300, t_301, t_302, pa_x, pb_y, gg_216, gg_217, gh_300, gh_301, \
                         hh_s_300, hh_s_301, hh_s_302, hg_215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_300[k] = f_7 * gg_216[k]
                   + pa_x[k] * gh_300[k]
                   + f_2 * hh_s_300[k];

        t_301[k] = f_7 * gg_217[k]
                   + pa_x[k] * gh_301[k]
                   + f_2 * hh_s_301[k];

        t_302[k] = f_2 * hh_s_302[k]
                   + pb_y[k] * hg_215[k];
    }

#pragma omp simd aligned(t_303, t_304, t_305, pa_x, pb_x, gg_219, gg_220, gg_221, gh_303, \
                         hh_s_303, hh_s_304, hh_s_305, hg_220, hg_221 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_303[k] = f_7 * gg_219[k]
                   + pa_x[k] * gh_303[k]
                   + f_2 * hh_s_303[k];

        t_304[k] = f_5 * gg_220[k]
                   + f_2 * hh_s_304[k]
                   + pb_x[k] * hg_220[k];

        t_305[k] = f_5 * gg_221[k]
                   + f_2 * hh_s_305[k]
                   + pb_x[k] * hg_221[k];
    }

#pragma omp simd aligned(t_306, t_307, t_308, pb_x, pb_y, gg_222, gg_224, hh_s_306, hh_s_307, \
                         hh_s_308, hg_219, hg_222, hg_224 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_306[k] = f_5 * gg_222[k]
                   + f_2 * hh_s_306[k]
                   + pb_x[k] * hg_222[k];

        t_307[k] = f_2 * hh_s_307[k]
                   + pb_y[k] * hg_219[k];

        t_308[k] = f_5 * gg_224[k]
                   + f_2 * hh_s_308[k]
                   + pb_x[k] * hg_224[k];
    }

#pragma omp simd aligned(t_309, t_310, t_311, t_312, pa_x, gh_309, gh_310, gh_311, gh_312, \
                         hh_s_309, hh_s_310, hh_s_311, hh_s_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_309[k] = pa_x[k] * gh_309[k]
                   + f_2 * hh_s_309[k];

        t_310[k] = pa_x[k] * gh_310[k]
                   + f_2 * hh_s_310[k];

        t_311[k] = pa_x[k] * gh_311[k]
                   + f_2 * hh_s_311[k];

        t_312[k] = pa_x[k] * gh_312[k]
                   + f_2 * hh_s_312[k];
    }

#pragma omp simd aligned(t_313, t_314, t_315, pa_x, pb_x, pb_y, gh_314, hf_s_150, hh_s_313, \
                         hh_s_314, hh_s_315, hf_150, hg_224, hg_225 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_313[k] = f_2 * hh_s_313[k]
                   + pb_y[k] * hg_224[k];

        t_314[k] = pa_x[k] * gh_314[k]
                   + f_2 * hh_s_314[k];

        t_315[k] = -f_1 * hf_s_150[k]
                   + f_2 * hh_s_315[k]
                   + f_3 * hf_150[k]
                   + pb_x[k] * hg_225[k];
    }

#pragma omp simd aligned(t_316, t_317, t_318, pb_x, pb_z, hf_s_151, hf_s_153, hh_s_316, \
                         hh_s_317, hh_s_318, hf_151, hf_153, hg_225, hg_226, \
                         hg_228 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_316[k] = -f_10 * hf_s_151[k]
                   + f_2 * hh_s_316[k]
                   + f_8 * hf_151[k]
                   + pb_x[k] * hg_226[k];

        t_317[k] = f_2 * hh_s_317[k]
                   + pb_z[k] * hg_225[k];

        t_318[k] = -f_6 * hf_s_153[k]
                   + f_2 * hh_s_318[k]
                   + f_7 * hf_153[k]
                   + pb_x[k] * hg_228[k];
    }
}

static auto
compute_prim_hh_kinetic_energy_0_piece3(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t fh_s,
                                        const size_t fh, const size_t gg, const size_t gh,
                                        const size_t hf_s, const size_t hh_s, const size_t hf,
                                        const size_t hg, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_9 = 3.0 * beta / p;
    const auto f_10 = 3.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 2.0 * beta / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *fh_s_141 = buffer.data(fh_s + 141);
    const auto *fh_s_162 = buffer.data(fh_s + 162);
    const auto *fh_s_167 = buffer.data(fh_s + 167);
    const auto *fh_s_188 = buffer.data(fh_s + 188);
    const auto *fh_s_209 = buffer.data(fh_s + 209);

    const auto *fh_141 = buffer.data(fh + 141);
    const auto *fh_162 = buffer.data(fh + 162);
    const auto *fh_167 = buffer.data(fh + 167);
    const auto *fh_188 = buffer.data(fh + 188);
    const auto *fh_209 = buffer.data(fh + 209);

    const auto *gg_151 = buffer.data(gg + 151);
    const auto *gg_153 = buffer.data(gg + 153);
    const auto *gg_154 = buffer.data(gg + 154);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_161 = buffer.data(gg + 161);
    const auto *gg_162 = buffer.data(gg + 162);
    const auto *gg_164 = buffer.data(gg + 164);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_193 = buffer.data(gg + 193);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_208 = buffer.data(gg + 208);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_210 = buffer.data(gg + 210);
    const auto *gg_211 = buffer.data(gg + 211);
    const auto *gg_212 = buffer.data(gg + 212);
    const auto *gg_213 = buffer.data(gg + 213);
    const auto *gg_214 = buffer.data(gg + 214);
    const auto *gg_215 = buffer.data(gg + 215);
    const auto *gg_220 = buffer.data(gg + 220);
    const auto *gg_222 = buffer.data(gg + 222);
    const auto *gg_223 = buffer.data(gg + 223);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *gh_210 = buffer.data(gh + 210);
    const auto *gh_211 = buffer.data(gh + 211);
    const auto *gh_213 = buffer.data(gh + 213);
    const auto *gh_214 = buffer.data(gh + 214);
    const auto *gh_216 = buffer.data(gh + 216);
    const auto *gh_217 = buffer.data(gh + 217);
    const auto *gh_218 = buffer.data(gh + 218);
    const auto *gh_225 = buffer.data(gh + 225);
    const auto *gh_227 = buffer.data(gh + 227);
    const auto *gh_228 = buffer.data(gh + 228);
    const auto *gh_246 = buffer.data(gh + 246);
    const auto *gh_251 = buffer.data(gh + 251);
    const auto *gh_267 = buffer.data(gh + 267);
    const auto *gh_272 = buffer.data(gh + 272);
    const auto *gh_293 = buffer.data(gh + 293);
    const auto *gh_294 = buffer.data(gh + 294);
    const auto *gh_295 = buffer.data(gh + 295);
    const auto *gh_296 = buffer.data(gh + 296);
    const auto *gh_297 = buffer.data(gh + 297);
    const auto *gh_298 = buffer.data(gh + 298);
    const auto *gh_299 = buffer.data(gh + 299);
    const auto *gh_300 = buffer.data(gh + 300);
    const auto *gh_301 = buffer.data(gh + 301);
    const auto *gh_302 = buffer.data(gh + 302);
    const auto *gh_303 = buffer.data(gh + 303);
    const auto *gh_309 = buffer.data(gh + 309);
    const auto *gh_311 = buffer.data(gh + 311);
    const auto *gh_312 = buffer.data(gh + 312);
    const auto *gh_314 = buffer.data(gh + 314);

    const auto *hf_s_155 = buffer.data(hf_s + 155);
    const auto *hf_s_156 = buffer.data(hf_s + 156);
    const auto *hf_s_157 = buffer.data(hf_s + 157);
    const auto *hf_s_158 = buffer.data(hf_s + 158);
    const auto *hf_s_159 = buffer.data(hf_s + 159);
    const auto *hf_s_162 = buffer.data(hf_s + 162);
    const auto *hf_s_165 = buffer.data(hf_s + 165);
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
    const auto *hf_s_200 = buffer.data(hf_s + 200);

    const auto *hh_s_319 = buffer.data(hh_s + 319);
    const auto *hh_s_320 = buffer.data(hh_s + 320);
    const auto *hh_s_321 = buffer.data(hh_s + 321);
    const auto *hh_s_322 = buffer.data(hh_s + 322);
    const auto *hh_s_323 = buffer.data(hh_s + 323);
    const auto *hh_s_324 = buffer.data(hh_s + 324);
    const auto *hh_s_325 = buffer.data(hh_s + 325);
    const auto *hh_s_326 = buffer.data(hh_s + 326);
    const auto *hh_s_327 = buffer.data(hh_s + 327);
    const auto *hh_s_328 = buffer.data(hh_s + 328);
    const auto *hh_s_329 = buffer.data(hh_s + 329);
    const auto *hh_s_330 = buffer.data(hh_s + 330);
    const auto *hh_s_331 = buffer.data(hh_s + 331);
    const auto *hh_s_332 = buffer.data(hh_s + 332);
    const auto *hh_s_333 = buffer.data(hh_s + 333);
    const auto *hh_s_334 = buffer.data(hh_s + 334);
    const auto *hh_s_335 = buffer.data(hh_s + 335);
    const auto *hh_s_336 = buffer.data(hh_s + 336);
    const auto *hh_s_337 = buffer.data(hh_s + 337);
    const auto *hh_s_338 = buffer.data(hh_s + 338);
    const auto *hh_s_339 = buffer.data(hh_s + 339);
    const auto *hh_s_340 = buffer.data(hh_s + 340);
    const auto *hh_s_341 = buffer.data(hh_s + 341);
    const auto *hh_s_342 = buffer.data(hh_s + 342);
    const auto *hh_s_343 = buffer.data(hh_s + 343);
    const auto *hh_s_344 = buffer.data(hh_s + 344);
    const auto *hh_s_345 = buffer.data(hh_s + 345);
    const auto *hh_s_346 = buffer.data(hh_s + 346);
    const auto *hh_s_347 = buffer.data(hh_s + 347);
    const auto *hh_s_348 = buffer.data(hh_s + 348);
    const auto *hh_s_349 = buffer.data(hh_s + 349);
    const auto *hh_s_350 = buffer.data(hh_s + 350);
    const auto *hh_s_351 = buffer.data(hh_s + 351);
    const auto *hh_s_352 = buffer.data(hh_s + 352);
    const auto *hh_s_353 = buffer.data(hh_s + 353);
    const auto *hh_s_354 = buffer.data(hh_s + 354);
    const auto *hh_s_355 = buffer.data(hh_s + 355);
    const auto *hh_s_356 = buffer.data(hh_s + 356);
    const auto *hh_s_357 = buffer.data(hh_s + 357);
    const auto *hh_s_358 = buffer.data(hh_s + 358);
    const auto *hh_s_359 = buffer.data(hh_s + 359);
    const auto *hh_s_360 = buffer.data(hh_s + 360);
    const auto *hh_s_361 = buffer.data(hh_s + 361);
    const auto *hh_s_362 = buffer.data(hh_s + 362);
    const auto *hh_s_363 = buffer.data(hh_s + 363);
    const auto *hh_s_364 = buffer.data(hh_s + 364);
    const auto *hh_s_365 = buffer.data(hh_s + 365);
    const auto *hh_s_366 = buffer.data(hh_s + 366);
    const auto *hh_s_367 = buffer.data(hh_s + 367);
    const auto *hh_s_368 = buffer.data(hh_s + 368);
    const auto *hh_s_369 = buffer.data(hh_s + 369);
    const auto *hh_s_370 = buffer.data(hh_s + 370);
    const auto *hh_s_371 = buffer.data(hh_s + 371);
    const auto *hh_s_372 = buffer.data(hh_s + 372);
    const auto *hh_s_373 = buffer.data(hh_s + 373);
    const auto *hh_s_374 = buffer.data(hh_s + 374);
    const auto *hh_s_375 = buffer.data(hh_s + 375);
    const auto *hh_s_376 = buffer.data(hh_s + 376);
    const auto *hh_s_377 = buffer.data(hh_s + 377);
    const auto *hh_s_378 = buffer.data(hh_s + 378);
    const auto *hh_s_379 = buffer.data(hh_s + 379);
    const auto *hh_s_380 = buffer.data(hh_s + 380);
    const auto *hh_s_381 = buffer.data(hh_s + 381);
    const auto *hh_s_382 = buffer.data(hh_s + 382);
    const auto *hh_s_383 = buffer.data(hh_s + 383);
    const auto *hh_s_384 = buffer.data(hh_s + 384);
    const auto *hh_s_385 = buffer.data(hh_s + 385);
    const auto *hh_s_386 = buffer.data(hh_s + 386);
    const auto *hh_s_387 = buffer.data(hh_s + 387);
    const auto *hh_s_388 = buffer.data(hh_s + 388);
    const auto *hh_s_389 = buffer.data(hh_s + 389);
    const auto *hh_s_390 = buffer.data(hh_s + 390);
    const auto *hh_s_391 = buffer.data(hh_s + 391);
    const auto *hh_s_392 = buffer.data(hh_s + 392);
    const auto *hh_s_393 = buffer.data(hh_s + 393);
    const auto *hh_s_394 = buffer.data(hh_s + 394);
    const auto *hh_s_395 = buffer.data(hh_s + 395);
    const auto *hh_s_396 = buffer.data(hh_s + 396);
    const auto *hh_s_397 = buffer.data(hh_s + 397);
    const auto *hh_s_398 = buffer.data(hh_s + 398);
    const auto *hh_s_399 = buffer.data(hh_s + 399);
    const auto *hh_s_400 = buffer.data(hh_s + 400);
    const auto *hh_s_401 = buffer.data(hh_s + 401);
    const auto *hh_s_402 = buffer.data(hh_s + 402);
    const auto *hh_s_403 = buffer.data(hh_s + 403);
    const auto *hh_s_404 = buffer.data(hh_s + 404);
    const auto *hh_s_405 = buffer.data(hh_s + 405);
    const auto *hh_s_406 = buffer.data(hh_s + 406);
    const auto *hh_s_407 = buffer.data(hh_s + 407);
    const auto *hh_s_408 = buffer.data(hh_s + 408);
    const auto *hh_s_409 = buffer.data(hh_s + 409);
    const auto *hh_s_410 = buffer.data(hh_s + 410);
    const auto *hh_s_411 = buffer.data(hh_s + 411);
    const auto *hh_s_412 = buffer.data(hh_s + 412);
    const auto *hh_s_413 = buffer.data(hh_s + 413);
    const auto *hh_s_414 = buffer.data(hh_s + 414);
    const auto *hh_s_415 = buffer.data(hh_s + 415);
    const auto *hh_s_416 = buffer.data(hh_s + 416);
    const auto *hh_s_417 = buffer.data(hh_s + 417);
    const auto *hh_s_418 = buffer.data(hh_s + 418);
    const auto *hh_s_419 = buffer.data(hh_s + 419);
    const auto *hh_s_420 = buffer.data(hh_s + 420);

    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_162 = buffer.data(hf + 162);
    const auto *hf_165 = buffer.data(hf + 165);
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
    const auto *hf_200 = buffer.data(hf + 200);

    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_231 = buffer.data(hg + 231);
    const auto *hg_233 = buffer.data(hg + 233);
    const auto *hg_234 = buffer.data(hg + 234);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_236 = buffer.data(hg + 236);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_242 = buffer.data(hg + 242);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_249 = buffer.data(hg + 249);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_256 = buffer.data(hg + 256);
    const auto *hg_257 = buffer.data(hg + 257);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_259 = buffer.data(hg + 259);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_261 = buffer.data(hg + 261);
    const auto *hg_262 = buffer.data(hg + 262);
    const auto *hg_263 = buffer.data(hg + 263);
    const auto *hg_264 = buffer.data(hg + 264);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_271 = buffer.data(hg + 271);
    const auto *hg_272 = buffer.data(hg + 272);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_274 = buffer.data(hg + 274);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_276 = buffer.data(hg + 276);
    const auto *hg_277 = buffer.data(hg + 277);
    const auto *hg_278 = buffer.data(hg + 278);
    const auto *hg_279 = buffer.data(hg + 279);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);

#pragma omp simd aligned(t_319, t_320, t_321, pb_x, pb_z, hf_s_155, hf_s_156, hh_s_319, \
                         hh_s_320, hh_s_321, hf_155, hf_156, hg_226, hg_230, \
                         hg_231 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_319[k] = f_2 * hh_s_319[k]
                   + pb_z[k] * hg_226[k];

        t_320[k] = -f_6 * hf_s_155[k]
                   + f_2 * hh_s_320[k]
                   + f_7 * hf_155[k]
                   + pb_x[k] * hg_230[k];

        t_321[k] = -f_4 * hf_s_156[k]
                   + f_2 * hh_s_321[k]
                   + f_5 * hf_156[k]
                   + pb_x[k] * hg_231[k];
    }

#pragma omp simd aligned(t_322, t_323, t_324, pb_x, pb_z, hf_s_158, hf_s_159, hh_s_322, \
                         hh_s_323, hh_s_324, hf_158, hf_159, hg_228, hg_233, \
                         hg_234 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_322[k] = f_2 * hh_s_322[k]
                   + pb_z[k] * hg_228[k];

        t_323[k] = -f_4 * hf_s_158[k]
                   + f_2 * hh_s_323[k]
                   + f_5 * hf_158[k]
                   + pb_x[k] * hg_233[k];

        t_324[k] = -f_4 * hf_s_159[k]
                   + f_2 * hh_s_324[k]
                   + f_5 * hf_159[k]
                   + pb_x[k] * hg_234[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, t_328, t_329, pb_x, hh_s_325, hh_s_326, \
                         hh_s_327, hh_s_328, hh_s_329, hg_235, hg_236, hg_237, hg_238, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_2 * hh_s_325[k]
                   + pb_x[k] * hg_235[k];

        t_326[k] = f_2 * hh_s_326[k]
                   + pb_x[k] * hg_236[k];

        t_327[k] = f_2 * hh_s_327[k]
                   + pb_x[k] * hg_237[k];

        t_328[k] = f_2 * hh_s_328[k]
                   + pb_x[k] * hg_238[k];

        t_329[k] = f_2 * hh_s_329[k]
                   + pb_x[k] * hg_239[k];
    }

#pragma omp simd aligned(t_330, t_331, t_332, pb_y, pb_z, gg_160, hf_s_156, hh_s_330, \
                         hh_s_331, hh_s_332, hf_156, hg_235, hg_236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_330[k] = f_0 * gg_160[k]
                   - f_1 * hf_s_156[k]
                   + f_2 * hh_s_330[k]
                   + f_3 * hf_156[k]
                   + pb_y[k] * hg_235[k];

        t_331[k] = f_2 * hh_s_331[k]
                   + pb_z[k] * hg_235[k];

        t_332[k] = -f_4 * hf_s_156[k]
                   + f_2 * hh_s_332[k]
                   + f_5 * hf_156[k]
                   + pb_z[k] * hg_236[k];
    }

#pragma omp simd aligned(t_333, t_334, t_335, pb_y, pb_z, gg_164, hf_s_157, hf_s_159, \
                         hh_s_333, hh_s_334, hh_s_335, hf_157, hf_159, hg_237, \
                         hg_239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_333[k] = -f_6 * hf_s_157[k]
                   + f_2 * hh_s_333[k]
                   + f_7 * hf_157[k]
                   + pb_z[k] * hg_237[k];

        t_334[k] = f_0 * gg_164[k]
                   + f_2 * hh_s_334[k]
                   + pb_y[k] * hg_239[k];

        t_335[k] = -f_1 * hf_s_159[k]
                   + f_2 * hh_s_335[k]
                   + f_3 * hf_159[k]
                   + pb_z[k] * hg_239[k];
    }

#pragma omp simd aligned(t_336, t_337, t_338, t_339, pa_z, pb_x, gh_210, gh_211, gh_213, \
                         hf_s_162, hh_s_336, hh_s_337, hh_s_338, hh_s_339, hf_162, \
                         hg_242 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_336[k] = pa_z[k] * gh_210[k]
                   + f_2 * hh_s_336[k];

        t_337[k] = pa_z[k] * gh_211[k]
                   + f_2 * hh_s_337[k];

        t_338[k] = -f_10 * hf_s_162[k]
                   + f_2 * hh_s_338[k]
                   + f_8 * hf_162[k]
                   + pb_x[k] * hg_242[k];

        t_339[k] = pa_z[k] * gh_213[k]
                   + f_2 * hh_s_339[k];
    }

#pragma omp simd aligned(t_340, t_341, t_342, pa_z, pb_x, gg_151, gh_214, gh_216, hf_s_165, \
                         hh_s_340, hh_s_341, hh_s_342, hf_165, hg_245 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_340[k] = f_5 * gg_151[k]
                   + pa_z[k] * gh_214[k]
                   + f_2 * hh_s_340[k];

        t_341[k] = -f_6 * hf_s_165[k]
                   + f_2 * hh_s_341[k]
                   + f_7 * hf_165[k]
                   + pb_x[k] * hg_245[k];

        t_342[k] = pa_z[k] * gh_216[k]
                   + f_2 * hh_s_342[k];
    }

#pragma omp simd aligned(t_343, t_344, t_345, pa_z, pb_x, gg_153, gg_154, gh_217, gh_218, \
                         hf_s_169, hh_s_343, hh_s_344, hh_s_345, hf_169, \
                         hg_249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_343[k] = f_5 * gg_153[k]
                   + pa_z[k] * gh_217[k]
                   + f_2 * hh_s_343[k];

        t_344[k] = f_7 * gg_154[k]
                   + pa_z[k] * gh_218[k]
                   + f_2 * hh_s_344[k];

        t_345[k] = -f_4 * hf_s_169[k]
                   + f_2 * hh_s_345[k]
                   + f_5 * hf_169[k]
                   + pb_x[k] * hg_249[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, t_349, t_350, pb_x, hh_s_346, hh_s_347, \
                         hh_s_348, hh_s_349, hh_s_350, hg_250, hg_251, hg_252, hg_253, \
                         hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = f_2 * hh_s_346[k]
                   + pb_x[k] * hg_250[k];

        t_347[k] = f_2 * hh_s_347[k]
                   + pb_x[k] * hg_251[k];

        t_348[k] = f_2 * hh_s_348[k]
                   + pb_x[k] * hg_252[k];

        t_349[k] = f_2 * hh_s_349[k]
                   + pb_x[k] * hg_253[k];

        t_350[k] = f_2 * hh_s_350[k]
                   + pb_x[k] * hg_254[k];
    }

#pragma omp simd aligned(t_351, t_352, t_353, pa_z, pb_z, gg_160, gg_161, gh_225, gh_227, \
                         hh_s_351, hh_s_352, hh_s_353, hg_250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_351[k] = pa_z[k] * gh_225[k]
                   + f_2 * hh_s_351[k];

        t_352[k] = f_5 * gg_160[k]
                   + f_2 * hh_s_352[k]
                   + pb_z[k] * hg_250[k];

        t_353[k] = f_7 * gg_161[k]
                   + pa_z[k] * gh_227[k]
                   + f_2 * hh_s_353[k];
    }

#pragma omp simd aligned(t_354, t_355, t_356, pa_y, pa_z, pb_y, fh_s_167, fh_167, gg_162, \
                         gg_179, gh_228, gh_251, hh_s_354, hh_s_355, hh_s_356, \
                         hg_254 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_354[k] = f_8 * gg_162[k]
                   + pa_z[k] * gh_228[k]
                   + f_2 * hh_s_354[k];

        t_355[k] = f_3 * gg_179[k]
                   + f_2 * hh_s_355[k]
                   + pb_y[k] * hg_254[k];

        t_356[k] = -f_9 * fh_s_167[k]
                   + f_8 * fh_167[k]
                   + pa_y[k] * gh_251[k]
                   + f_2 * hh_s_356[k];
    }

#pragma omp simd aligned(t_357, t_358, t_359, pb_x, hf_s_170, hf_s_171, hf_s_172, hh_s_357, \
                         hh_s_358, hh_s_359, hf_170, hf_171, hf_172, hg_255, hg_256, \
                         hg_257 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_357[k] = -f_1 * hf_s_170[k]
                   + f_2 * hh_s_357[k]
                   + f_3 * hf_170[k]
                   + pb_x[k] * hg_255[k];

        t_358[k] = -f_10 * hf_s_171[k]
                   + f_2 * hh_s_358[k]
                   + f_8 * hf_171[k]
                   + pb_x[k] * hg_256[k];

        t_359[k] = -f_10 * hf_s_172[k]
                   + f_2 * hh_s_359[k]
                   + f_8 * hf_172[k]
                   + pb_x[k] * hg_257[k];
    }

#pragma omp simd aligned(t_360, t_361, t_362, pb_x, hf_s_173, hf_s_174, hf_s_175, hh_s_360, \
                         hh_s_361, hh_s_362, hf_173, hf_174, hf_175, hg_258, hg_259, \
                         hg_260 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_360[k] = -f_6 * hf_s_173[k]
                   + f_2 * hh_s_360[k]
                   + f_7 * hf_173[k]
                   + pb_x[k] * hg_258[k];

        t_361[k] = -f_6 * hf_s_174[k]
                   + f_2 * hh_s_361[k]
                   + f_7 * hf_174[k]
                   + pb_x[k] * hg_259[k];

        t_362[k] = -f_6 * hf_s_175[k]
                   + f_2 * hh_s_362[k]
                   + f_7 * hf_175[k]
                   + pb_x[k] * hg_260[k];
    }

#pragma omp simd aligned(t_363, t_364, t_365, pb_x, hf_s_176, hf_s_177, hf_s_178, hh_s_363, \
                         hh_s_364, hh_s_365, hf_176, hf_177, hf_178, hg_261, hg_262, \
                         hg_263 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_363[k] = -f_4 * hf_s_176[k]
                   + f_2 * hh_s_363[k]
                   + f_5 * hf_176[k]
                   + pb_x[k] * hg_261[k];

        t_364[k] = -f_4 * hf_s_177[k]
                   + f_2 * hh_s_364[k]
                   + f_5 * hf_177[k]
                   + pb_x[k] * hg_262[k];

        t_365[k] = -f_4 * hf_s_178[k]
                   + f_2 * hh_s_365[k]
                   + f_5 * hf_178[k]
                   + pb_x[k] * hg_263[k];
    }

#pragma omp simd aligned(t_366, t_367, t_368, t_369, pb_x, hf_s_179, hh_s_366, hh_s_367, \
                         hh_s_368, hh_s_369, hf_179, hg_264, hg_265, hg_266, \
                         hg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_366[k] = -f_4 * hf_s_179[k]
                   + f_2 * hh_s_366[k]
                   + f_5 * hf_179[k]
                   + pb_x[k] * hg_264[k];

        t_367[k] = f_2 * hh_s_367[k]
                   + pb_x[k] * hg_265[k];

        t_368[k] = f_2 * hh_s_368[k]
                   + pb_x[k] * hg_266[k];

        t_369[k] = f_2 * hh_s_369[k]
                   + pb_x[k] * hg_267[k];
    }

#pragma omp simd aligned(t_370, t_371, t_372, pa_z, pb_x, fh_s_141, fh_141, gh_246, hh_s_370, \
                         hh_s_371, hh_s_372, hg_268, hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_370[k] = f_2 * hh_s_370[k]
                   + pb_x[k] * hg_268[k];

        t_371[k] = f_2 * hh_s_371[k]
                   + pb_x[k] * hg_269[k];

        t_372[k] = -f_11 * fh_s_141[k]
                   + f_5 * fh_141[k]
                   + pa_z[k] * gh_246[k]
                   + f_2 * hh_s_372[k];
    }

#pragma omp simd aligned(t_373, t_374, pb_y, pb_z, gg_175, gg_192, hf_s_178, hh_s_373, \
                         hh_s_374, hf_178, hg_265, hg_267 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_373[k] = f_7 * gg_175[k]
                   + f_2 * hh_s_373[k]
                   + pb_z[k] * hg_265[k];

        t_374[k] = f_8 * gg_192[k]
                   - f_6 * hf_s_178[k]
                   + f_2 * hh_s_374[k]
                   + f_7 * hf_178[k]
                   + pb_y[k] * hg_267[k];
    }

#pragma omp simd aligned(t_375, t_376, pb_y, gg_193, gg_194, hf_s_179, hh_s_375, hh_s_376, \
                         hf_179, hg_268, hg_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_375[k] = f_8 * gg_193[k]
                   - f_4 * hf_s_179[k]
                   + f_2 * hh_s_375[k]
                   + f_5 * hf_179[k]
                   + pb_y[k] * hg_268[k];

        t_376[k] = f_8 * gg_194[k]
                   + f_2 * hh_s_376[k]
                   + pb_y[k] * hg_269[k];
    }

#pragma omp simd aligned(t_377, t_378, pa_y, pb_x, fh_s_188, fh_188, gh_272, hf_s_180, \
                         hh_s_377, hh_s_378, hf_180, hg_270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_377[k] = -f_12 * fh_s_188[k]
                   + f_7 * fh_188[k]
                   + pa_y[k] * gh_272[k]
                   + f_2 * hh_s_377[k];

        t_378[k] = -f_1 * hf_s_180[k]
                   + f_2 * hh_s_378[k]
                   + f_3 * hf_180[k]
                   + pb_x[k] * hg_270[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pb_x, hf_s_181, hf_s_182, hf_s_183, hh_s_379, \
                         hh_s_380, hh_s_381, hf_181, hf_182, hf_183, hg_271, hg_272, \
                         hg_273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -f_10 * hf_s_181[k]
                   + f_2 * hh_s_379[k]
                   + f_8 * hf_181[k]
                   + pb_x[k] * hg_271[k];

        t_380[k] = -f_10 * hf_s_182[k]
                   + f_2 * hh_s_380[k]
                   + f_8 * hf_182[k]
                   + pb_x[k] * hg_272[k];

        t_381[k] = -f_6 * hf_s_183[k]
                   + f_2 * hh_s_381[k]
                   + f_7 * hf_183[k]
                   + pb_x[k] * hg_273[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, pb_x, hf_s_184, hf_s_185, hf_s_186, hh_s_382, \
                         hh_s_383, hh_s_384, hf_184, hf_185, hf_186, hg_274, hg_275, \
                         hg_276 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = -f_6 * hf_s_184[k]
                   + f_2 * hh_s_382[k]
                   + f_7 * hf_184[k]
                   + pb_x[k] * hg_274[k];

        t_383[k] = -f_6 * hf_s_185[k]
                   + f_2 * hh_s_383[k]
                   + f_7 * hf_185[k]
                   + pb_x[k] * hg_275[k];

        t_384[k] = -f_4 * hf_s_186[k]
                   + f_2 * hh_s_384[k]
                   + f_5 * hf_186[k]
                   + pb_x[k] * hg_276[k];
    }

#pragma omp simd aligned(t_385, t_386, t_387, pb_x, hf_s_187, hf_s_188, hf_s_189, hh_s_385, \
                         hh_s_386, hh_s_387, hf_187, hf_188, hf_189, hg_277, hg_278, \
                         hg_279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_385[k] = -f_4 * hf_s_187[k]
                   + f_2 * hh_s_385[k]
                   + f_5 * hf_187[k]
                   + pb_x[k] * hg_277[k];

        t_386[k] = -f_4 * hf_s_188[k]
                   + f_2 * hh_s_386[k]
                   + f_5 * hf_188[k]
                   + pb_x[k] * hg_278[k];

        t_387[k] = -f_4 * hf_s_189[k]
                   + f_2 * hh_s_387[k]
                   + f_5 * hf_189[k]
                   + pb_x[k] * hg_279[k];
    }

#pragma omp simd aligned(t_388, t_389, t_390, t_391, t_392, pb_x, hh_s_388, hh_s_389, \
                         hh_s_390, hh_s_391, hh_s_392, hg_280, hg_281, hg_282, hg_283, \
                         hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_388[k] = f_2 * hh_s_388[k]
                   + pb_x[k] * hg_280[k];

        t_389[k] = f_2 * hh_s_389[k]
                   + pb_x[k] * hg_281[k];

        t_390[k] = f_2 * hh_s_390[k]
                   + pb_x[k] * hg_282[k];

        t_391[k] = f_2 * hh_s_391[k]
                   + pb_x[k] * hg_283[k];

        t_392[k] = f_2 * hh_s_392[k]
                   + pb_x[k] * hg_284[k];
    }

#pragma omp simd aligned(t_393, t_394, pa_z, pb_z, fh_s_162, fh_162, gg_190, gh_267, hh_s_393, \
                         hh_s_394, hg_280 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = -f_12 * fh_s_162[k]
                   + f_7 * fh_162[k]
                   + pa_z[k] * gh_267[k]
                   + f_2 * hh_s_393[k];

        t_394[k] = f_8 * gg_190[k]
                   + f_2 * hh_s_394[k]
                   + pb_z[k] * hg_280[k];
    }

#pragma omp simd aligned(t_395, t_396, pb_y, gg_207, gg_208, hf_s_188, hf_s_189, hh_s_395, \
                         hh_s_396, hf_188, hf_189, hg_282, hg_283 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_395[k] = f_7 * gg_207[k]
                   - f_6 * hf_s_188[k]
                   + f_2 * hh_s_395[k]
                   + f_7 * hf_188[k]
                   + pb_y[k] * hg_282[k];

        t_396[k] = f_7 * gg_208[k]
                   - f_4 * hf_s_189[k]
                   + f_2 * hh_s_396[k]
                   + f_5 * hf_189[k]
                   + pb_y[k] * hg_283[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, pa_y, pb_y, fh_s_209, fh_209, gg_209, gh_293, \
                         gh_294, hh_s_397, hh_s_398, hh_s_399, hg_284 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_7 * gg_209[k]
                   + f_2 * hh_s_397[k]
                   + pb_y[k] * hg_284[k];

        t_398[k] = -f_11 * fh_s_209[k]
                   + f_5 * fh_209[k]
                   + pa_y[k] * gh_293[k]
                   + f_2 * hh_s_398[k];

        t_399[k] = pa_y[k] * gh_294[k]
                   + f_2 * hh_s_399[k];
    }

#pragma omp simd aligned(t_400, t_401, t_402, t_403, pa_y, gg_210, gg_211, gg_212, gh_295, \
                         gh_296, gh_297, gh_298, hh_s_400, hh_s_401, hh_s_402, \
                         hh_s_403 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_400[k] = f_5 * gg_210[k]
                   + pa_y[k] * gh_295[k]
                   + f_2 * hh_s_400[k];

        t_401[k] = pa_y[k] * gh_296[k]
                   + f_2 * hh_s_401[k];

        t_402[k] = f_7 * gg_211[k]
                   + pa_y[k] * gh_297[k]
                   + f_2 * hh_s_402[k];

        t_403[k] = f_5 * gg_212[k]
                   + pa_y[k] * gh_298[k]
                   + f_2 * hh_s_403[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, t_407, pa_y, gg_213, gg_214, gg_215, gh_299, \
                         gh_300, gh_301, gh_302, hh_s_404, hh_s_405, hh_s_406, \
                         hh_s_407 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * gh_299[k]
                   + f_2 * hh_s_404[k];

        t_405[k] = f_8 * gg_213[k]
                   + pa_y[k] * gh_300[k]
                   + f_2 * hh_s_405[k];

        t_406[k] = f_7 * gg_214[k]
                   + pa_y[k] * gh_301[k]
                   + f_2 * hh_s_406[k];

        t_407[k] = f_5 * gg_215[k]
                   + pa_y[k] * gh_302[k]
                   + f_2 * hh_s_407[k];
    }

#pragma omp simd aligned(t_408, t_409, t_410, t_411, pa_y, pb_x, gh_303, hh_s_408, hh_s_409, \
                         hh_s_410, hh_s_411, hg_295, hg_296, hg_297 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_408[k] = pa_y[k] * gh_303[k]
                   + f_2 * hh_s_408[k];

        t_409[k] = f_2 * hh_s_409[k]
                   + pb_x[k] * hg_295[k];

        t_410[k] = f_2 * hh_s_410[k]
                   + pb_x[k] * hg_296[k];

        t_411[k] = f_2 * hh_s_411[k]
                   + pb_x[k] * hg_297[k];
    }

#pragma omp simd aligned(t_412, t_413, t_414, pa_y, pb_x, gg_220, gh_309, hh_s_412, hh_s_413, \
                         hh_s_414, hg_298, hg_299 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_412[k] = f_2 * hh_s_412[k]
                   + pb_x[k] * hg_298[k];

        t_413[k] = f_2 * hh_s_413[k]
                   + pb_x[k] * hg_299[k];

        t_414[k] = f_0 * gg_220[k]
                   + pa_y[k] * gh_309[k]
                   + f_2 * hh_s_414[k];
    }

#pragma omp simd aligned(t_415, t_416, t_417, pa_y, pb_z, gg_205, gg_222, gg_223, gh_311, \
                         gh_312, hh_s_415, hh_s_416, hh_s_417, hg_295 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_415[k] = f_3 * gg_205[k]
                   + f_2 * hh_s_415[k]
                   + pb_z[k] * hg_295[k];

        t_416[k] = f_8 * gg_222[k]
                   + pa_y[k] * gh_311[k]
                   + f_2 * hh_s_416[k];

        t_417[k] = f_7 * gg_223[k]
                   + pa_y[k] * gh_312[k]
                   + f_2 * hh_s_417[k];
    }

#pragma omp simd aligned(t_418, t_419, t_420, pa_y, pb_x, pb_y, gg_224, gh_314, hf_s_200, \
                         hh_s_418, hh_s_419, hh_s_420, hf_200, hg_299, \
                         hg_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_418[k] = f_5 * gg_224[k]
                   + f_2 * hh_s_418[k]
                   + pb_y[k] * hg_299[k];

        t_419[k] = pa_y[k] * gh_314[k]
                   + f_2 * hh_s_419[k];

        t_420[k] = -f_1 * hf_s_200[k]
                   + f_2 * hh_s_420[k]
                   + f_3 * hf_200[k]
                   + pb_x[k] * hg_300[k];
    }
}

static auto
compute_prim_hh_kinetic_energy_0_piece4(CSimdMatrix &buffer, const size_t target,
                                        const size_t pb, const size_t gg, const size_t hf_s,
                                        const size_t hh_s, const size_t hf, const size_t hg,
                                        const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 2.5 / p;
    const auto f_1 = 4.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 2.0 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 2.0 * alpha / p;
    const auto f_7 = 1.0 / p;
    const auto f_8 = 1.5 / p;
    const auto f_10 = 3.0 * alpha / p;

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

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_224 = buffer.data(gg + 224);

    const auto *hf_s_202 = buffer.data(hf_s + 202);
    const auto *hf_s_203 = buffer.data(hf_s + 203);
    const auto *hf_s_205 = buffer.data(hf_s + 205);
    const auto *hf_s_206 = buffer.data(hf_s + 206);
    const auto *hf_s_207 = buffer.data(hf_s + 207);
    const auto *hf_s_208 = buffer.data(hf_s + 208);
    const auto *hf_s_209 = buffer.data(hf_s + 209);

    const auto *hh_s_421 = buffer.data(hh_s + 421);
    const auto *hh_s_422 = buffer.data(hh_s + 422);
    const auto *hh_s_423 = buffer.data(hh_s + 423);
    const auto *hh_s_424 = buffer.data(hh_s + 424);
    const auto *hh_s_425 = buffer.data(hh_s + 425);
    const auto *hh_s_426 = buffer.data(hh_s + 426);
    const auto *hh_s_427 = buffer.data(hh_s + 427);
    const auto *hh_s_428 = buffer.data(hh_s + 428);
    const auto *hh_s_429 = buffer.data(hh_s + 429);
    const auto *hh_s_430 = buffer.data(hh_s + 430);
    const auto *hh_s_431 = buffer.data(hh_s + 431);
    const auto *hh_s_432 = buffer.data(hh_s + 432);
    const auto *hh_s_433 = buffer.data(hh_s + 433);
    const auto *hh_s_434 = buffer.data(hh_s + 434);
    const auto *hh_s_435 = buffer.data(hh_s + 435);
    const auto *hh_s_436 = buffer.data(hh_s + 436);
    const auto *hh_s_437 = buffer.data(hh_s + 437);
    const auto *hh_s_438 = buffer.data(hh_s + 438);
    const auto *hh_s_439 = buffer.data(hh_s + 439);
    const auto *hh_s_440 = buffer.data(hh_s + 440);

    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_306 = buffer.data(hg + 306);
    const auto *hg_307 = buffer.data(hg + 307);
    const auto *hg_309 = buffer.data(hg + 309);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_313 = buffer.data(hg + 313);
    const auto *hg_314 = buffer.data(hg + 314);

#pragma omp simd aligned(t_421, t_422, t_423, pb_x, pb_y, hf_s_202, hf_s_203, hh_s_421, \
                         hh_s_422, hh_s_423, hf_202, hf_203, hg_300, hg_302, \
                         hg_303 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_421[k] = f_2 * hh_s_421[k]
                   + pb_y[k] * hg_300[k];

        t_422[k] = -f_10 * hf_s_202[k]
                   + f_2 * hh_s_422[k]
                   + f_8 * hf_202[k]
                   + pb_x[k] * hg_302[k];

        t_423[k] = -f_6 * hf_s_203[k]
                   + f_2 * hh_s_423[k]
                   + f_7 * hf_203[k]
                   + pb_x[k] * hg_303[k];
    }

#pragma omp simd aligned(t_424, t_425, t_426, pb_x, pb_y, hf_s_205, hf_s_206, hh_s_424, \
                         hh_s_425, hh_s_426, hf_205, hf_206, hg_302, hg_305, \
                         hg_306 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_424[k] = f_2 * hh_s_424[k]
                   + pb_y[k] * hg_302[k];

        t_425[k] = -f_6 * hf_s_205[k]
                   + f_2 * hh_s_425[k]
                   + f_7 * hf_205[k]
                   + pb_x[k] * hg_305[k];

        t_426[k] = -f_4 * hf_s_206[k]
                   + f_2 * hh_s_426[k]
                   + f_5 * hf_206[k]
                   + pb_x[k] * hg_306[k];
    }

#pragma omp simd aligned(t_427, t_428, t_429, pb_x, pb_y, hf_s_207, hf_s_209, hh_s_427, \
                         hh_s_428, hh_s_429, hf_207, hf_209, hg_305, hg_307, \
                         hg_309 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_427[k] = -f_4 * hf_s_207[k]
                   + f_2 * hh_s_427[k]
                   + f_5 * hf_207[k]
                   + pb_x[k] * hg_307[k];

        t_428[k] = f_2 * hh_s_428[k]
                   + pb_y[k] * hg_305[k];

        t_429[k] = -f_4 * hf_s_209[k]
                   + f_2 * hh_s_429[k]
                   + f_5 * hf_209[k]
                   + pb_x[k] * hg_309[k];
    }

#pragma omp simd aligned(t_430, t_431, t_432, t_433, t_434, pb_x, hh_s_430, hh_s_431, \
                         hh_s_432, hh_s_433, hh_s_434, hg_310, hg_311, hg_312, hg_313, \
                         hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_430[k] = f_2 * hh_s_430[k]
                   + pb_x[k] * hg_310[k];

        t_431[k] = f_2 * hh_s_431[k]
                   + pb_x[k] * hg_311[k];

        t_432[k] = f_2 * hh_s_432[k]
                   + pb_x[k] * hg_312[k];

        t_433[k] = f_2 * hh_s_433[k]
                   + pb_x[k] * hg_313[k];

        t_434[k] = f_2 * hh_s_434[k]
                   + pb_x[k] * hg_314[k];
    }

#pragma omp simd aligned(t_435, t_436, t_437, pb_y, hf_s_206, hf_s_207, hf_s_208, hh_s_435, \
                         hh_s_436, hh_s_437, hf_206, hf_207, hf_208, hg_310, hg_311, \
                         hg_312 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_435[k] = -f_1 * hf_s_206[k]
                   + f_2 * hh_s_435[k]
                   + f_3 * hf_206[k]
                   + pb_y[k] * hg_310[k];

        t_436[k] = -f_10 * hf_s_207[k]
                   + f_2 * hh_s_436[k]
                   + f_8 * hf_207[k]
                   + pb_y[k] * hg_311[k];

        t_437[k] = -f_6 * hf_s_208[k]
                   + f_2 * hh_s_437[k]
                   + f_7 * hf_208[k]
                   + pb_y[k] * hg_312[k];
    }

#pragma omp simd aligned(t_438, t_439, t_440, pb_y, pb_z, gg_224, hf_s_209, hh_s_438, \
                         hh_s_439, hh_s_440, hf_209, hg_313, hg_314 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_438[k] = -f_4 * hf_s_209[k]
                   + f_2 * hh_s_438[k]
                   + f_5 * hf_209[k]
                   + pb_y[k] * hg_313[k];

        t_439[k] = f_2 * hh_s_439[k]
                   + pb_y[k] * hg_314[k];

        t_440[k] = f_0 * gg_224[k]
                   - f_1 * hf_s_209[k]
                   + f_2 * hh_s_440[k]
                   + f_3 * hf_209[k]
                   + pb_z[k] * hg_314[k];
    }
}

auto
compute_prim_hh_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t fh_s, const size_t fh,
                                 const size_t gg, const size_t gh, const size_t hf_s,
                                 const size_t hh_s, const size_t hf, const size_t hg,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_hh_kinetic_energy_0_piece0(buffer, target, pa, pb, fh_s, fh, gg, gh, hf_s,
                                            hh_s, hf, hg, ncols, alpha, beta, p);

    compute_prim_hh_kinetic_energy_0_piece1(buffer, target, pa, pb, fh_s, fh, gg, gh, hf_s,
                                            hh_s, hf, hg, ncols, alpha, beta, p);

    compute_prim_hh_kinetic_energy_0_piece2(buffer, target, pa, pb, fh_s, fh, gg, gh, hf_s,
                                            hh_s, hf, hg, ncols, alpha, beta, p);

    compute_prim_hh_kinetic_energy_0_piece3(buffer, target, pa, pb, fh_s, fh, gg, gh, hf_s,
                                            hh_s, hf, hg, ncols, alpha, beta, p);

    compute_prim_hh_kinetic_energy_0_piece4(buffer, target, pb, gg, hf_s, hh_s, hf, hg, ncols,
                                            alpha, beta, p);
}

}  // namespace simdkin
