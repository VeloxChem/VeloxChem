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


#include "SimdKineticEnergyVrrRecIG.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_ig_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gg_s,
                                        const size_t gg, const size_t hf, const size_t hg,
                                        const size_t id_s, const size_t ig_s, const size_t id,
                                        const size_t if_, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.5 / p;
    const auto f_8 = 4.0 * beta / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;
    const auto f_13 = 2.0 * beta / p;

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

    const auto *gg_s_0 = buffer.data(gg_s + 0);
    const auto *gg_s_15 = buffer.data(gg_s + 15);
    const auto *gg_s_25 = buffer.data(gg_s + 25);
    const auto *gg_s_44 = buffer.data(gg_s + 44);
    const auto *gg_s_55 = buffer.data(gg_s + 55);
    const auto *gg_s_72 = buffer.data(gg_s + 72);
    const auto *gg_s_89 = buffer.data(gg_s + 89);
    const auto *gg_s_100 = buffer.data(gg_s + 100);

    const auto *gg_0 = buffer.data(gg + 0);
    const auto *gg_15 = buffer.data(gg + 15);
    const auto *gg_25 = buffer.data(gg + 25);
    const auto *gg_44 = buffer.data(gg + 44);
    const auto *gg_55 = buffer.data(gg + 55);
    const auto *gg_72 = buffer.data(gg + 72);
    const auto *gg_89 = buffer.data(gg + 89);
    const auto *gg_100 = buffer.data(gg + 100);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_1 = buffer.data(hf + 1);
    const auto *hf_2 = buffer.data(hf + 2);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_18 = buffer.data(hf + 18);
    const auto *hf_19 = buffer.data(hf + 19);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_27 = buffer.data(hf + 27);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_38 = buffer.data(hf + 38);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_47 = buffer.data(hf + 47);
    const auto *hf_48 = buffer.data(hf + 48);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_56 = buffer.data(hf + 56);
    const auto *hf_57 = buffer.data(hf + 57);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_68 = buffer.data(hf + 68);
    const auto *hf_69 = buffer.data(hf + 69);

    const auto *hg_0 = buffer.data(hg + 0);
    const auto *hg_3 = buffer.data(hg + 3);
    const auto *hg_5 = buffer.data(hg + 5);
    const auto *hg_6 = buffer.data(hg + 6);
    const auto *hg_9 = buffer.data(hg + 9);
    const auto *hg_10 = buffer.data(hg + 10);
    const auto *hg_14 = buffer.data(hg + 14);
    const auto *hg_15 = buffer.data(hg + 15);
    const auto *hg_16 = buffer.data(hg + 16);
    const auto *hg_18 = buffer.data(hg + 18);
    const auto *hg_21 = buffer.data(hg + 21);
    const auto *hg_25 = buffer.data(hg + 25);
    const auto *hg_30 = buffer.data(hg + 30);
    const auto *hg_32 = buffer.data(hg + 32);
    const auto *hg_35 = buffer.data(hg + 35);
    const auto *hg_39 = buffer.data(hg + 39);
    const auto *hg_44 = buffer.data(hg + 44);
    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_72 = buffer.data(hg + 72);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_100 = buffer.data(hg + 100);

    const auto *id_s_0 = buffer.data(id_s + 0);
    const auto *id_s_3 = buffer.data(id_s + 3);
    const auto *id_s_5 = buffer.data(id_s + 5);
    const auto *id_s_9 = buffer.data(id_s + 9);
    const auto *id_s_16 = buffer.data(id_s + 16);
    const auto *id_s_17 = buffer.data(id_s + 17);
    const auto *id_s_18 = buffer.data(id_s + 18);
    const auto *id_s_21 = buffer.data(id_s + 21);
    const auto *id_s_23 = buffer.data(id_s + 23);
    const auto *id_s_30 = buffer.data(id_s + 30);
    const auto *id_s_33 = buffer.data(id_s + 33);
    const auto *id_s_34 = buffer.data(id_s + 34);
    const auto *id_s_35 = buffer.data(id_s + 35);
    const auto *id_s_36 = buffer.data(id_s + 36);
    const auto *id_s_39 = buffer.data(id_s + 39);

    const auto *ig_s_0 = buffer.data(ig_s + 0);
    const auto *ig_s_1 = buffer.data(ig_s + 1);
    const auto *ig_s_2 = buffer.data(ig_s + 2);
    const auto *ig_s_3 = buffer.data(ig_s + 3);
    const auto *ig_s_4 = buffer.data(ig_s + 4);
    const auto *ig_s_5 = buffer.data(ig_s + 5);
    const auto *ig_s_6 = buffer.data(ig_s + 6);
    const auto *ig_s_7 = buffer.data(ig_s + 7);
    const auto *ig_s_8 = buffer.data(ig_s + 8);
    const auto *ig_s_9 = buffer.data(ig_s + 9);
    const auto *ig_s_10 = buffer.data(ig_s + 10);
    const auto *ig_s_11 = buffer.data(ig_s + 11);
    const auto *ig_s_12 = buffer.data(ig_s + 12);
    const auto *ig_s_13 = buffer.data(ig_s + 13);
    const auto *ig_s_14 = buffer.data(ig_s + 14);
    const auto *ig_s_15 = buffer.data(ig_s + 15);
    const auto *ig_s_16 = buffer.data(ig_s + 16);
    const auto *ig_s_17 = buffer.data(ig_s + 17);
    const auto *ig_s_18 = buffer.data(ig_s + 18);
    const auto *ig_s_19 = buffer.data(ig_s + 19);
    const auto *ig_s_20 = buffer.data(ig_s + 20);
    const auto *ig_s_21 = buffer.data(ig_s + 21);
    const auto *ig_s_22 = buffer.data(ig_s + 22);
    const auto *ig_s_23 = buffer.data(ig_s + 23);
    const auto *ig_s_24 = buffer.data(ig_s + 24);
    const auto *ig_s_25 = buffer.data(ig_s + 25);
    const auto *ig_s_26 = buffer.data(ig_s + 26);
    const auto *ig_s_27 = buffer.data(ig_s + 27);
    const auto *ig_s_28 = buffer.data(ig_s + 28);
    const auto *ig_s_29 = buffer.data(ig_s + 29);
    const auto *ig_s_30 = buffer.data(ig_s + 30);
    const auto *ig_s_31 = buffer.data(ig_s + 31);
    const auto *ig_s_32 = buffer.data(ig_s + 32);
    const auto *ig_s_33 = buffer.data(ig_s + 33);
    const auto *ig_s_34 = buffer.data(ig_s + 34);
    const auto *ig_s_35 = buffer.data(ig_s + 35);
    const auto *ig_s_36 = buffer.data(ig_s + 36);
    const auto *ig_s_37 = buffer.data(ig_s + 37);
    const auto *ig_s_38 = buffer.data(ig_s + 38);
    const auto *ig_s_39 = buffer.data(ig_s + 39);
    const auto *ig_s_40 = buffer.data(ig_s + 40);
    const auto *ig_s_41 = buffer.data(ig_s + 41);
    const auto *ig_s_42 = buffer.data(ig_s + 42);
    const auto *ig_s_43 = buffer.data(ig_s + 43);
    const auto *ig_s_44 = buffer.data(ig_s + 44);
    const auto *ig_s_45 = buffer.data(ig_s + 45);
    const auto *ig_s_46 = buffer.data(ig_s + 46);
    const auto *ig_s_47 = buffer.data(ig_s + 47);
    const auto *ig_s_48 = buffer.data(ig_s + 48);
    const auto *ig_s_49 = buffer.data(ig_s + 49);
    const auto *ig_s_50 = buffer.data(ig_s + 50);
    const auto *ig_s_51 = buffer.data(ig_s + 51);
    const auto *ig_s_52 = buffer.data(ig_s + 52);
    const auto *ig_s_53 = buffer.data(ig_s + 53);
    const auto *ig_s_54 = buffer.data(ig_s + 54);
    const auto *ig_s_55 = buffer.data(ig_s + 55);
    const auto *ig_s_56 = buffer.data(ig_s + 56);
    const auto *ig_s_57 = buffer.data(ig_s + 57);
    const auto *ig_s_58 = buffer.data(ig_s + 58);
    const auto *ig_s_59 = buffer.data(ig_s + 59);
    const auto *ig_s_60 = buffer.data(ig_s + 60);
    const auto *ig_s_61 = buffer.data(ig_s + 61);
    const auto *ig_s_62 = buffer.data(ig_s + 62);
    const auto *ig_s_63 = buffer.data(ig_s + 63);
    const auto *ig_s_64 = buffer.data(ig_s + 64);
    const auto *ig_s_65 = buffer.data(ig_s + 65);
    const auto *ig_s_66 = buffer.data(ig_s + 66);
    const auto *ig_s_67 = buffer.data(ig_s + 67);
    const auto *ig_s_68 = buffer.data(ig_s + 68);
    const auto *ig_s_69 = buffer.data(ig_s + 69);
    const auto *ig_s_70 = buffer.data(ig_s + 70);
    const auto *ig_s_71 = buffer.data(ig_s + 71);
    const auto *ig_s_72 = buffer.data(ig_s + 72);
    const auto *ig_s_73 = buffer.data(ig_s + 73);
    const auto *ig_s_74 = buffer.data(ig_s + 74);
    const auto *ig_s_75 = buffer.data(ig_s + 75);
    const auto *ig_s_76 = buffer.data(ig_s + 76);
    const auto *ig_s_77 = buffer.data(ig_s + 77);
    const auto *ig_s_78 = buffer.data(ig_s + 78);
    const auto *ig_s_79 = buffer.data(ig_s + 79);
    const auto *ig_s_80 = buffer.data(ig_s + 80);
    const auto *ig_s_81 = buffer.data(ig_s + 81);
    const auto *ig_s_82 = buffer.data(ig_s + 82);
    const auto *ig_s_83 = buffer.data(ig_s + 83);
    const auto *ig_s_84 = buffer.data(ig_s + 84);
    const auto *ig_s_85 = buffer.data(ig_s + 85);
    const auto *ig_s_86 = buffer.data(ig_s + 86);
    const auto *ig_s_87 = buffer.data(ig_s + 87);
    const auto *ig_s_88 = buffer.data(ig_s + 88);
    const auto *ig_s_89 = buffer.data(ig_s + 89);
    const auto *ig_s_90 = buffer.data(ig_s + 90);
    const auto *ig_s_91 = buffer.data(ig_s + 91);
    const auto *ig_s_92 = buffer.data(ig_s + 92);
    const auto *ig_s_93 = buffer.data(ig_s + 93);
    const auto *ig_s_94 = buffer.data(ig_s + 94);
    const auto *ig_s_95 = buffer.data(ig_s + 95);
    const auto *ig_s_96 = buffer.data(ig_s + 96);
    const auto *ig_s_97 = buffer.data(ig_s + 97);
    const auto *ig_s_98 = buffer.data(ig_s + 98);
    const auto *ig_s_99 = buffer.data(ig_s + 99);
    const auto *ig_s_100 = buffer.data(ig_s + 100);
    const auto *ig_s_101 = buffer.data(ig_s + 101);
    const auto *ig_s_102 = buffer.data(ig_s + 102);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_39 = buffer.data(id + 39);

    const auto *if__0 = buffer.data(if_ + 0);
    const auto *if__1 = buffer.data(if_ + 1);
    const auto *if__2 = buffer.data(if_ + 2);
    const auto *if__3 = buffer.data(if_ + 3);
    const auto *if__5 = buffer.data(if_ + 5);
    const auto *if__6 = buffer.data(if_ + 6);
    const auto *if__8 = buffer.data(if_ + 8);
    const auto *if__9 = buffer.data(if_ + 9);
    const auto *if__10 = buffer.data(if_ + 10);
    const auto *if__11 = buffer.data(if_ + 11);
    const auto *if__13 = buffer.data(if_ + 13);
    const auto *if__16 = buffer.data(if_ + 16);
    const auto *if__17 = buffer.data(if_ + 17);
    const auto *if__18 = buffer.data(if_ + 18);
    const auto *if__19 = buffer.data(if_ + 19);
    const auto *if__20 = buffer.data(if_ + 20);
    const auto *if__22 = buffer.data(if_ + 22);
    const auto *if__25 = buffer.data(if_ + 25);
    const auto *if__27 = buffer.data(if_ + 27);
    const auto *if__28 = buffer.data(if_ + 28);
    const auto *if__29 = buffer.data(if_ + 29);
    const auto *if__30 = buffer.data(if_ + 30);
    const auto *if__31 = buffer.data(if_ + 31);
    const auto *if__32 = buffer.data(if_ + 32);
    const auto *if__33 = buffer.data(if_ + 33);
    const auto *if__36 = buffer.data(if_ + 36);
    const auto *if__37 = buffer.data(if_ + 37);
    const auto *if__38 = buffer.data(if_ + 38);
    const auto *if__39 = buffer.data(if_ + 39);
    const auto *if__42 = buffer.data(if_ + 42);
    const auto *if__46 = buffer.data(if_ + 46);
    const auto *if__47 = buffer.data(if_ + 47);
    const auto *if__48 = buffer.data(if_ + 48);
    const auto *if__49 = buffer.data(if_ + 49);
    const auto *if__50 = buffer.data(if_ + 50);
    const auto *if__51 = buffer.data(if_ + 51);
    const auto *if__52 = buffer.data(if_ + 52);
    const auto *if__55 = buffer.data(if_ + 55);
    const auto *if__56 = buffer.data(if_ + 56);
    const auto *if__57 = buffer.data(if_ + 57);
    const auto *if__58 = buffer.data(if_ + 58);
    const auto *if__59 = buffer.data(if_ + 59);
    const auto *if__60 = buffer.data(if_ + 60);
    const auto *if__61 = buffer.data(if_ + 61);
    const auto *if__62 = buffer.data(if_ + 62);
    const auto *if__63 = buffer.data(if_ + 63);
    const auto *if__66 = buffer.data(if_ + 66);
    const auto *if__67 = buffer.data(if_ + 67);
    const auto *if__68 = buffer.data(if_ + 68);
    const auto *if__69 = buffer.data(if_ + 69);

#pragma omp simd aligned(t_0, t_1, t_2, t_3, pb_x, pb_y, pb_z, hf_0, id_s_0, ig_s_0, ig_s_1, \
                         ig_s_2, ig_s_3, id_0, if__0, if__1 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hf_0[k]
                 - f_1 * id_s_0[k]
                 + f_2 * ig_s_0[k]
                 + f_3 * id_0[k]
                 + pb_x[k] * if__0[k];

        t_1[k] = f_2 * ig_s_1[k]
                 + pb_y[k] * if__0[k];

        t_2[k] = f_2 * ig_s_2[k]
                 + pb_z[k] * if__0[k];

        t_3[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_3[k]
                 + f_5 * id_0[k]
                 + pb_y[k] * if__1[k];
    }

#pragma omp simd aligned(t_4, t_5, t_6, pb_x, pb_y, pb_z, hf_6, id_s_0, ig_s_4, ig_s_5, \
                         ig_s_6, id_0, if__2, if__6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_4[k] = f_2 * ig_s_4[k]
                 + pb_y[k] * if__2[k];

        t_5[k] = -f_4 * id_s_0[k]
                 + f_2 * ig_s_5[k]
                 + f_5 * id_0[k]
                 + pb_z[k] * if__2[k];

        t_6[k] = f_0 * hf_6[k]
                 + f_2 * ig_s_6[k]
                 + pb_x[k] * if__6[k];
    }

#pragma omp simd aligned(t_7, t_8, t_9, pb_x, pb_y, pb_z, hf_9, ig_s_7, ig_s_8, ig_s_9, if__3, \
                         if__5, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_7[k] = f_2 * ig_s_7[k]
                 + pb_z[k] * if__3[k];

        t_8[k] = f_2 * ig_s_8[k]
                 + pb_y[k] * if__5[k];

        t_9[k] = f_0 * hf_9[k]
                 + f_2 * ig_s_9[k]
                 + pb_x[k] * if__9[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pb_y, pb_z, id_s_3, id_s_5, ig_s_10, ig_s_11, \
                         ig_s_12, id_3, id_5, if__6, if__8 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = -f_1 * id_s_3[k]
                  + f_2 * ig_s_10[k]
                  + f_3 * id_3[k]
                  + pb_y[k] * if__6[k];

        t_11[k] = f_2 * ig_s_11[k]
                  + pb_z[k] * if__6[k];

        t_12[k] = -f_4 * id_s_5[k]
                  + f_2 * ig_s_12[k]
                  + f_5 * id_5[k]
                  + pb_y[k] * if__8[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_y, pb_z, hg_0, id_s_5, ig_s_13, ig_s_14, \
                         ig_s_15, id_5, if__9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_2 * ig_s_13[k]
                  + pb_y[k] * if__9[k];

        t_14[k] = -f_1 * id_s_5[k]
                  + f_2 * ig_s_14[k]
                  + f_3 * id_5[k]
                  + pb_z[k] * if__9[k];

        t_15[k] = pa_y[k] * hg_0[k]
                  + f_2 * ig_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, t_19, pa_y, pb_y, pb_z, hf_0, hf_1, hg_3, ig_s_16, \
                         ig_s_17, ig_s_18, ig_s_19, if__10, if__11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_16[k]
                  + pb_y[k] * if__10[k];

        t_17[k] = f_2 * ig_s_17[k]
                  + pb_z[k] * if__10[k];

        t_18[k] = f_6 * hf_1[k]
                  + pa_y[k] * hg_3[k]
                  + f_2 * ig_s_18[k];

        t_19[k] = f_2 * ig_s_19[k]
                  + pb_z[k] * if__11[k];
    }

#pragma omp simd aligned(t_20, t_21, t_22, pa_y, pb_x, pb_z, hf_16, hg_5, ig_s_20, ig_s_21, \
                         ig_s_22, if__13, if__16 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_20[k] = pa_y[k] * hg_5[k]
                  + f_2 * ig_s_20[k];

        t_21[k] = f_7 * hf_16[k]
                  + f_2 * ig_s_21[k]
                  + pb_x[k] * if__16[k];

        t_22[k] = f_2 * ig_s_22[k]
                  + pb_z[k] * if__13[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, pa_x, pa_y, pb_x, gg_s_25, gg_25, hf_18, hg_9, \
                         hg_25, ig_s_23, ig_s_24, ig_s_25, if__18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = f_7 * hf_18[k]
                  + f_2 * ig_s_23[k]
                  + pb_x[k] * if__18[k];

        t_24[k] = pa_y[k] * hg_9[k]
                  + f_2 * ig_s_24[k];

        t_25[k] = -f_8 * gg_s_25[k]
                  + f_9 * gg_25[k]
                  + pa_x[k] * hg_25[k]
                  + f_2 * ig_s_25[k];
    }

#pragma omp simd aligned(t_26, t_27, t_28, pb_y, pb_z, hf_9, id_s_9, ig_s_26, ig_s_27, \
                         ig_s_28, id_9, if__16, if__17, if__19 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_26[k] = f_2 * ig_s_26[k]
                  + pb_z[k] * if__16[k];

        t_27[k] = -f_4 * id_s_9[k]
                  + f_2 * ig_s_27[k]
                  + f_5 * id_9[k]
                  + pb_z[k] * if__17[k];

        t_28[k] = f_5 * hf_9[k]
                  + f_2 * ig_s_28[k]
                  + pb_y[k] * if__19[k];
    }

#pragma omp simd aligned(t_29, t_30, t_31, t_32, pa_y, pa_z, pb_y, pb_z, hf_0, hg_0, hg_14, \
                         ig_s_29, ig_s_30, ig_s_31, ig_s_32, if__20 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_29[k] = pa_y[k] * hg_14[k]
                  + f_2 * ig_s_29[k];

        t_30[k] = pa_z[k] * hg_0[k]
                  + f_2 * ig_s_30[k];

        t_31[k] = f_2 * ig_s_31[k]
                  + pb_y[k] * if__20[k];

        t_32[k] = f_5 * hf_0[k]
                  + f_2 * ig_s_32[k]
                  + pb_z[k] * if__20[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, t_36, pa_z, pb_y, hf_2, hg_3, hg_5, hg_6, ig_s_33, \
                         ig_s_34, ig_s_35, ig_s_36, if__22 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = pa_z[k] * hg_3[k]
                  + f_2 * ig_s_33[k];

        t_34[k] = f_2 * ig_s_34[k]
                  + pb_y[k] * if__22[k];

        t_35[k] = f_6 * hf_2[k]
                  + pa_z[k] * hg_5[k]
                  + f_2 * ig_s_35[k];

        t_36[k] = pa_z[k] * hg_6[k]
                  + f_2 * ig_s_36[k];
    }

#pragma omp simd aligned(t_37, t_38, t_39, pb_x, pb_y, hf_27, hf_29, ig_s_37, ig_s_38, \
                         ig_s_39, if__25, if__27, if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_37[k] = f_7 * hf_27[k]
                  + f_2 * ig_s_37[k]
                  + pb_x[k] * if__27[k];

        t_38[k] = f_2 * ig_s_38[k]
                  + pb_y[k] * if__25[k];

        t_39[k] = f_7 * hf_29[k]
                  + f_2 * ig_s_39[k]
                  + pb_x[k] * if__29[k];
    }

#pragma omp simd aligned(t_40, t_41, t_42, pa_z, pb_y, hg_10, id_s_16, id_s_17, ig_s_40, \
                         ig_s_41, ig_s_42, id_16, id_17, if__27, \
                         if__28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_40[k] = pa_z[k] * hg_10[k]
                  + f_2 * ig_s_40[k];

        t_41[k] = -f_10 * id_s_16[k]
                  + f_2 * ig_s_41[k]
                  + f_6 * id_16[k]
                  + pb_y[k] * if__27[k];

        t_42[k] = -f_4 * id_s_17[k]
                  + f_2 * ig_s_42[k]
                  + f_5 * id_17[k]
                  + pb_y[k] * if__28[k];
    }

#pragma omp simd aligned(t_43, t_44, t_45, pa_x, pa_y, pb_y, gg_s_0, gg_s_44, gg_0, gg_44, \
                         hg_15, hg_44, ig_s_43, ig_s_44, ig_s_45, \
                         if__29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_43[k] = f_2 * ig_s_43[k]
                  + pb_y[k] * if__29[k];

        t_44[k] = -f_8 * gg_s_44[k]
                  + f_9 * gg_44[k]
                  + pa_x[k] * hg_44[k]
                  + f_2 * ig_s_44[k];

        t_45[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_y[k] * hg_15[k]
                  + f_2 * ig_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pb_x, pb_y, pb_z, hf_10, hf_33, id_s_21, ig_s_46, \
                         ig_s_47, ig_s_48, id_21, if__30, if__33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = f_6 * hf_10[k]
                  + f_2 * ig_s_46[k]
                  + pb_y[k] * if__30[k];

        t_47[k] = f_2 * ig_s_47[k]
                  + pb_z[k] * if__30[k];

        t_48[k] = f_9 * hf_33[k]
                  - f_4 * id_s_21[k]
                  + f_2 * ig_s_48[k]
                  + f_5 * id_21[k]
                  + pb_x[k] * if__33[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pb_x, pb_z, hf_36, id_s_18, ig_s_49, ig_s_50, \
                         ig_s_51, id_18, if__31, if__32, if__36 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = f_2 * ig_s_49[k]
                  + pb_z[k] * if__31[k];

        t_50[k] = -f_4 * id_s_18[k]
                  + f_2 * ig_s_50[k]
                  + f_5 * id_18[k]
                  + pb_z[k] * if__32[k];

        t_51[k] = f_9 * hf_36[k]
                  + f_2 * ig_s_51[k]
                  + pb_x[k] * if__36[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_z, hf_38, hf_39, ig_s_52, ig_s_53, \
                         ig_s_54, if__33, if__38, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_2 * ig_s_52[k]
                  + pb_z[k] * if__33[k];

        t_53[k] = f_9 * hf_38[k]
                  + f_2 * ig_s_53[k]
                  + pb_x[k] * if__38[k];

        t_54[k] = f_9 * hf_39[k]
                  + f_2 * ig_s_54[k]
                  + pb_x[k] * if__39[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pa_x, pb_z, gg_s_55, gg_55, hg_55, id_s_21, \
                         ig_s_55, ig_s_56, ig_s_57, id_21, if__36, \
                         if__37 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = -f_12 * gg_s_55[k]
                  + f_3 * gg_55[k]
                  + pa_x[k] * hg_55[k]
                  + f_2 * ig_s_55[k];

        t_56[k] = f_2 * ig_s_56[k]
                  + pb_z[k] * if__36[k];

        t_57[k] = -f_4 * id_s_21[k]
                  + f_2 * ig_s_57[k]
                  + f_5 * id_21[k]
                  + pb_z[k] * if__37[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_y, pb_y, pb_z, hf_19, hg_30, id_s_23, ig_s_58, \
                         ig_s_59, ig_s_60, id_23, if__39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_6 * hf_19[k]
                  + f_2 * ig_s_58[k]
                  + pb_y[k] * if__39[k];

        t_59[k] = -f_1 * id_s_23[k]
                  + f_2 * ig_s_59[k]
                  + f_3 * id_23[k]
                  + pb_z[k] * if__39[k];

        t_60[k] = pa_y[k] * hg_30[k]
                  + f_2 * ig_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pa_y, pa_z, pb_y, hf_22, hg_16, hg_18, hg_32, \
                         ig_s_61, ig_s_62, ig_s_63, ig_s_64, if__42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = pa_z[k] * hg_16[k]
                  + f_2 * ig_s_61[k];

        t_62[k] = pa_y[k] * hg_32[k]
                  + f_2 * ig_s_62[k];

        t_63[k] = pa_z[k] * hg_18[k]
                  + f_2 * ig_s_63[k];

        t_64[k] = f_5 * hf_22[k]
                  + f_2 * ig_s_64[k]
                  + pb_y[k] * if__42[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_y, pa_z, pb_x, hf_47, hg_21, hg_35, ig_s_65, \
                         ig_s_66, ig_s_67, if__47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = pa_y[k] * hg_35[k]
                  + f_2 * ig_s_65[k];

        t_66[k] = pa_z[k] * hg_21[k]
                  + f_2 * ig_s_66[k];

        t_67[k] = f_9 * hf_47[k]
                  + f_2 * ig_s_67[k]
                  + pb_x[k] * if__47[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_y, pa_z, pb_x, hf_48, hg_25, hg_39, ig_s_68, \
                         ig_s_69, ig_s_70, if__48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_9 * hf_48[k]
                  + f_2 * ig_s_68[k]
                  + pb_x[k] * if__48[k];

        t_69[k] = pa_y[k] * hg_39[k]
                  + f_2 * ig_s_69[k];

        t_70[k] = pa_z[k] * hg_25[k]
                  + f_2 * ig_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_x, pb_y, pb_z, gg_s_72, gg_72, hf_16, hf_29, \
                         hg_72, ig_s_71, ig_s_72, ig_s_73, if__46, \
                         if__49 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = f_5 * hf_16[k]
                  + f_2 * ig_s_71[k]
                  + pb_z[k] * if__46[k];

        t_72[k] = -f_12 * gg_s_72[k]
                  + f_3 * gg_72[k]
                  + pa_x[k] * hg_72[k]
                  + f_2 * ig_s_72[k];

        t_73[k] = f_5 * hf_29[k]
                  + f_2 * ig_s_73[k]
                  + pb_y[k] * if__49[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_y, pa_z, pb_y, gg_s_0, gg_0, hg_30, hg_44, \
                         ig_s_74, ig_s_75, ig_s_76, if__50 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = pa_y[k] * hg_44[k]
                  + f_2 * ig_s_74[k];

        t_75[k] = -f_11 * gg_s_0[k]
                  + f_5 * gg_0[k]
                  + pa_z[k] * hg_30[k]
                  + f_2 * ig_s_75[k];

        t_76[k] = f_2 * ig_s_76[k]
                  + pb_y[k] * if__50[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pb_y, pb_z, hf_20, id_s_30, ig_s_77, ig_s_78, \
                         ig_s_79, id_30, if__50, if__51, if__52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_6 * hf_20[k]
                  + f_2 * ig_s_77[k]
                  + pb_z[k] * if__50[k];

        t_78[k] = -f_4 * id_s_30[k]
                  + f_2 * ig_s_78[k]
                  + f_5 * id_30[k]
                  + pb_y[k] * if__51[k];

        t_79[k] = f_2 * ig_s_79[k]
                  + pb_y[k] * if__52[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pb_x, hf_55, hf_56, hf_57, id_s_35, ig_s_80, \
                         ig_s_81, ig_s_82, id_35, if__55, if__56, \
                         if__57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = f_9 * hf_55[k]
                  - f_4 * id_s_35[k]
                  + f_2 * ig_s_80[k]
                  + f_5 * id_35[k]
                  + pb_x[k] * if__55[k];

        t_81[k] = f_9 * hf_56[k]
                  + f_2 * ig_s_81[k]
                  + pb_x[k] * if__56[k];

        t_82[k] = f_9 * hf_57[k]
                  + f_2 * ig_s_82[k]
                  + pb_x[k] * if__57[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pb_x, pb_y, hf_59, id_s_33, ig_s_83, ig_s_84, \
                         ig_s_85, id_33, if__55, if__56, if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_2 * ig_s_83[k]
                  + pb_y[k] * if__55[k];

        t_84[k] = f_9 * hf_59[k]
                  + f_2 * ig_s_84[k]
                  + pb_x[k] * if__59[k];

        t_85[k] = -f_1 * id_s_33[k]
                  + f_2 * ig_s_85[k]
                  + f_3 * id_33[k]
                  + pb_y[k] * if__56[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pb_y, id_s_34, id_s_35, ig_s_86, ig_s_87, ig_s_88, \
                         id_34, id_35, if__57, if__58, if__59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_10 * id_s_34[k]
                  + f_2 * ig_s_86[k]
                  + f_6 * id_34[k]
                  + pb_y[k] * if__57[k];

        t_87[k] = -f_4 * id_s_35[k]
                  + f_2 * ig_s_87[k]
                  + f_5 * id_35[k]
                  + pb_y[k] * if__58[k];

        t_88[k] = f_2 * ig_s_88[k]
                  + pb_y[k] * if__59[k];
    }

#pragma omp simd aligned(t_89, t_90, pa_x, pa_y, gg_s_15, gg_s_89, gg_15, gg_89, hg_45, hg_89, \
                         ig_s_89, ig_s_90 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = -f_12 * gg_s_89[k]
                  + f_3 * gg_89[k]
                  + pa_x[k] * hg_89[k]
                  + f_2 * ig_s_89[k];

        t_90[k] = -f_13 * gg_s_15[k]
                  + f_6 * gg_15[k]
                  + pa_y[k] * hg_45[k]
                  + f_2 * ig_s_90[k];
    }

#pragma omp simd aligned(t_91, t_92, t_93, pb_x, pb_y, pb_z, hf_30, hf_63, id_s_39, ig_s_91, \
                         ig_s_92, ig_s_93, id_39, if__60, if__63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_91[k] = f_3 * hf_30[k]
                  + f_2 * ig_s_91[k]
                  + pb_y[k] * if__60[k];

        t_92[k] = f_2 * ig_s_92[k]
                  + pb_z[k] * if__60[k];

        t_93[k] = f_3 * hf_63[k]
                  - f_4 * id_s_39[k]
                  + f_2 * ig_s_93[k]
                  + f_5 * id_39[k]
                  + pb_x[k] * if__63[k];
    }

#pragma omp simd aligned(t_94, t_95, t_96, pb_x, pb_z, hf_66, id_s_36, ig_s_94, ig_s_95, \
                         ig_s_96, id_36, if__61, if__62, if__66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_94[k] = f_2 * ig_s_94[k]
                  + pb_z[k] * if__61[k];

        t_95[k] = -f_4 * id_s_36[k]
                  + f_2 * ig_s_95[k]
                  + f_5 * id_36[k]
                  + pb_z[k] * if__62[k];

        t_96[k] = f_3 * hf_66[k]
                  + f_2 * ig_s_96[k]
                  + pb_x[k] * if__66[k];
    }

#pragma omp simd aligned(t_97, t_98, t_99, pb_x, pb_z, hf_68, hf_69, ig_s_97, ig_s_98, \
                         ig_s_99, if__63, if__68, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_97[k] = f_2 * ig_s_97[k]
                  + pb_z[k] * if__63[k];

        t_98[k] = f_3 * hf_68[k]
                  + f_2 * ig_s_98[k]
                  + pb_x[k] * if__68[k];

        t_99[k] = f_3 * hf_69[k]
                  + f_2 * ig_s_99[k]
                  + pb_x[k] * if__69[k];
    }

#pragma omp simd aligned(t_100, t_101, t_102, pa_x, pb_z, gg_s_100, gg_100, hg_100, id_s_39, \
                         ig_s_100, ig_s_101, ig_s_102, id_39, if__66, \
                         if__67 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_100[k] = -f_13 * gg_s_100[k]
                   + f_6 * gg_100[k]
                   + pa_x[k] * hg_100[k]
                   + f_2 * ig_s_100[k];

        t_101[k] = f_2 * ig_s_101[k]
                   + pb_z[k] * if__66[k];

        t_102[k] = -f_4 * id_s_39[k]
                   + f_2 * ig_s_102[k]
                   + f_5 * id_39[k]
                   + pb_z[k] * if__67[k];
    }
}

static auto
compute_prim_ig_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gg_s,
                                        const size_t gg, const size_t hf, const size_t hg,
                                        const size_t id_s, const size_t ig_s, const size_t id,
                                        const size_t if_, const size_t ncols, const double alpha,
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
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;
    const auto f_13 = 2.0 * beta / p;

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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_s_30 = buffer.data(gg_s + 30);
    const auto *gg_s_35 = buffer.data(gg_s + 35);
    const auto *gg_s_45 = buffer.data(gg_s + 45);
    const auto *gg_s_48 = buffer.data(gg_s + 48);
    const auto *gg_s_65 = buffer.data(gg_s + 65);
    const auto *gg_s_75 = buffer.data(gg_s + 75);
    const auto *gg_s_80 = buffer.data(gg_s + 80);
    const auto *gg_s_117 = buffer.data(gg_s + 117);
    const auto *gg_s_119 = buffer.data(gg_s + 119);
    const auto *gg_s_130 = buffer.data(gg_s + 130);
    const auto *gg_s_132 = buffer.data(gg_s + 132);
    const auto *gg_s_149 = buffer.data(gg_s + 149);
    const auto *gg_s_160 = buffer.data(gg_s + 160);
    const auto *gg_s_177 = buffer.data(gg_s + 177);
    const auto *gg_s_179 = buffer.data(gg_s + 179);
    const auto *gg_s_190 = buffer.data(gg_s + 190);
    const auto *gg_s_192 = buffer.data(gg_s + 192);
    const auto *gg_s_194 = buffer.data(gg_s + 194);

    const auto *gg_30 = buffer.data(gg + 30);
    const auto *gg_35 = buffer.data(gg + 35);
    const auto *gg_45 = buffer.data(gg + 45);
    const auto *gg_48 = buffer.data(gg + 48);
    const auto *gg_65 = buffer.data(gg + 65);
    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_80 = buffer.data(gg + 80);
    const auto *gg_117 = buffer.data(gg + 117);
    const auto *gg_119 = buffer.data(gg + 119);
    const auto *gg_130 = buffer.data(gg + 130);
    const auto *gg_132 = buffer.data(gg + 132);
    const auto *gg_149 = buffer.data(gg + 149);
    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_177 = buffer.data(gg + 177);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_192 = buffer.data(gg + 192);
    const auto *gg_194 = buffer.data(gg + 194);

    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_39 = buffer.data(hf + 39);
    const auto *hf_42 = buffer.data(hf + 42);
    const auto *hf_46 = buffer.data(hf + 46);
    const auto *hf_49 = buffer.data(hf + 49);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_51 = buffer.data(hf + 51);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_66 = buffer.data(hf + 66);
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
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_108 = buffer.data(hf + 108);
    const auto *hf_109 = buffer.data(hf + 109);
    const auto *hf_117 = buffer.data(hf + 117);
    const auto *hf_118 = buffer.data(hf + 118);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_127 = buffer.data(hf + 127);
    const auto *hf_128 = buffer.data(hf + 128);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_137 = buffer.data(hf + 137);
    const auto *hf_138 = buffer.data(hf + 138);

    const auto *hg_45 = buffer.data(hg + 45);
    const auto *hg_46 = buffer.data(hg + 46);
    const auto *hg_48 = buffer.data(hg + 48);
    const auto *hg_51 = buffer.data(hg + 51);
    const auto *hg_55 = buffer.data(hg + 55);
    const auto *hg_65 = buffer.data(hg + 65);
    const auto *hg_75 = buffer.data(hg + 75);
    const auto *hg_77 = buffer.data(hg + 77);
    const auto *hg_78 = buffer.data(hg + 78);
    const auto *hg_80 = buffer.data(hg + 80);
    const auto *hg_84 = buffer.data(hg + 84);
    const auto *hg_89 = buffer.data(hg + 89);
    const auto *hg_90 = buffer.data(hg + 90);
    const auto *hg_91 = buffer.data(hg + 91);
    const auto *hg_93 = buffer.data(hg + 93);
    const auto *hg_96 = buffer.data(hg + 96);
    const auto *hg_100 = buffer.data(hg + 100);
    const auto *hg_108 = buffer.data(hg + 108);
    const auto *hg_110 = buffer.data(hg + 110);
    const auto *hg_117 = buffer.data(hg + 117);
    const auto *hg_119 = buffer.data(hg + 119);
    const auto *hg_120 = buffer.data(hg + 120);
    const auto *hg_125 = buffer.data(hg + 125);
    const auto *hg_130 = buffer.data(hg + 130);
    const auto *hg_132 = buffer.data(hg + 132);
    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_137 = buffer.data(hg + 137);
    const auto *hg_138 = buffer.data(hg + 138);
    const auto *hg_140 = buffer.data(hg + 140);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_160 = buffer.data(hg + 160);
    const auto *hg_177 = buffer.data(hg + 177);
    const auto *hg_179 = buffer.data(hg + 179);
    const auto *hg_190 = buffer.data(hg + 190);
    const auto *hg_192 = buffer.data(hg + 192);
    const auto *hg_194 = buffer.data(hg + 194);

    const auto *id_s_41 = buffer.data(id_s + 41);
    const auto *id_s_54 = buffer.data(id_s + 54);
    const auto *id_s_57 = buffer.data(id_s + 57);
    const auto *id_s_58 = buffer.data(id_s + 58);
    const auto *id_s_59 = buffer.data(id_s + 59);
    const auto *id_s_60 = buffer.data(id_s + 60);
    const auto *id_s_63 = buffer.data(id_s + 63);
    const auto *id_s_65 = buffer.data(id_s + 65);

    const auto *ig_s_103 = buffer.data(ig_s + 103);
    const auto *ig_s_104 = buffer.data(ig_s + 104);
    const auto *ig_s_105 = buffer.data(ig_s + 105);
    const auto *ig_s_106 = buffer.data(ig_s + 106);
    const auto *ig_s_107 = buffer.data(ig_s + 107);
    const auto *ig_s_108 = buffer.data(ig_s + 108);
    const auto *ig_s_109 = buffer.data(ig_s + 109);
    const auto *ig_s_110 = buffer.data(ig_s + 110);
    const auto *ig_s_111 = buffer.data(ig_s + 111);
    const auto *ig_s_112 = buffer.data(ig_s + 112);
    const auto *ig_s_113 = buffer.data(ig_s + 113);
    const auto *ig_s_114 = buffer.data(ig_s + 114);
    const auto *ig_s_115 = buffer.data(ig_s + 115);
    const auto *ig_s_116 = buffer.data(ig_s + 116);
    const auto *ig_s_117 = buffer.data(ig_s + 117);
    const auto *ig_s_118 = buffer.data(ig_s + 118);
    const auto *ig_s_119 = buffer.data(ig_s + 119);
    const auto *ig_s_120 = buffer.data(ig_s + 120);
    const auto *ig_s_121 = buffer.data(ig_s + 121);
    const auto *ig_s_122 = buffer.data(ig_s + 122);
    const auto *ig_s_123 = buffer.data(ig_s + 123);
    const auto *ig_s_124 = buffer.data(ig_s + 124);
    const auto *ig_s_125 = buffer.data(ig_s + 125);
    const auto *ig_s_126 = buffer.data(ig_s + 126);
    const auto *ig_s_127 = buffer.data(ig_s + 127);
    const auto *ig_s_128 = buffer.data(ig_s + 128);
    const auto *ig_s_129 = buffer.data(ig_s + 129);
    const auto *ig_s_130 = buffer.data(ig_s + 130);
    const auto *ig_s_131 = buffer.data(ig_s + 131);
    const auto *ig_s_132 = buffer.data(ig_s + 132);
    const auto *ig_s_133 = buffer.data(ig_s + 133);
    const auto *ig_s_134 = buffer.data(ig_s + 134);
    const auto *ig_s_135 = buffer.data(ig_s + 135);
    const auto *ig_s_136 = buffer.data(ig_s + 136);
    const auto *ig_s_137 = buffer.data(ig_s + 137);
    const auto *ig_s_138 = buffer.data(ig_s + 138);
    const auto *ig_s_139 = buffer.data(ig_s + 139);
    const auto *ig_s_140 = buffer.data(ig_s + 140);
    const auto *ig_s_141 = buffer.data(ig_s + 141);
    const auto *ig_s_142 = buffer.data(ig_s + 142);
    const auto *ig_s_143 = buffer.data(ig_s + 143);
    const auto *ig_s_144 = buffer.data(ig_s + 144);
    const auto *ig_s_145 = buffer.data(ig_s + 145);
    const auto *ig_s_146 = buffer.data(ig_s + 146);
    const auto *ig_s_147 = buffer.data(ig_s + 147);
    const auto *ig_s_148 = buffer.data(ig_s + 148);
    const auto *ig_s_149 = buffer.data(ig_s + 149);
    const auto *ig_s_150 = buffer.data(ig_s + 150);
    const auto *ig_s_151 = buffer.data(ig_s + 151);
    const auto *ig_s_152 = buffer.data(ig_s + 152);
    const auto *ig_s_153 = buffer.data(ig_s + 153);
    const auto *ig_s_154 = buffer.data(ig_s + 154);
    const auto *ig_s_155 = buffer.data(ig_s + 155);
    const auto *ig_s_156 = buffer.data(ig_s + 156);
    const auto *ig_s_157 = buffer.data(ig_s + 157);
    const auto *ig_s_158 = buffer.data(ig_s + 158);
    const auto *ig_s_159 = buffer.data(ig_s + 159);
    const auto *ig_s_160 = buffer.data(ig_s + 160);
    const auto *ig_s_161 = buffer.data(ig_s + 161);
    const auto *ig_s_162 = buffer.data(ig_s + 162);
    const auto *ig_s_163 = buffer.data(ig_s + 163);
    const auto *ig_s_164 = buffer.data(ig_s + 164);
    const auto *ig_s_165 = buffer.data(ig_s + 165);
    const auto *ig_s_166 = buffer.data(ig_s + 166);
    const auto *ig_s_167 = buffer.data(ig_s + 167);
    const auto *ig_s_168 = buffer.data(ig_s + 168);
    const auto *ig_s_169 = buffer.data(ig_s + 169);
    const auto *ig_s_170 = buffer.data(ig_s + 170);
    const auto *ig_s_171 = buffer.data(ig_s + 171);
    const auto *ig_s_172 = buffer.data(ig_s + 172);
    const auto *ig_s_173 = buffer.data(ig_s + 173);
    const auto *ig_s_174 = buffer.data(ig_s + 174);
    const auto *ig_s_175 = buffer.data(ig_s + 175);
    const auto *ig_s_176 = buffer.data(ig_s + 176);
    const auto *ig_s_177 = buffer.data(ig_s + 177);
    const auto *ig_s_178 = buffer.data(ig_s + 178);
    const auto *ig_s_179 = buffer.data(ig_s + 179);
    const auto *ig_s_180 = buffer.data(ig_s + 180);
    const auto *ig_s_181 = buffer.data(ig_s + 181);
    const auto *ig_s_182 = buffer.data(ig_s + 182);
    const auto *ig_s_183 = buffer.data(ig_s + 183);
    const auto *ig_s_184 = buffer.data(ig_s + 184);
    const auto *ig_s_185 = buffer.data(ig_s + 185);
    const auto *ig_s_186 = buffer.data(ig_s + 186);
    const auto *ig_s_187 = buffer.data(ig_s + 187);
    const auto *ig_s_188 = buffer.data(ig_s + 188);
    const auto *ig_s_189 = buffer.data(ig_s + 189);
    const auto *ig_s_190 = buffer.data(ig_s + 190);
    const auto *ig_s_191 = buffer.data(ig_s + 191);
    const auto *ig_s_192 = buffer.data(ig_s + 192);
    const auto *ig_s_193 = buffer.data(ig_s + 193);
    const auto *ig_s_194 = buffer.data(ig_s + 194);
    const auto *ig_s_195 = buffer.data(ig_s + 195);
    const auto *ig_s_196 = buffer.data(ig_s + 196);
    const auto *ig_s_197 = buffer.data(ig_s + 197);
    const auto *ig_s_198 = buffer.data(ig_s + 198);
    const auto *ig_s_199 = buffer.data(ig_s + 199);
    const auto *ig_s_200 = buffer.data(ig_s + 200);
    const auto *ig_s_201 = buffer.data(ig_s + 201);
    const auto *ig_s_202 = buffer.data(ig_s + 202);
    const auto *ig_s_203 = buffer.data(ig_s + 203);

    const auto *id_41 = buffer.data(id + 41);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);

    const auto *if__69 = buffer.data(if_ + 69);
    const auto *if__70 = buffer.data(if_ + 70);
    const auto *if__72 = buffer.data(if_ + 72);
    const auto *if__76 = buffer.data(if_ + 76);
    const auto *if__77 = buffer.data(if_ + 77);
    const auto *if__78 = buffer.data(if_ + 78);
    const auto *if__79 = buffer.data(if_ + 79);
    const auto *if__80 = buffer.data(if_ + 80);
    const auto *if__82 = buffer.data(if_ + 82);
    const auto *if__86 = buffer.data(if_ + 86);
    const auto *if__87 = buffer.data(if_ + 87);
    const auto *if__88 = buffer.data(if_ + 88);
    const auto *if__89 = buffer.data(if_ + 89);
    const auto *if__90 = buffer.data(if_ + 90);
    const auto *if__91 = buffer.data(if_ + 91);
    const auto *if__92 = buffer.data(if_ + 92);
    const auto *if__95 = buffer.data(if_ + 95);
    const auto *if__96 = buffer.data(if_ + 96);
    const auto *if__97 = buffer.data(if_ + 97);
    const auto *if__98 = buffer.data(if_ + 98);
    const auto *if__99 = buffer.data(if_ + 99);
    const auto *if__100 = buffer.data(if_ + 100);
    const auto *if__101 = buffer.data(if_ + 101);
    const auto *if__102 = buffer.data(if_ + 102);
    const auto *if__103 = buffer.data(if_ + 103);
    const auto *if__106 = buffer.data(if_ + 106);
    const auto *if__107 = buffer.data(if_ + 107);
    const auto *if__108 = buffer.data(if_ + 108);
    const auto *if__109 = buffer.data(if_ + 109);
    const auto *if__110 = buffer.data(if_ + 110);
    const auto *if__112 = buffer.data(if_ + 112);
    const auto *if__116 = buffer.data(if_ + 116);
    const auto *if__117 = buffer.data(if_ + 117);
    const auto *if__118 = buffer.data(if_ + 118);
    const auto *if__119 = buffer.data(if_ + 119);
    const auto *if__120 = buffer.data(if_ + 120);
    const auto *if__122 = buffer.data(if_ + 122);
    const auto *if__126 = buffer.data(if_ + 126);
    const auto *if__127 = buffer.data(if_ + 127);
    const auto *if__128 = buffer.data(if_ + 128);
    const auto *if__129 = buffer.data(if_ + 129);
    const auto *if__130 = buffer.data(if_ + 130);
    const auto *if__132 = buffer.data(if_ + 132);
    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__137 = buffer.data(if_ + 137);
    const auto *if__138 = buffer.data(if_ + 138);

#pragma omp simd aligned(t_103, t_104, t_105, pa_z, pb_y, pb_z, hf_39, hg_45, id_s_41, \
                         ig_s_103, ig_s_104, ig_s_105, id_41, if__69 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_103[k] = f_3 * hf_39[k]
                   + f_2 * ig_s_103[k]
                   + pb_y[k] * if__69[k];

        t_104[k] = -f_1 * id_s_41[k]
                   + f_2 * ig_s_104[k]
                   + f_3 * id_41[k]
                   + pb_z[k] * if__69[k];

        t_105[k] = pa_z[k] * hg_45[k]
                   + f_2 * ig_s_105[k];
    }

#pragma omp simd aligned(t_106, t_107, t_108, pa_z, pb_z, hf_30, hg_46, hg_48, ig_s_106, \
                         ig_s_107, ig_s_108, if__70 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_106[k] = pa_z[k] * hg_46[k]
                   + f_2 * ig_s_106[k];

        t_107[k] = f_5 * hf_30[k]
                   + f_2 * ig_s_107[k]
                   + pb_z[k] * if__70[k];

        t_108[k] = pa_z[k] * hg_48[k]
                   + f_2 * ig_s_108[k];
    }

#pragma omp simd aligned(t_109, t_110, t_111, pa_y, pa_z, pb_y, gg_s_35, gg_35, hf_42, hg_51, \
                         hg_65, ig_s_109, ig_s_110, ig_s_111, if__72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_109[k] = f_6 * hf_42[k]
                   + f_2 * ig_s_109[k]
                   + pb_y[k] * if__72[k];

        t_110[k] = -f_11 * gg_s_35[k]
                   + f_5 * gg_35[k]
                   + pa_y[k] * hg_65[k]
                   + f_2 * ig_s_110[k];

        t_111[k] = pa_z[k] * hg_51[k]
                   + f_2 * ig_s_111[k];
    }

#pragma omp simd aligned(t_112, t_113, t_114, pb_x, hf_77, hf_78, hf_79, ig_s_112, ig_s_113, \
                         ig_s_114, if__77, if__78, if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_112[k] = f_3 * hf_77[k]
                   + f_2 * ig_s_112[k]
                   + pb_x[k] * if__77[k];

        t_113[k] = f_3 * hf_78[k]
                   + f_2 * ig_s_113[k]
                   + pb_x[k] * if__78[k];

        t_114[k] = f_3 * hf_79[k]
                   + f_2 * ig_s_114[k]
                   + pb_x[k] * if__79[k];
    }

#pragma omp simd aligned(t_115, t_116, t_117, pa_x, pa_z, pb_z, gg_s_117, gg_117, hf_36, \
                         hg_55, hg_117, ig_s_115, ig_s_116, ig_s_117, \
                         if__76 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_115[k] = pa_z[k] * hg_55[k]
                   + f_2 * ig_s_115[k];

        t_116[k] = f_5 * hf_36[k]
                   + f_2 * ig_s_116[k]
                   + pb_z[k] * if__76[k];

        t_117[k] = -f_13 * gg_s_117[k]
                   + f_6 * gg_117[k]
                   + pa_x[k] * hg_117[k]
                   + f_2 * ig_s_117[k];
    }

#pragma omp simd aligned(t_118, t_119, t_120, pa_x, pa_y, pb_y, gg_s_119, gg_119, hf_49, \
                         hg_75, hg_119, ig_s_118, ig_s_119, ig_s_120, \
                         if__79 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_118[k] = f_6 * hf_49[k]
                   + f_2 * ig_s_118[k]
                   + pb_y[k] * if__79[k];

        t_119[k] = -f_13 * gg_s_119[k]
                   + f_6 * gg_119[k]
                   + pa_x[k] * hg_119[k]
                   + f_2 * ig_s_119[k];

        t_120[k] = pa_y[k] * hg_75[k]
                   + f_2 * ig_s_120[k];
    }

#pragma omp simd aligned(t_121, t_122, t_123, pa_y, pb_y, hf_50, hf_51, hg_77, hg_78, \
                         ig_s_121, ig_s_122, ig_s_123, if__80 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_121[k] = f_5 * hf_50[k]
                   + f_2 * ig_s_121[k]
                   + pb_y[k] * if__80[k];

        t_122[k] = pa_y[k] * hg_77[k]
                   + f_2 * ig_s_122[k];

        t_123[k] = f_6 * hf_51[k]
                   + pa_y[k] * hg_78[k]
                   + f_2 * ig_s_123[k];
    }

#pragma omp simd aligned(t_124, t_125, t_126, pa_y, pb_x, pb_y, hf_52, hf_86, hg_80, ig_s_124, \
                         ig_s_125, ig_s_126, if__82, if__86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_124[k] = f_5 * hf_52[k]
                   + f_2 * ig_s_124[k]
                   + pb_y[k] * if__82[k];

        t_125[k] = pa_y[k] * hg_80[k]
                   + f_2 * ig_s_125[k];

        t_126[k] = f_3 * hf_86[k]
                   + f_2 * ig_s_126[k]
                   + pb_x[k] * if__86[k];
    }

#pragma omp simd aligned(t_127, t_128, t_129, pa_y, pb_x, hf_87, hf_88, hg_84, ig_s_127, \
                         ig_s_128, ig_s_129, if__87, if__88 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_127[k] = f_3 * hf_87[k]
                   + f_2 * ig_s_127[k]
                   + pb_x[k] * if__87[k];

        t_128[k] = f_3 * hf_88[k]
                   + f_2 * ig_s_128[k]
                   + pb_x[k] * if__88[k];

        t_129[k] = pa_y[k] * hg_84[k]
                   + f_2 * ig_s_129[k];
    }

#pragma omp simd aligned(t_130, t_131, t_132, pa_x, pb_z, gg_s_130, gg_s_132, gg_130, gg_132, \
                         hf_46, hg_130, hg_132, ig_s_130, ig_s_131, ig_s_132, \
                         if__86 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_130[k] = -f_13 * gg_s_130[k]
                   + f_6 * gg_130[k]
                   + pa_x[k] * hg_130[k]
                   + f_2 * ig_s_130[k];

        t_131[k] = f_6 * hf_46[k]
                   + f_2 * ig_s_131[k]
                   + pb_z[k] * if__86[k];

        t_132[k] = -f_13 * gg_s_132[k]
                   + f_6 * gg_132[k]
                   + pa_x[k] * hg_132[k]
                   + f_2 * ig_s_132[k];
    }

#pragma omp simd aligned(t_133, t_134, t_135, pa_y, pa_z, pb_y, gg_s_30, gg_30, hf_59, hg_75, \
                         hg_89, ig_s_133, ig_s_134, ig_s_135, if__89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_133[k] = f_5 * hf_59[k]
                   + f_2 * ig_s_133[k]
                   + pb_y[k] * if__89[k];

        t_134[k] = pa_y[k] * hg_89[k]
                   + f_2 * ig_s_134[k];

        t_135[k] = -f_13 * gg_s_30[k]
                   + f_6 * gg_30[k]
                   + pa_z[k] * hg_75[k]
                   + f_2 * ig_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, t_139, pb_y, pb_z, hf_50, id_s_54, ig_s_136, \
                         ig_s_137, ig_s_138, ig_s_139, id_54, if__90, if__91, \
                         if__92 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = f_2 * ig_s_136[k]
                   + pb_y[k] * if__90[k];

        t_137[k] = f_3 * hf_50[k]
                   + f_2 * ig_s_137[k]
                   + pb_z[k] * if__90[k];

        t_138[k] = -f_4 * id_s_54[k]
                   + f_2 * ig_s_138[k]
                   + f_5 * id_54[k]
                   + pb_y[k] * if__91[k];

        t_139[k] = f_2 * ig_s_139[k]
                   + pb_y[k] * if__92[k];
    }

#pragma omp simd aligned(t_140, t_141, t_142, pb_x, hf_95, hf_96, hf_97, id_s_59, ig_s_140, \
                         ig_s_141, ig_s_142, id_59, if__95, if__96, \
                         if__97 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_140[k] = f_3 * hf_95[k]
                   - f_4 * id_s_59[k]
                   + f_2 * ig_s_140[k]
                   + f_5 * id_59[k]
                   + pb_x[k] * if__95[k];

        t_141[k] = f_3 * hf_96[k]
                   + f_2 * ig_s_141[k]
                   + pb_x[k] * if__96[k];

        t_142[k] = f_3 * hf_97[k]
                   + f_2 * ig_s_142[k]
                   + pb_x[k] * if__97[k];
    }

#pragma omp simd aligned(t_143, t_144, t_145, pb_x, pb_y, hf_99, id_s_57, ig_s_143, ig_s_144, \
                         ig_s_145, id_57, if__95, if__96, if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_143[k] = f_2 * ig_s_143[k]
                   + pb_y[k] * if__95[k];

        t_144[k] = f_3 * hf_99[k]
                   + f_2 * ig_s_144[k]
                   + pb_x[k] * if__99[k];

        t_145[k] = -f_1 * id_s_57[k]
                   + f_2 * ig_s_145[k]
                   + f_3 * id_57[k]
                   + pb_y[k] * if__96[k];
    }

#pragma omp simd aligned(t_146, t_147, t_148, pb_y, id_s_58, id_s_59, ig_s_146, ig_s_147, \
                         ig_s_148, id_58, id_59, if__97, if__98, \
                         if__99 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_146[k] = -f_10 * id_s_58[k]
                   + f_2 * ig_s_146[k]
                   + f_6 * id_58[k]
                   + pb_y[k] * if__97[k];

        t_147[k] = -f_4 * id_s_59[k]
                   + f_2 * ig_s_147[k]
                   + f_5 * id_59[k]
                   + pb_y[k] * if__98[k];

        t_148[k] = f_2 * ig_s_148[k]
                   + pb_y[k] * if__99[k];
    }

#pragma omp simd aligned(t_149, t_150, pa_x, pa_y, gg_s_45, gg_s_149, gg_45, gg_149, hg_90, \
                         hg_149, ig_s_149, ig_s_150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_149[k] = -f_13 * gg_s_149[k]
                   + f_6 * gg_149[k]
                   + pa_x[k] * hg_149[k]
                   + f_2 * ig_s_149[k];

        t_150[k] = -f_12 * gg_s_45[k]
                   + f_3 * gg_45[k]
                   + pa_y[k] * hg_90[k]
                   + f_2 * ig_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, pb_x, pb_y, pb_z, hf_60, hf_103, id_s_63, \
                         ig_s_151, ig_s_152, ig_s_153, id_63, if__100, \
                         if__103 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_9 * hf_60[k]
                   + f_2 * ig_s_151[k]
                   + pb_y[k] * if__100[k];

        t_152[k] = f_2 * ig_s_152[k]
                   + pb_z[k] * if__100[k];

        t_153[k] = f_6 * hf_103[k]
                   - f_4 * id_s_63[k]
                   + f_2 * ig_s_153[k]
                   + f_5 * id_63[k]
                   + pb_x[k] * if__103[k];
    }

#pragma omp simd aligned(t_154, t_155, t_156, pb_x, pb_z, hf_106, id_s_60, ig_s_154, ig_s_155, \
                         ig_s_156, id_60, if__101, if__102, if__106 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_154[k] = f_2 * ig_s_154[k]
                   + pb_z[k] * if__101[k];

        t_155[k] = -f_4 * id_s_60[k]
                   + f_2 * ig_s_155[k]
                   + f_5 * id_60[k]
                   + pb_z[k] * if__102[k];

        t_156[k] = f_6 * hf_106[k]
                   + f_2 * ig_s_156[k]
                   + pb_x[k] * if__106[k];
    }

#pragma omp simd aligned(t_157, t_158, t_159, pb_x, pb_z, hf_108, hf_109, ig_s_157, ig_s_158, \
                         ig_s_159, if__103, if__108, if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_157[k] = f_2 * ig_s_157[k]
                   + pb_z[k] * if__103[k];

        t_158[k] = f_6 * hf_108[k]
                   + f_2 * ig_s_158[k]
                   + pb_x[k] * if__108[k];

        t_159[k] = f_6 * hf_109[k]
                   + f_2 * ig_s_159[k]
                   + pb_x[k] * if__109[k];
    }

#pragma omp simd aligned(t_160, t_161, t_162, pa_x, pb_z, gg_s_160, gg_160, hg_160, id_s_63, \
                         ig_s_160, ig_s_161, ig_s_162, id_63, if__106, \
                         if__107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_160[k] = -f_11 * gg_s_160[k]
                   + f_5 * gg_160[k]
                   + pa_x[k] * hg_160[k]
                   + f_2 * ig_s_160[k];

        t_161[k] = f_2 * ig_s_161[k]
                   + pb_z[k] * if__106[k];

        t_162[k] = -f_4 * id_s_63[k]
                   + f_2 * ig_s_162[k]
                   + f_5 * id_63[k]
                   + pb_z[k] * if__107[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_z, pb_y, pb_z, hf_69, hg_90, id_s_65, \
                         ig_s_163, ig_s_164, ig_s_165, id_65, if__109 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = f_9 * hf_69[k]
                   + f_2 * ig_s_163[k]
                   + pb_y[k] * if__109[k];

        t_164[k] = -f_1 * id_s_65[k]
                   + f_2 * ig_s_164[k]
                   + f_3 * id_65[k]
                   + pb_z[k] * if__109[k];

        t_165[k] = pa_z[k] * hg_90[k]
                   + f_2 * ig_s_165[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, pa_z, pb_z, hf_60, hg_91, hg_93, ig_s_166, \
                         ig_s_167, ig_s_168, if__110 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_z[k] * hg_91[k]
                   + f_2 * ig_s_166[k];

        t_167[k] = f_5 * hf_60[k]
                   + f_2 * ig_s_167[k]
                   + pb_z[k] * if__110[k];

        t_168[k] = pa_z[k] * hg_93[k]
                   + f_2 * ig_s_168[k];
    }

#pragma omp simd aligned(t_169, t_170, t_171, pa_y, pa_z, pb_y, gg_s_65, gg_65, hf_72, hg_96, \
                         hg_110, ig_s_169, ig_s_170, ig_s_171, \
                         if__112 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_169[k] = f_3 * hf_72[k]
                   + f_2 * ig_s_169[k]
                   + pb_y[k] * if__112[k];

        t_170[k] = -f_13 * gg_s_65[k]
                   + f_6 * gg_65[k]
                   + pa_y[k] * hg_110[k]
                   + f_2 * ig_s_170[k];

        t_171[k] = pa_z[k] * hg_96[k]
                   + f_2 * ig_s_171[k];
    }

#pragma omp simd aligned(t_172, t_173, t_174, pb_x, hf_117, hf_118, hf_119, ig_s_172, \
                         ig_s_173, ig_s_174, if__117, if__118, \
                         if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_172[k] = f_6 * hf_117[k]
                   + f_2 * ig_s_172[k]
                   + pb_x[k] * if__117[k];

        t_173[k] = f_6 * hf_118[k]
                   + f_2 * ig_s_173[k]
                   + pb_x[k] * if__118[k];

        t_174[k] = f_6 * hf_119[k]
                   + f_2 * ig_s_174[k]
                   + pb_x[k] * if__119[k];
    }

#pragma omp simd aligned(t_175, t_176, t_177, pa_x, pa_z, pb_z, gg_s_177, gg_177, hf_66, \
                         hg_100, hg_177, ig_s_175, ig_s_176, ig_s_177, \
                         if__116 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_175[k] = pa_z[k] * hg_100[k]
                   + f_2 * ig_s_175[k];

        t_176[k] = f_5 * hf_66[k]
                   + f_2 * ig_s_176[k]
                   + pb_z[k] * if__116[k];

        t_177[k] = -f_11 * gg_s_177[k]
                   + f_5 * gg_177[k]
                   + pa_x[k] * hg_177[k]
                   + f_2 * ig_s_177[k];
    }

#pragma omp simd aligned(t_178, t_179, pa_x, pb_y, gg_s_179, gg_179, hf_79, hg_179, ig_s_178, \
                         ig_s_179, if__119 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_178[k] = f_3 * hf_79[k]
                   + f_2 * ig_s_178[k]
                   + pb_y[k] * if__119[k];

        t_179[k] = -f_11 * gg_s_179[k]
                   + f_5 * gg_179[k]
                   + pa_x[k] * hg_179[k]
                   + f_2 * ig_s_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_y, pb_y, pb_z, gg_s_75, gg_75, hf_70, hf_80, \
                         hg_120, ig_s_180, ig_s_181, ig_s_182, \
                         if__120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = -f_11 * gg_s_75[k]
                   + f_5 * gg_75[k]
                   + pa_y[k] * hg_120[k]
                   + f_2 * ig_s_180[k];

        t_181[k] = f_6 * hf_80[k]
                   + f_2 * ig_s_181[k]
                   + pb_y[k] * if__120[k];

        t_182[k] = f_6 * hf_70[k]
                   + f_2 * ig_s_182[k]
                   + pb_z[k] * if__120[k];
    }

#pragma omp simd aligned(t_183, t_184, pa_z, pb_y, gg_s_48, gg_48, hf_82, hg_108, ig_s_183, \
                         ig_s_184, if__122 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = -f_11 * gg_s_48[k]
                   + f_5 * gg_48[k]
                   + pa_z[k] * hg_108[k]
                   + f_2 * ig_s_183[k];

        t_184[k] = f_6 * hf_82[k]
                   + f_2 * ig_s_184[k]
                   + pb_y[k] * if__122[k];
    }

#pragma omp simd aligned(t_185, t_186, t_187, pa_y, pb_x, gg_s_80, gg_80, hf_126, hf_127, \
                         hg_125, ig_s_185, ig_s_186, ig_s_187, if__126, \
                         if__127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_185[k] = -f_11 * gg_s_80[k]
                   + f_5 * gg_80[k]
                   + pa_y[k] * hg_125[k]
                   + f_2 * ig_s_185[k];

        t_186[k] = f_6 * hf_126[k]
                   + f_2 * ig_s_186[k]
                   + pb_x[k] * if__126[k];

        t_187[k] = f_6 * hf_127[k]
                   + f_2 * ig_s_187[k]
                   + pb_x[k] * if__127[k];
    }

#pragma omp simd aligned(t_188, t_189, t_190, pa_x, pb_x, gg_s_190, gg_190, hf_128, hf_129, \
                         hg_190, ig_s_188, ig_s_189, ig_s_190, if__128, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_188[k] = f_6 * hf_128[k]
                   + f_2 * ig_s_188[k]
                   + pb_x[k] * if__128[k];

        t_189[k] = f_6 * hf_129[k]
                   + f_2 * ig_s_189[k]
                   + pb_x[k] * if__129[k];

        t_190[k] = -f_11 * gg_s_190[k]
                   + f_5 * gg_190[k]
                   + pa_x[k] * hg_190[k]
                   + f_2 * ig_s_190[k];
    }

#pragma omp simd aligned(t_191, t_192, t_193, pa_x, pb_y, pb_z, gg_s_192, gg_192, hf_76, \
                         hf_89, hg_192, ig_s_191, ig_s_192, ig_s_193, if__126, \
                         if__129 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_191[k] = f_6 * hf_76[k]
                   + f_2 * ig_s_191[k]
                   + pb_z[k] * if__126[k];

        t_192[k] = -f_11 * gg_s_192[k]
                   + f_5 * gg_192[k]
                   + pa_x[k] * hg_192[k]
                   + f_2 * ig_s_192[k];

        t_193[k] = f_6 * hf_89[k]
                   + f_2 * ig_s_193[k]
                   + pb_y[k] * if__129[k];
    }

#pragma omp simd aligned(t_194, t_195, t_196, pa_x, pa_y, pb_y, gg_s_194, gg_194, hf_90, \
                         hg_135, hg_194, ig_s_194, ig_s_195, ig_s_196, \
                         if__130 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_194[k] = -f_11 * gg_s_194[k]
                   + f_5 * gg_194[k]
                   + pa_x[k] * hg_194[k]
                   + f_2 * ig_s_194[k];

        t_195[k] = pa_y[k] * hg_135[k]
                   + f_2 * ig_s_195[k];

        t_196[k] = f_5 * hf_90[k]
                   + f_2 * ig_s_196[k]
                   + pb_y[k] * if__130[k];
    }

#pragma omp simd aligned(t_197, t_198, t_199, t_200, pa_y, pb_y, hf_91, hf_92, hg_137, hg_138, \
                         hg_140, ig_s_197, ig_s_198, ig_s_199, ig_s_200, \
                         if__132 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_197[k] = pa_y[k] * hg_137[k]
                   + f_2 * ig_s_197[k];

        t_198[k] = f_6 * hf_91[k]
                   + pa_y[k] * hg_138[k]
                   + f_2 * ig_s_198[k];

        t_199[k] = f_5 * hf_92[k]
                   + f_2 * ig_s_199[k]
                   + pb_y[k] * if__132[k];

        t_200[k] = pa_y[k] * hg_140[k]
                   + f_2 * ig_s_200[k];
    }

#pragma omp simd aligned(t_201, t_202, t_203, pb_x, hf_136, hf_137, hf_138, ig_s_201, \
                         ig_s_202, ig_s_203, if__136, if__137, \
                         if__138 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_201[k] = f_6 * hf_136[k]
                   + f_2 * ig_s_201[k]
                   + pb_x[k] * if__136[k];

        t_202[k] = f_6 * hf_137[k]
                   + f_2 * ig_s_202[k]
                   + pb_x[k] * if__137[k];

        t_203[k] = f_6 * hf_138[k]
                   + f_2 * ig_s_203[k]
                   + pb_x[k] * if__138[k];
    }
}

static auto
compute_prim_ig_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gg_s,
                                        const size_t gg, const size_t hf, const size_t hg,
                                        const size_t id_s, const size_t ig_s, const size_t id,
                                        const size_t if_, const size_t ncols, const double alpha,
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
    const auto f_7 = 2.5 / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;

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

    const auto *gg_s_75 = buffer.data(gg_s + 75);
    const auto *gg_s_205 = buffer.data(gg_s + 205);
    const auto *gg_s_207 = buffer.data(gg_s + 207);
    const auto *gg_s_224 = buffer.data(gg_s + 224);

    const auto *gg_75 = buffer.data(gg + 75);
    const auto *gg_205 = buffer.data(gg + 205);
    const auto *gg_207 = buffer.data(gg + 207);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_110 = buffer.data(hf + 110);
    const auto *hf_112 = buffer.data(hf + 112);
    const auto *hf_120 = buffer.data(hf + 120);
    const auto *hf_122 = buffer.data(hf + 122);
    const auto *hf_130 = buffer.data(hf + 130);
    const auto *hf_132 = buffer.data(hf + 132);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_146 = buffer.data(hf + 146);
    const auto *hf_147 = buffer.data(hf + 147);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_153 = buffer.data(hf + 153);
    const auto *hf_155 = buffer.data(hf + 155);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_165 = buffer.data(hf + 165);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_173 = buffer.data(hf + 173);
    const auto *hf_175 = buffer.data(hf + 175);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_183 = buffer.data(hf + 183);
    const auto *hf_185 = buffer.data(hf + 185);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_193 = buffer.data(hf + 193);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_203 = buffer.data(hf + 203);
    const auto *hf_205 = buffer.data(hf + 205);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *hg_135 = buffer.data(hg + 135);
    const auto *hg_144 = buffer.data(hg + 144);
    const auto *hg_149 = buffer.data(hg + 149);
    const auto *hg_150 = buffer.data(hg + 150);
    const auto *hg_151 = buffer.data(hg + 151);
    const auto *hg_153 = buffer.data(hg + 153);
    const auto *hg_156 = buffer.data(hg + 156);
    const auto *hg_205 = buffer.data(hg + 205);
    const auto *hg_207 = buffer.data(hg + 207);
    const auto *hg_210 = buffer.data(hg + 210);
    const auto *hg_212 = buffer.data(hg + 212);
    const auto *hg_215 = buffer.data(hg + 215);
    const auto *hg_219 = buffer.data(hg + 219);
    const auto *hg_224 = buffer.data(hg + 224);
    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_230 = buffer.data(hg + 230);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_238 = buffer.data(hg + 238);
    const auto *hg_239 = buffer.data(hg + 239);
    const auto *hg_245 = buffer.data(hg + 245);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_251 = buffer.data(hg + 251);
    const auto *hg_252 = buffer.data(hg + 252);
    const auto *hg_253 = buffer.data(hg + 253);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_255 = buffer.data(hg + 255);
    const auto *hg_258 = buffer.data(hg + 258);
    const auto *hg_260 = buffer.data(hg + 260);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_266 = buffer.data(hg + 266);
    const auto *hg_267 = buffer.data(hg + 267);
    const auto *hg_268 = buffer.data(hg + 268);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_270 = buffer.data(hg + 270);
    const auto *hg_273 = buffer.data(hg + 273);
    const auto *hg_275 = buffer.data(hg + 275);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_281 = buffer.data(hg + 281);
    const auto *hg_282 = buffer.data(hg + 282);
    const auto *hg_283 = buffer.data(hg + 283);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_288 = buffer.data(hg + 288);
    const auto *hg_295 = buffer.data(hg + 295);
    const auto *hg_296 = buffer.data(hg + 296);
    const auto *hg_297 = buffer.data(hg + 297);
    const auto *hg_298 = buffer.data(hg + 298);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_311 = buffer.data(hg + 311);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_314 = buffer.data(hg + 314);

    const auto *id_s_84 = buffer.data(id_s + 84);
    const auto *id_s_87 = buffer.data(id_s + 87);
    const auto *id_s_88 = buffer.data(id_s + 88);
    const auto *id_s_89 = buffer.data(id_s + 89);

    const auto *ig_s_204 = buffer.data(ig_s + 204);
    const auto *ig_s_205 = buffer.data(ig_s + 205);
    const auto *ig_s_206 = buffer.data(ig_s + 206);
    const auto *ig_s_207 = buffer.data(ig_s + 207);
    const auto *ig_s_208 = buffer.data(ig_s + 208);
    const auto *ig_s_209 = buffer.data(ig_s + 209);
    const auto *ig_s_210 = buffer.data(ig_s + 210);
    const auto *ig_s_211 = buffer.data(ig_s + 211);
    const auto *ig_s_212 = buffer.data(ig_s + 212);
    const auto *ig_s_213 = buffer.data(ig_s + 213);
    const auto *ig_s_214 = buffer.data(ig_s + 214);
    const auto *ig_s_215 = buffer.data(ig_s + 215);
    const auto *ig_s_216 = buffer.data(ig_s + 216);
    const auto *ig_s_217 = buffer.data(ig_s + 217);
    const auto *ig_s_218 = buffer.data(ig_s + 218);
    const auto *ig_s_219 = buffer.data(ig_s + 219);
    const auto *ig_s_220 = buffer.data(ig_s + 220);
    const auto *ig_s_221 = buffer.data(ig_s + 221);
    const auto *ig_s_222 = buffer.data(ig_s + 222);
    const auto *ig_s_223 = buffer.data(ig_s + 223);
    const auto *ig_s_224 = buffer.data(ig_s + 224);
    const auto *ig_s_225 = buffer.data(ig_s + 225);
    const auto *ig_s_226 = buffer.data(ig_s + 226);
    const auto *ig_s_227 = buffer.data(ig_s + 227);
    const auto *ig_s_228 = buffer.data(ig_s + 228);
    const auto *ig_s_229 = buffer.data(ig_s + 229);
    const auto *ig_s_230 = buffer.data(ig_s + 230);
    const auto *ig_s_231 = buffer.data(ig_s + 231);
    const auto *ig_s_232 = buffer.data(ig_s + 232);
    const auto *ig_s_233 = buffer.data(ig_s + 233);
    const auto *ig_s_234 = buffer.data(ig_s + 234);
    const auto *ig_s_235 = buffer.data(ig_s + 235);
    const auto *ig_s_236 = buffer.data(ig_s + 236);
    const auto *ig_s_237 = buffer.data(ig_s + 237);
    const auto *ig_s_238 = buffer.data(ig_s + 238);
    const auto *ig_s_239 = buffer.data(ig_s + 239);
    const auto *ig_s_240 = buffer.data(ig_s + 240);
    const auto *ig_s_241 = buffer.data(ig_s + 241);
    const auto *ig_s_242 = buffer.data(ig_s + 242);
    const auto *ig_s_243 = buffer.data(ig_s + 243);
    const auto *ig_s_244 = buffer.data(ig_s + 244);
    const auto *ig_s_245 = buffer.data(ig_s + 245);
    const auto *ig_s_246 = buffer.data(ig_s + 246);
    const auto *ig_s_247 = buffer.data(ig_s + 247);
    const auto *ig_s_248 = buffer.data(ig_s + 248);
    const auto *ig_s_249 = buffer.data(ig_s + 249);
    const auto *ig_s_250 = buffer.data(ig_s + 250);
    const auto *ig_s_251 = buffer.data(ig_s + 251);
    const auto *ig_s_252 = buffer.data(ig_s + 252);
    const auto *ig_s_253 = buffer.data(ig_s + 253);
    const auto *ig_s_254 = buffer.data(ig_s + 254);
    const auto *ig_s_255 = buffer.data(ig_s + 255);
    const auto *ig_s_256 = buffer.data(ig_s + 256);
    const auto *ig_s_257 = buffer.data(ig_s + 257);
    const auto *ig_s_258 = buffer.data(ig_s + 258);
    const auto *ig_s_259 = buffer.data(ig_s + 259);
    const auto *ig_s_260 = buffer.data(ig_s + 260);
    const auto *ig_s_261 = buffer.data(ig_s + 261);
    const auto *ig_s_262 = buffer.data(ig_s + 262);
    const auto *ig_s_263 = buffer.data(ig_s + 263);
    const auto *ig_s_264 = buffer.data(ig_s + 264);
    const auto *ig_s_265 = buffer.data(ig_s + 265);
    const auto *ig_s_266 = buffer.data(ig_s + 266);
    const auto *ig_s_267 = buffer.data(ig_s + 267);
    const auto *ig_s_268 = buffer.data(ig_s + 268);
    const auto *ig_s_269 = buffer.data(ig_s + 269);
    const auto *ig_s_270 = buffer.data(ig_s + 270);
    const auto *ig_s_271 = buffer.data(ig_s + 271);
    const auto *ig_s_272 = buffer.data(ig_s + 272);
    const auto *ig_s_273 = buffer.data(ig_s + 273);
    const auto *ig_s_274 = buffer.data(ig_s + 274);
    const auto *ig_s_275 = buffer.data(ig_s + 275);
    const auto *ig_s_276 = buffer.data(ig_s + 276);
    const auto *ig_s_277 = buffer.data(ig_s + 277);
    const auto *ig_s_278 = buffer.data(ig_s + 278);
    const auto *ig_s_279 = buffer.data(ig_s + 279);
    const auto *ig_s_280 = buffer.data(ig_s + 280);
    const auto *ig_s_281 = buffer.data(ig_s + 281);
    const auto *ig_s_282 = buffer.data(ig_s + 282);
    const auto *ig_s_283 = buffer.data(ig_s + 283);
    const auto *ig_s_284 = buffer.data(ig_s + 284);
    const auto *ig_s_285 = buffer.data(ig_s + 285);
    const auto *ig_s_286 = buffer.data(ig_s + 286);
    const auto *ig_s_287 = buffer.data(ig_s + 287);
    const auto *ig_s_288 = buffer.data(ig_s + 288);
    const auto *ig_s_289 = buffer.data(ig_s + 289);
    const auto *ig_s_290 = buffer.data(ig_s + 290);
    const auto *ig_s_291 = buffer.data(ig_s + 291);
    const auto *ig_s_292 = buffer.data(ig_s + 292);
    const auto *ig_s_293 = buffer.data(ig_s + 293);
    const auto *ig_s_294 = buffer.data(ig_s + 294);
    const auto *ig_s_295 = buffer.data(ig_s + 295);
    const auto *ig_s_296 = buffer.data(ig_s + 296);
    const auto *ig_s_297 = buffer.data(ig_s + 297);
    const auto *ig_s_298 = buffer.data(ig_s + 298);
    const auto *ig_s_299 = buffer.data(ig_s + 299);
    const auto *ig_s_300 = buffer.data(ig_s + 300);
    const auto *ig_s_301 = buffer.data(ig_s + 301);
    const auto *ig_s_302 = buffer.data(ig_s + 302);
    const auto *ig_s_303 = buffer.data(ig_s + 303);
    const auto *ig_s_304 = buffer.data(ig_s + 304);
    const auto *ig_s_305 = buffer.data(ig_s + 305);
    const auto *ig_s_306 = buffer.data(ig_s + 306);
    const auto *ig_s_307 = buffer.data(ig_s + 307);
    const auto *ig_s_308 = buffer.data(ig_s + 308);
    const auto *ig_s_309 = buffer.data(ig_s + 309);
    const auto *ig_s_310 = buffer.data(ig_s + 310);
    const auto *ig_s_311 = buffer.data(ig_s + 311);
    const auto *ig_s_312 = buffer.data(ig_s + 312);
    const auto *ig_s_313 = buffer.data(ig_s + 313);
    const auto *ig_s_314 = buffer.data(ig_s + 314);

    const auto *id_84 = buffer.data(id + 84);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);

    const auto *if__136 = buffer.data(if_ + 136);
    const auto *if__139 = buffer.data(if_ + 139);
    const auto *if__140 = buffer.data(if_ + 140);
    const auto *if__141 = buffer.data(if_ + 141);
    const auto *if__142 = buffer.data(if_ + 142);
    const auto *if__145 = buffer.data(if_ + 145);
    const auto *if__146 = buffer.data(if_ + 146);
    const auto *if__147 = buffer.data(if_ + 147);
    const auto *if__148 = buffer.data(if_ + 148);
    const auto *if__149 = buffer.data(if_ + 149);
    const auto *if__150 = buffer.data(if_ + 150);
    const auto *if__151 = buffer.data(if_ + 151);
    const auto *if__153 = buffer.data(if_ + 153);
    const auto *if__156 = buffer.data(if_ + 156);
    const auto *if__158 = buffer.data(if_ + 158);
    const auto *if__159 = buffer.data(if_ + 159);
    const auto *if__160 = buffer.data(if_ + 160);
    const auto *if__162 = buffer.data(if_ + 162);
    const auto *if__167 = buffer.data(if_ + 167);
    const auto *if__168 = buffer.data(if_ + 168);
    const auto *if__169 = buffer.data(if_ + 169);
    const auto *if__170 = buffer.data(if_ + 170);
    const auto *if__172 = buffer.data(if_ + 172);
    const auto *if__176 = buffer.data(if_ + 176);
    const auto *if__177 = buffer.data(if_ + 177);
    const auto *if__178 = buffer.data(if_ + 178);
    const auto *if__179 = buffer.data(if_ + 179);
    const auto *if__180 = buffer.data(if_ + 180);
    const auto *if__182 = buffer.data(if_ + 182);
    const auto *if__186 = buffer.data(if_ + 186);
    const auto *if__187 = buffer.data(if_ + 187);
    const auto *if__188 = buffer.data(if_ + 188);
    const auto *if__189 = buffer.data(if_ + 189);
    const auto *if__190 = buffer.data(if_ + 190);
    const auto *if__192 = buffer.data(if_ + 192);
    const auto *if__196 = buffer.data(if_ + 196);
    const auto *if__197 = buffer.data(if_ + 197);
    const auto *if__198 = buffer.data(if_ + 198);
    const auto *if__200 = buffer.data(if_ + 200);
    const auto *if__202 = buffer.data(if_ + 202);
    const auto *if__205 = buffer.data(if_ + 205);
    const auto *if__206 = buffer.data(if_ + 206);
    const auto *if__207 = buffer.data(if_ + 207);
    const auto *if__209 = buffer.data(if_ + 209);

#pragma omp simd aligned(t_204, t_205, t_206, pa_x, pa_y, pb_z, gg_s_205, gg_205, hf_86, \
                         hg_144, hg_205, ig_s_204, ig_s_205, ig_s_206, \
                         if__136 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_204[k] = pa_y[k] * hg_144[k]
                   + f_2 * ig_s_204[k];

        t_205[k] = -f_11 * gg_s_205[k]
                   + f_5 * gg_205[k]
                   + pa_x[k] * hg_205[k]
                   + f_2 * ig_s_205[k];

        t_206[k] = f_3 * hf_86[k]
                   + f_2 * ig_s_206[k]
                   + pb_z[k] * if__136[k];
    }

#pragma omp simd aligned(t_207, t_208, t_209, pa_x, pa_y, pb_y, gg_s_207, gg_207, hf_99, \
                         hg_149, hg_207, ig_s_207, ig_s_208, ig_s_209, \
                         if__139 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_207[k] = -f_11 * gg_s_207[k]
                   + f_5 * gg_207[k]
                   + pa_x[k] * hg_207[k]
                   + f_2 * ig_s_207[k];

        t_208[k] = f_5 * hf_99[k]
                   + f_2 * ig_s_208[k]
                   + pb_y[k] * if__139[k];

        t_209[k] = pa_y[k] * hg_149[k]
                   + f_2 * ig_s_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pa_z, pb_y, pb_z, gg_s_75, gg_75, hf_90, hg_135, \
                         ig_s_210, ig_s_211, ig_s_212, if__140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_12 * gg_s_75[k]
                   + f_3 * gg_75[k]
                   + pa_z[k] * hg_135[k]
                   + f_2 * ig_s_210[k];

        t_211[k] = f_2 * ig_s_211[k]
                   + pb_y[k] * if__140[k];

        t_212[k] = f_9 * hf_90[k]
                   + f_2 * ig_s_212[k]
                   + pb_z[k] * if__140[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, pb_x, pb_y, hf_145, id_s_84, id_s_89, ig_s_213, \
                         ig_s_214, ig_s_215, id_84, id_89, if__141, if__142, \
                         if__145 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = -f_4 * id_s_84[k]
                   + f_2 * ig_s_213[k]
                   + f_5 * id_84[k]
                   + pb_y[k] * if__141[k];

        t_214[k] = f_2 * ig_s_214[k]
                   + pb_y[k] * if__142[k];

        t_215[k] = f_6 * hf_145[k]
                   - f_4 * id_s_89[k]
                   + f_2 * ig_s_215[k]
                   + f_5 * id_89[k]
                   + pb_x[k] * if__145[k];
    }

#pragma omp simd aligned(t_216, t_217, t_218, pb_x, pb_y, hf_146, hf_147, ig_s_216, ig_s_217, \
                         ig_s_218, if__145, if__146, if__147 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_216[k] = f_6 * hf_146[k]
                   + f_2 * ig_s_216[k]
                   + pb_x[k] * if__146[k];

        t_217[k] = f_6 * hf_147[k]
                   + f_2 * ig_s_217[k]
                   + pb_x[k] * if__147[k];

        t_218[k] = f_2 * ig_s_218[k]
                   + pb_y[k] * if__145[k];
    }

#pragma omp simd aligned(t_219, t_220, t_221, pb_x, pb_y, hf_149, id_s_87, id_s_88, ig_s_219, \
                         ig_s_220, ig_s_221, id_87, id_88, if__146, if__147, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_219[k] = f_6 * hf_149[k]
                   + f_2 * ig_s_219[k]
                   + pb_x[k] * if__149[k];

        t_220[k] = -f_1 * id_s_87[k]
                   + f_2 * ig_s_220[k]
                   + f_3 * id_87[k]
                   + pb_y[k] * if__146[k];

        t_221[k] = -f_10 * id_s_88[k]
                   + f_2 * ig_s_221[k]
                   + f_6 * id_88[k]
                   + pb_y[k] * if__147[k];
    }

#pragma omp simd aligned(t_222, t_223, t_224, pa_x, pb_y, gg_s_224, gg_224, hg_224, id_s_89, \
                         ig_s_222, ig_s_223, ig_s_224, id_89, if__148, \
                         if__149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_222[k] = -f_4 * id_s_89[k]
                   + f_2 * ig_s_222[k]
                   + f_5 * id_89[k]
                   + pb_y[k] * if__148[k];

        t_223[k] = f_2 * ig_s_223[k]
                   + pb_y[k] * if__149[k];

        t_224[k] = -f_11 * gg_s_224[k]
                   + f_5 * gg_224[k]
                   + pa_x[k] * hg_224[k]
                   + f_2 * ig_s_224[k];
    }

#pragma omp simd aligned(t_225, t_226, t_227, pa_x, pb_y, pb_z, hf_100, hf_150, hg_225, \
                         ig_s_225, ig_s_226, ig_s_227, if__150 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_225[k] = f_9 * hf_150[k]
                   + pa_x[k] * hg_225[k]
                   + f_2 * ig_s_225[k];

        t_226[k] = f_7 * hf_100[k]
                   + f_2 * ig_s_226[k]
                   + pb_y[k] * if__150[k];

        t_227[k] = f_2 * ig_s_227[k]
                   + pb_z[k] * if__150[k];
    }

#pragma omp simd aligned(t_228, t_229, t_230, pa_x, pb_z, hf_153, hf_155, hg_228, hg_230, \
                         ig_s_228, ig_s_229, ig_s_230, if__151 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_6 * hf_153[k]
                   + pa_x[k] * hg_228[k]
                   + f_2 * ig_s_228[k];

        t_229[k] = f_2 * ig_s_229[k]
                   + pb_z[k] * if__151[k];

        t_230[k] = f_6 * hf_155[k]
                   + pa_x[k] * hg_230[k]
                   + f_2 * ig_s_230[k];
    }

#pragma omp simd aligned(t_231, t_232, t_233, pb_x, pb_z, hf_156, hf_158, ig_s_231, ig_s_232, \
                         ig_s_233, if__153, if__156, if__158 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_231[k] = f_5 * hf_156[k]
                   + f_2 * ig_s_231[k]
                   + pb_x[k] * if__156[k];

        t_232[k] = f_2 * ig_s_232[k]
                   + pb_z[k] * if__153[k];

        t_233[k] = f_5 * hf_158[k]
                   + f_2 * ig_s_233[k]
                   + pb_x[k] * if__158[k];
    }

#pragma omp simd aligned(t_234, t_235, t_236, t_237, pa_x, pb_x, pb_z, hf_159, hg_235, hg_237, \
                         ig_s_234, ig_s_235, ig_s_236, ig_s_237, if__156, \
                         if__159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_234[k] = f_5 * hf_159[k]
                   + f_2 * ig_s_234[k]
                   + pb_x[k] * if__159[k];

        t_235[k] = pa_x[k] * hg_235[k]
                   + f_2 * ig_s_235[k];

        t_236[k] = f_2 * ig_s_236[k]
                   + pb_z[k] * if__156[k];

        t_237[k] = pa_x[k] * hg_237[k]
                   + f_2 * ig_s_237[k];
    }

#pragma omp simd aligned(t_238, t_239, t_240, t_241, pa_x, pa_z, hg_150, hg_151, hg_238, \
                         hg_239, ig_s_238, ig_s_239, ig_s_240, \
                         ig_s_241 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_238[k] = pa_x[k] * hg_238[k]
                   + f_2 * ig_s_238[k];

        t_239[k] = pa_x[k] * hg_239[k]
                   + f_2 * ig_s_239[k];

        t_240[k] = pa_z[k] * hg_150[k]
                   + f_2 * ig_s_240[k];

        t_241[k] = pa_z[k] * hg_151[k]
                   + f_2 * ig_s_241[k];
    }

#pragma omp simd aligned(t_242, t_243, t_244, pa_z, pb_y, pb_z, hf_100, hf_112, hg_153, \
                         ig_s_242, ig_s_243, ig_s_244, if__160, \
                         if__162 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_242[k] = f_5 * hf_100[k]
                   + f_2 * ig_s_242[k]
                   + pb_z[k] * if__160[k];

        t_243[k] = pa_z[k] * hg_153[k]
                   + f_2 * ig_s_243[k];

        t_244[k] = f_9 * hf_112[k]
                   + f_2 * ig_s_244[k]
                   + pb_y[k] * if__162[k];
    }

#pragma omp simd aligned(t_245, t_246, t_247, pa_x, pa_z, pb_x, hf_165, hf_167, hg_156, \
                         hg_245, ig_s_245, ig_s_246, ig_s_247, \
                         if__167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_245[k] = f_6 * hf_165[k]
                   + pa_x[k] * hg_245[k]
                   + f_2 * ig_s_245[k];

        t_246[k] = pa_z[k] * hg_156[k]
                   + f_2 * ig_s_246[k];

        t_247[k] = f_5 * hf_167[k]
                   + f_2 * ig_s_247[k]
                   + pb_x[k] * if__167[k];
    }

#pragma omp simd aligned(t_248, t_249, t_250, t_251, pa_x, pb_x, hf_168, hf_169, hg_250, \
                         hg_251, ig_s_248, ig_s_249, ig_s_250, ig_s_251, if__168, \
                         if__169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_248[k] = f_5 * hf_168[k]
                   + f_2 * ig_s_248[k]
                   + pb_x[k] * if__168[k];

        t_249[k] = f_5 * hf_169[k]
                   + f_2 * ig_s_249[k]
                   + pb_x[k] * if__169[k];

        t_250[k] = pa_x[k] * hg_250[k]
                   + f_2 * ig_s_250[k];

        t_251[k] = pa_x[k] * hg_251[k]
                   + f_2 * ig_s_251[k];
    }

#pragma omp simd aligned(t_252, t_253, t_254, t_255, pa_x, hf_170, hg_252, hg_253, hg_254, \
                         hg_255, ig_s_252, ig_s_253, ig_s_254, \
                         ig_s_255 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_252[k] = pa_x[k] * hg_252[k]
                   + f_2 * ig_s_252[k];

        t_253[k] = pa_x[k] * hg_253[k]
                   + f_2 * ig_s_253[k];

        t_254[k] = pa_x[k] * hg_254[k]
                   + f_2 * ig_s_254[k];

        t_255[k] = f_9 * hf_170[k]
                   + pa_x[k] * hg_255[k]
                   + f_2 * ig_s_255[k];
    }

#pragma omp simd aligned(t_256, t_257, t_258, pa_x, pb_y, pb_z, hf_110, hf_120, hf_173, \
                         hg_258, ig_s_256, ig_s_257, ig_s_258, \
                         if__170 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_256[k] = f_3 * hf_120[k]
                   + f_2 * ig_s_256[k]
                   + pb_y[k] * if__170[k];

        t_257[k] = f_6 * hf_110[k]
                   + f_2 * ig_s_257[k]
                   + pb_z[k] * if__170[k];

        t_258[k] = f_6 * hf_173[k]
                   + pa_x[k] * hg_258[k]
                   + f_2 * ig_s_258[k];
    }

#pragma omp simd aligned(t_259, t_260, t_261, pa_x, pb_x, pb_y, hf_122, hf_175, hf_176, \
                         hg_260, ig_s_259, ig_s_260, ig_s_261, if__172, \
                         if__176 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_259[k] = f_3 * hf_122[k]
                   + f_2 * ig_s_259[k]
                   + pb_y[k] * if__172[k];

        t_260[k] = f_6 * hf_175[k]
                   + pa_x[k] * hg_260[k]
                   + f_2 * ig_s_260[k];

        t_261[k] = f_5 * hf_176[k]
                   + f_2 * ig_s_261[k]
                   + pb_x[k] * if__176[k];
    }

#pragma omp simd aligned(t_262, t_263, t_264, pb_x, hf_177, hf_178, hf_179, ig_s_262, \
                         ig_s_263, ig_s_264, if__177, if__178, \
                         if__179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_262[k] = f_5 * hf_177[k]
                   + f_2 * ig_s_262[k]
                   + pb_x[k] * if__177[k];

        t_263[k] = f_5 * hf_178[k]
                   + f_2 * ig_s_263[k]
                   + pb_x[k] * if__178[k];

        t_264[k] = f_5 * hf_179[k]
                   + f_2 * ig_s_264[k]
                   + pb_x[k] * if__179[k];
    }

#pragma omp simd aligned(t_265, t_266, t_267, t_268, t_269, pa_x, hg_265, hg_266, hg_267, \
                         hg_268, hg_269, ig_s_265, ig_s_266, ig_s_267, ig_s_268, \
                         ig_s_269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_265[k] = pa_x[k] * hg_265[k]
                   + f_2 * ig_s_265[k];

        t_266[k] = pa_x[k] * hg_266[k]
                   + f_2 * ig_s_266[k];

        t_267[k] = pa_x[k] * hg_267[k]
                   + f_2 * ig_s_267[k];

        t_268[k] = pa_x[k] * hg_268[k]
                   + f_2 * ig_s_268[k];

        t_269[k] = pa_x[k] * hg_269[k]
                   + f_2 * ig_s_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pa_x, pb_y, pb_z, hf_120, hf_130, hf_180, \
                         hg_270, ig_s_270, ig_s_271, ig_s_272, \
                         if__180 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = f_9 * hf_180[k]
                   + pa_x[k] * hg_270[k]
                   + f_2 * ig_s_270[k];

        t_271[k] = f_6 * hf_130[k]
                   + f_2 * ig_s_271[k]
                   + pb_y[k] * if__180[k];

        t_272[k] = f_3 * hf_120[k]
                   + f_2 * ig_s_272[k]
                   + pb_z[k] * if__180[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, pa_x, pb_y, hf_132, hf_183, hf_185, hg_273, \
                         hg_275, ig_s_273, ig_s_274, ig_s_275, \
                         if__182 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_6 * hf_183[k]
                   + pa_x[k] * hg_273[k]
                   + f_2 * ig_s_273[k];

        t_274[k] = f_6 * hf_132[k]
                   + f_2 * ig_s_274[k]
                   + pb_y[k] * if__182[k];

        t_275[k] = f_6 * hf_185[k]
                   + pa_x[k] * hg_275[k]
                   + f_2 * ig_s_275[k];
    }

#pragma omp simd aligned(t_276, t_277, t_278, pb_x, hf_186, hf_187, hf_188, ig_s_276, \
                         ig_s_277, ig_s_278, if__186, if__187, \
                         if__188 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_276[k] = f_5 * hf_186[k]
                   + f_2 * ig_s_276[k]
                   + pb_x[k] * if__186[k];

        t_277[k] = f_5 * hf_187[k]
                   + f_2 * ig_s_277[k]
                   + pb_x[k] * if__187[k];

        t_278[k] = f_5 * hf_188[k]
                   + f_2 * ig_s_278[k]
                   + pb_x[k] * if__188[k];
    }

#pragma omp simd aligned(t_279, t_280, t_281, t_282, pa_x, pb_x, hf_189, hg_280, hg_281, \
                         hg_282, ig_s_279, ig_s_280, ig_s_281, ig_s_282, \
                         if__189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_279[k] = f_5 * hf_189[k]
                   + f_2 * ig_s_279[k]
                   + pb_x[k] * if__189[k];

        t_280[k] = pa_x[k] * hg_280[k]
                   + f_2 * ig_s_280[k];

        t_281[k] = pa_x[k] * hg_281[k]
                   + f_2 * ig_s_281[k];

        t_282[k] = pa_x[k] * hg_282[k]
                   + f_2 * ig_s_282[k];
    }

#pragma omp simd aligned(t_283, t_284, t_285, t_286, pa_x, pa_y, pb_y, hf_140, hg_210, hg_283, \
                         hg_284, ig_s_283, ig_s_284, ig_s_285, ig_s_286, \
                         if__190 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_283[k] = pa_x[k] * hg_283[k]
                   + f_2 * ig_s_283[k];

        t_284[k] = pa_x[k] * hg_284[k]
                   + f_2 * ig_s_284[k];

        t_285[k] = pa_y[k] * hg_210[k]
                   + f_2 * ig_s_285[k];

        t_286[k] = f_5 * hf_140[k]
                   + f_2 * ig_s_286[k]
                   + pb_y[k] * if__190[k];
    }

#pragma omp simd aligned(t_287, t_288, t_289, pa_x, pa_y, pb_y, hf_142, hf_193, hg_212, \
                         hg_288, ig_s_287, ig_s_288, ig_s_289, \
                         if__192 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_287[k] = pa_y[k] * hg_212[k]
                   + f_2 * ig_s_287[k];

        t_288[k] = f_6 * hf_193[k]
                   + pa_x[k] * hg_288[k]
                   + f_2 * ig_s_288[k];

        t_289[k] = f_5 * hf_142[k]
                   + f_2 * ig_s_289[k]
                   + pb_y[k] * if__192[k];
    }

#pragma omp simd aligned(t_290, t_291, t_292, pa_y, pb_x, hf_196, hf_197, hg_215, ig_s_290, \
                         ig_s_291, ig_s_292, if__196, if__197 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_290[k] = pa_y[k] * hg_215[k]
                   + f_2 * ig_s_290[k];

        t_291[k] = f_5 * hf_196[k]
                   + f_2 * ig_s_291[k]
                   + pb_x[k] * if__196[k];

        t_292[k] = f_5 * hf_197[k]
                   + f_2 * ig_s_292[k]
                   + pb_x[k] * if__197[k];
    }

#pragma omp simd aligned(t_293, t_294, t_295, t_296, pa_x, pa_y, pb_x, hf_198, hg_219, hg_295, \
                         hg_296, ig_s_293, ig_s_294, ig_s_295, ig_s_296, \
                         if__198 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_293[k] = f_5 * hf_198[k]
                   + f_2 * ig_s_293[k]
                   + pb_x[k] * if__198[k];

        t_294[k] = pa_y[k] * hg_219[k]
                   + f_2 * ig_s_294[k];

        t_295[k] = pa_x[k] * hg_295[k]
                   + f_2 * ig_s_295[k];

        t_296[k] = pa_x[k] * hg_296[k]
                   + f_2 * ig_s_296[k];
    }

#pragma omp simd aligned(t_297, t_298, t_299, t_300, pa_x, hf_200, hg_297, hg_298, hg_299, \
                         hg_300, ig_s_297, ig_s_298, ig_s_299, \
                         ig_s_300 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_297[k] = pa_x[k] * hg_297[k]
                   + f_2 * ig_s_297[k];

        t_298[k] = pa_x[k] * hg_298[k]
                   + f_2 * ig_s_298[k];

        t_299[k] = pa_x[k] * hg_299[k]
                   + f_2 * ig_s_299[k];

        t_300[k] = f_9 * hf_200[k]
                   + pa_x[k] * hg_300[k]
                   + f_2 * ig_s_300[k];
    }

#pragma omp simd aligned(t_301, t_302, t_303, t_304, pa_x, pb_y, pb_z, hf_140, hf_203, hg_303, \
                         ig_s_301, ig_s_302, ig_s_303, ig_s_304, if__200, \
                         if__202 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_301[k] = f_2 * ig_s_301[k]
                   + pb_y[k] * if__200[k];

        t_302[k] = f_7 * hf_140[k]
                   + f_2 * ig_s_302[k]
                   + pb_z[k] * if__200[k];

        t_303[k] = f_6 * hf_203[k]
                   + pa_x[k] * hg_303[k]
                   + f_2 * ig_s_303[k];

        t_304[k] = f_2 * ig_s_304[k]
                   + pb_y[k] * if__202[k];
    }

#pragma omp simd aligned(t_305, t_306, t_307, pa_x, pb_x, hf_205, hf_206, hf_207, hg_305, \
                         ig_s_305, ig_s_306, ig_s_307, if__206, \
                         if__207 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_305[k] = f_6 * hf_205[k]
                   + pa_x[k] * hg_305[k]
                   + f_2 * ig_s_305[k];

        t_306[k] = f_5 * hf_206[k]
                   + f_2 * ig_s_306[k]
                   + pb_x[k] * if__206[k];

        t_307[k] = f_5 * hf_207[k]
                   + f_2 * ig_s_307[k]
                   + pb_x[k] * if__207[k];
    }

#pragma omp simd aligned(t_308, t_309, t_310, t_311, pa_x, pb_x, pb_y, hf_209, hg_310, hg_311, \
                         ig_s_308, ig_s_309, ig_s_310, ig_s_311, if__205, \
                         if__209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_308[k] = f_2 * ig_s_308[k]
                   + pb_y[k] * if__205[k];

        t_309[k] = f_5 * hf_209[k]
                   + f_2 * ig_s_309[k]
                   + pb_x[k] * if__209[k];

        t_310[k] = pa_x[k] * hg_310[k]
                   + f_2 * ig_s_310[k];

        t_311[k] = pa_x[k] * hg_311[k]
                   + f_2 * ig_s_311[k];
    }

#pragma omp simd aligned(t_312, t_313, t_314, pa_x, pb_y, hg_312, hg_314, ig_s_312, ig_s_313, \
                         ig_s_314, if__209 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_312[k] = pa_x[k] * hg_312[k]
                   + f_2 * ig_s_312[k];

        t_313[k] = f_2 * ig_s_313[k]
                   + pb_y[k] * if__209[k];

        t_314[k] = pa_x[k] * hg_314[k]
                   + f_2 * ig_s_314[k];
    }
}

static auto
compute_prim_ig_kinetic_energy_0_piece3(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gg_s,
                                        const size_t gg, const size_t hf, const size_t hg,
                                        const size_t id_s, const size_t ig_s, const size_t id,
                                        const size_t if_, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;
    const auto f_6 = 1.0 / p;
    const auto f_7 = 2.5 / p;
    const auto f_8 = 4.0 * beta / p;
    const auto f_9 = 2.0 / p;
    const auto f_10 = 2.0 * alpha / p;
    const auto f_11 = beta / p;
    const auto f_12 = 3.0 * beta / p;
    const auto f_13 = 2.0 * beta / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gg_s_160 = buffer.data(gg_s + 160);
    const auto *gg_s_175 = buffer.data(gg_s + 175);
    const auto *gg_s_179 = buffer.data(gg_s + 179);
    const auto *gg_s_190 = buffer.data(gg_s + 190);
    const auto *gg_s_194 = buffer.data(gg_s + 194);
    const auto *gg_s_209 = buffer.data(gg_s + 209);
    const auto *gg_s_224 = buffer.data(gg_s + 224);

    const auto *gg_160 = buffer.data(gg + 160);
    const auto *gg_175 = buffer.data(gg + 175);
    const auto *gg_179 = buffer.data(gg + 179);
    const auto *gg_190 = buffer.data(gg + 190);
    const auto *gg_194 = buffer.data(gg + 194);
    const auto *gg_209 = buffer.data(gg + 209);
    const auto *gg_224 = buffer.data(gg + 224);

    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_157 = buffer.data(hf + 157);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_208 = buffer.data(hf + 208);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *hg_225 = buffer.data(hg + 225);
    const auto *hg_226 = buffer.data(hg + 226);
    const auto *hg_228 = buffer.data(hg + 228);
    const auto *hg_229 = buffer.data(hg + 229);
    const auto *hg_235 = buffer.data(hg + 235);
    const auto *hg_237 = buffer.data(hg + 237);
    const auto *hg_250 = buffer.data(hg + 250);
    const auto *hg_254 = buffer.data(hg + 254);
    const auto *hg_265 = buffer.data(hg + 265);
    const auto *hg_269 = buffer.data(hg + 269);
    const auto *hg_280 = buffer.data(hg + 280);
    const auto *hg_284 = buffer.data(hg + 284);
    const auto *hg_299 = buffer.data(hg + 299);
    const auto *hg_300 = buffer.data(hg + 300);
    const auto *hg_301 = buffer.data(hg + 301);
    const auto *hg_302 = buffer.data(hg + 302);
    const auto *hg_303 = buffer.data(hg + 303);
    const auto *hg_304 = buffer.data(hg + 304);
    const auto *hg_305 = buffer.data(hg + 305);
    const auto *hg_310 = buffer.data(hg + 310);
    const auto *hg_312 = buffer.data(hg + 312);
    const auto *hg_314 = buffer.data(hg + 314);

    const auto *id_s_126 = buffer.data(id_s + 126);
    const auto *id_s_127 = buffer.data(id_s + 127);
    const auto *id_s_129 = buffer.data(id_s + 129);
    const auto *id_s_131 = buffer.data(id_s + 131);
    const auto *id_s_134 = buffer.data(id_s + 134);
    const auto *id_s_137 = buffer.data(id_s + 137);
    const auto *id_s_138 = buffer.data(id_s + 138);
    const auto *id_s_139 = buffer.data(id_s + 139);
    const auto *id_s_140 = buffer.data(id_s + 140);
    const auto *id_s_141 = buffer.data(id_s + 141);
    const auto *id_s_142 = buffer.data(id_s + 142);
    const auto *id_s_143 = buffer.data(id_s + 143);
    const auto *id_s_144 = buffer.data(id_s + 144);
    const auto *id_s_145 = buffer.data(id_s + 145);
    const auto *id_s_146 = buffer.data(id_s + 146);
    const auto *id_s_147 = buffer.data(id_s + 147);
    const auto *id_s_148 = buffer.data(id_s + 148);
    const auto *id_s_149 = buffer.data(id_s + 149);
    const auto *id_s_150 = buffer.data(id_s + 150);
    const auto *id_s_151 = buffer.data(id_s + 151);
    const auto *id_s_152 = buffer.data(id_s + 152);
    const auto *id_s_153 = buffer.data(id_s + 153);
    const auto *id_s_154 = buffer.data(id_s + 154);
    const auto *id_s_155 = buffer.data(id_s + 155);
    const auto *id_s_162 = buffer.data(id_s + 162);
    const auto *id_s_164 = buffer.data(id_s + 164);
    const auto *id_s_165 = buffer.data(id_s + 165);
    const auto *id_s_166 = buffer.data(id_s + 166);
    const auto *id_s_167 = buffer.data(id_s + 167);

    const auto *ig_s_315 = buffer.data(ig_s + 315);
    const auto *ig_s_316 = buffer.data(ig_s + 316);
    const auto *ig_s_317 = buffer.data(ig_s + 317);
    const auto *ig_s_318 = buffer.data(ig_s + 318);
    const auto *ig_s_319 = buffer.data(ig_s + 319);
    const auto *ig_s_320 = buffer.data(ig_s + 320);
    const auto *ig_s_321 = buffer.data(ig_s + 321);
    const auto *ig_s_322 = buffer.data(ig_s + 322);
    const auto *ig_s_323 = buffer.data(ig_s + 323);
    const auto *ig_s_324 = buffer.data(ig_s + 324);
    const auto *ig_s_325 = buffer.data(ig_s + 325);
    const auto *ig_s_326 = buffer.data(ig_s + 326);
    const auto *ig_s_327 = buffer.data(ig_s + 327);
    const auto *ig_s_328 = buffer.data(ig_s + 328);
    const auto *ig_s_329 = buffer.data(ig_s + 329);
    const auto *ig_s_330 = buffer.data(ig_s + 330);
    const auto *ig_s_331 = buffer.data(ig_s + 331);
    const auto *ig_s_332 = buffer.data(ig_s + 332);
    const auto *ig_s_333 = buffer.data(ig_s + 333);
    const auto *ig_s_334 = buffer.data(ig_s + 334);
    const auto *ig_s_335 = buffer.data(ig_s + 335);
    const auto *ig_s_336 = buffer.data(ig_s + 336);
    const auto *ig_s_337 = buffer.data(ig_s + 337);
    const auto *ig_s_338 = buffer.data(ig_s + 338);
    const auto *ig_s_339 = buffer.data(ig_s + 339);
    const auto *ig_s_340 = buffer.data(ig_s + 340);
    const auto *ig_s_341 = buffer.data(ig_s + 341);
    const auto *ig_s_342 = buffer.data(ig_s + 342);
    const auto *ig_s_343 = buffer.data(ig_s + 343);
    const auto *ig_s_344 = buffer.data(ig_s + 344);
    const auto *ig_s_345 = buffer.data(ig_s + 345);
    const auto *ig_s_346 = buffer.data(ig_s + 346);
    const auto *ig_s_347 = buffer.data(ig_s + 347);
    const auto *ig_s_348 = buffer.data(ig_s + 348);
    const auto *ig_s_349 = buffer.data(ig_s + 349);
    const auto *ig_s_350 = buffer.data(ig_s + 350);
    const auto *ig_s_351 = buffer.data(ig_s + 351);
    const auto *ig_s_352 = buffer.data(ig_s + 352);
    const auto *ig_s_353 = buffer.data(ig_s + 353);
    const auto *ig_s_354 = buffer.data(ig_s + 354);
    const auto *ig_s_355 = buffer.data(ig_s + 355);
    const auto *ig_s_356 = buffer.data(ig_s + 356);
    const auto *ig_s_357 = buffer.data(ig_s + 357);
    const auto *ig_s_358 = buffer.data(ig_s + 358);
    const auto *ig_s_359 = buffer.data(ig_s + 359);
    const auto *ig_s_360 = buffer.data(ig_s + 360);
    const auto *ig_s_361 = buffer.data(ig_s + 361);
    const auto *ig_s_362 = buffer.data(ig_s + 362);
    const auto *ig_s_363 = buffer.data(ig_s + 363);
    const auto *ig_s_364 = buffer.data(ig_s + 364);
    const auto *ig_s_365 = buffer.data(ig_s + 365);
    const auto *ig_s_366 = buffer.data(ig_s + 366);
    const auto *ig_s_367 = buffer.data(ig_s + 367);
    const auto *ig_s_368 = buffer.data(ig_s + 368);
    const auto *ig_s_369 = buffer.data(ig_s + 369);
    const auto *ig_s_370 = buffer.data(ig_s + 370);
    const auto *ig_s_371 = buffer.data(ig_s + 371);
    const auto *ig_s_372 = buffer.data(ig_s + 372);
    const auto *ig_s_373 = buffer.data(ig_s + 373);
    const auto *ig_s_374 = buffer.data(ig_s + 374);
    const auto *ig_s_375 = buffer.data(ig_s + 375);
    const auto *ig_s_376 = buffer.data(ig_s + 376);
    const auto *ig_s_377 = buffer.data(ig_s + 377);
    const auto *ig_s_378 = buffer.data(ig_s + 378);
    const auto *ig_s_379 = buffer.data(ig_s + 379);
    const auto *ig_s_380 = buffer.data(ig_s + 380);
    const auto *ig_s_381 = buffer.data(ig_s + 381);
    const auto *ig_s_382 = buffer.data(ig_s + 382);
    const auto *ig_s_383 = buffer.data(ig_s + 383);
    const auto *ig_s_384 = buffer.data(ig_s + 384);
    const auto *ig_s_385 = buffer.data(ig_s + 385);
    const auto *ig_s_386 = buffer.data(ig_s + 386);
    const auto *ig_s_387 = buffer.data(ig_s + 387);
    const auto *ig_s_388 = buffer.data(ig_s + 388);
    const auto *ig_s_389 = buffer.data(ig_s + 389);
    const auto *ig_s_390 = buffer.data(ig_s + 390);
    const auto *ig_s_391 = buffer.data(ig_s + 391);
    const auto *ig_s_392 = buffer.data(ig_s + 392);
    const auto *ig_s_393 = buffer.data(ig_s + 393);
    const auto *ig_s_394 = buffer.data(ig_s + 394);
    const auto *ig_s_395 = buffer.data(ig_s + 395);
    const auto *ig_s_396 = buffer.data(ig_s + 396);
    const auto *ig_s_397 = buffer.data(ig_s + 397);
    const auto *ig_s_398 = buffer.data(ig_s + 398);
    const auto *ig_s_399 = buffer.data(ig_s + 399);
    const auto *ig_s_400 = buffer.data(ig_s + 400);
    const auto *ig_s_401 = buffer.data(ig_s + 401);
    const auto *ig_s_402 = buffer.data(ig_s + 402);
    const auto *ig_s_403 = buffer.data(ig_s + 403);
    const auto *ig_s_404 = buffer.data(ig_s + 404);
    const auto *ig_s_405 = buffer.data(ig_s + 405);
    const auto *ig_s_406 = buffer.data(ig_s + 406);
    const auto *ig_s_407 = buffer.data(ig_s + 407);
    const auto *ig_s_408 = buffer.data(ig_s + 408);
    const auto *ig_s_409 = buffer.data(ig_s + 409);
    const auto *ig_s_410 = buffer.data(ig_s + 410);
    const auto *ig_s_411 = buffer.data(ig_s + 411);
    const auto *ig_s_412 = buffer.data(ig_s + 412);
    const auto *ig_s_413 = buffer.data(ig_s + 413);
    const auto *ig_s_414 = buffer.data(ig_s + 414);
    const auto *ig_s_415 = buffer.data(ig_s + 415);
    const auto *ig_s_416 = buffer.data(ig_s + 416);

    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_137 = buffer.data(id + 137);
    const auto *id_138 = buffer.data(id + 138);
    const auto *id_139 = buffer.data(id + 139);
    const auto *id_140 = buffer.data(id + 140);
    const auto *id_141 = buffer.data(id + 141);
    const auto *id_142 = buffer.data(id + 142);
    const auto *id_143 = buffer.data(id + 143);
    const auto *id_144 = buffer.data(id + 144);
    const auto *id_145 = buffer.data(id + 145);
    const auto *id_146 = buffer.data(id + 146);
    const auto *id_147 = buffer.data(id + 147);
    const auto *id_148 = buffer.data(id + 148);
    const auto *id_149 = buffer.data(id + 149);
    const auto *id_150 = buffer.data(id + 150);
    const auto *id_151 = buffer.data(id + 151);
    const auto *id_152 = buffer.data(id + 152);
    const auto *id_153 = buffer.data(id + 153);
    const auto *id_154 = buffer.data(id + 154);
    const auto *id_155 = buffer.data(id + 155);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

    const auto *if__210 = buffer.data(if_ + 210);
    const auto *if__211 = buffer.data(if_ + 211);
    const auto *if__213 = buffer.data(if_ + 213);
    const auto *if__215 = buffer.data(if_ + 215);
    const auto *if__216 = buffer.data(if_ + 216);
    const auto *if__217 = buffer.data(if_ + 217);
    const auto *if__218 = buffer.data(if_ + 218);
    const auto *if__219 = buffer.data(if_ + 219);
    const auto *if__222 = buffer.data(if_ + 222);
    const auto *if__225 = buffer.data(if_ + 225);
    const auto *if__226 = buffer.data(if_ + 226);
    const auto *if__227 = buffer.data(if_ + 227);
    const auto *if__228 = buffer.data(if_ + 228);
    const auto *if__229 = buffer.data(if_ + 229);
    const auto *if__230 = buffer.data(if_ + 230);
    const auto *if__231 = buffer.data(if_ + 231);
    const auto *if__232 = buffer.data(if_ + 232);
    const auto *if__233 = buffer.data(if_ + 233);
    const auto *if__234 = buffer.data(if_ + 234);
    const auto *if__235 = buffer.data(if_ + 235);
    const auto *if__236 = buffer.data(if_ + 236);
    const auto *if__237 = buffer.data(if_ + 237);
    const auto *if__238 = buffer.data(if_ + 238);
    const auto *if__239 = buffer.data(if_ + 239);
    const auto *if__240 = buffer.data(if_ + 240);
    const auto *if__241 = buffer.data(if_ + 241);
    const auto *if__242 = buffer.data(if_ + 242);
    const auto *if__243 = buffer.data(if_ + 243);
    const auto *if__244 = buffer.data(if_ + 244);
    const auto *if__245 = buffer.data(if_ + 245);
    const auto *if__246 = buffer.data(if_ + 246);
    const auto *if__247 = buffer.data(if_ + 247);
    const auto *if__248 = buffer.data(if_ + 248);
    const auto *if__249 = buffer.data(if_ + 249);
    const auto *if__250 = buffer.data(if_ + 250);
    const auto *if__251 = buffer.data(if_ + 251);
    const auto *if__252 = buffer.data(if_ + 252);
    const auto *if__253 = buffer.data(if_ + 253);
    const auto *if__254 = buffer.data(if_ + 254);
    const auto *if__255 = buffer.data(if_ + 255);
    const auto *if__256 = buffer.data(if_ + 256);
    const auto *if__257 = buffer.data(if_ + 257);
    const auto *if__258 = buffer.data(if_ + 258);
    const auto *if__259 = buffer.data(if_ + 259);
    const auto *if__266 = buffer.data(if_ + 266);
    const auto *if__267 = buffer.data(if_ + 267);
    const auto *if__268 = buffer.data(if_ + 268);
    const auto *if__269 = buffer.data(if_ + 269);
    const auto *if__270 = buffer.data(if_ + 270);
    const auto *if__272 = buffer.data(if_ + 272);
    const auto *if__273 = buffer.data(if_ + 273);
    const auto *if__275 = buffer.data(if_ + 275);
    const auto *if__276 = buffer.data(if_ + 276);
    const auto *if__277 = buffer.data(if_ + 277);
    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(t_315, t_316, t_317, pb_x, pb_z, id_s_126, id_s_127, ig_s_315, \
                         ig_s_316, ig_s_317, id_126, id_127, if__210, \
                         if__211 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_315[k] = -f_1 * id_s_126[k]
                   + f_2 * ig_s_315[k]
                   + f_3 * id_126[k]
                   + pb_x[k] * if__210[k];

        t_316[k] = -f_10 * id_s_127[k]
                   + f_2 * ig_s_316[k]
                   + f_6 * id_127[k]
                   + pb_x[k] * if__211[k];

        t_317[k] = f_2 * ig_s_317[k]
                   + pb_z[k] * if__210[k];
    }

#pragma omp simd aligned(t_318, t_319, t_320, pb_x, pb_z, id_s_129, id_s_131, ig_s_318, \
                         ig_s_319, ig_s_320, id_129, id_131, if__211, if__213, \
                         if__215 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_318[k] = -f_4 * id_s_129[k]
                   + f_2 * ig_s_318[k]
                   + f_5 * id_129[k]
                   + pb_x[k] * if__213[k];

        t_319[k] = f_2 * ig_s_319[k]
                   + pb_z[k] * if__211[k];

        t_320[k] = -f_4 * id_s_131[k]
                   + f_2 * ig_s_320[k]
                   + f_5 * id_131[k]
                   + pb_x[k] * if__215[k];
    }

#pragma omp simd aligned(t_321, t_322, t_323, t_324, pb_x, ig_s_321, ig_s_322, ig_s_323, \
                         ig_s_324, if__216, if__217, if__218, if__219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_321[k] = f_2 * ig_s_321[k]
                   + pb_x[k] * if__216[k];

        t_322[k] = f_2 * ig_s_322[k]
                   + pb_x[k] * if__217[k];

        t_323[k] = f_2 * ig_s_323[k]
                   + pb_x[k] * if__218[k];

        t_324[k] = f_2 * ig_s_324[k]
                   + pb_x[k] * if__219[k];
    }

#pragma omp simd aligned(t_325, t_326, t_327, pb_y, pb_z, hf_156, id_s_129, ig_s_325, \
                         ig_s_326, ig_s_327, id_129, if__216, if__217 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_325[k] = f_0 * hf_156[k]
                   - f_1 * id_s_129[k]
                   + f_2 * ig_s_325[k]
                   + f_3 * id_129[k]
                   + pb_y[k] * if__216[k];

        t_326[k] = f_2 * ig_s_326[k]
                   + pb_z[k] * if__216[k];

        t_327[k] = -f_4 * id_s_129[k]
                   + f_2 * ig_s_327[k]
                   + f_5 * id_129[k]
                   + pb_z[k] * if__217[k];
    }

#pragma omp simd aligned(t_328, t_329, t_330, pa_z, pb_y, pb_z, hf_159, hg_225, id_s_131, \
                         ig_s_328, ig_s_329, ig_s_330, id_131, \
                         if__219 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_328[k] = f_0 * hf_159[k]
                   + f_2 * ig_s_328[k]
                   + pb_y[k] * if__219[k];

        t_329[k] = -f_1 * id_s_131[k]
                   + f_2 * ig_s_329[k]
                   + f_3 * id_131[k]
                   + pb_z[k] * if__219[k];

        t_330[k] = pa_z[k] * hg_225[k]
                   + f_2 * ig_s_330[k];
    }

#pragma omp simd aligned(t_331, t_332, t_333, pa_z, pb_x, hg_226, hg_228, id_s_134, ig_s_331, \
                         ig_s_332, ig_s_333, id_134, if__222 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_331[k] = pa_z[k] * hg_226[k]
                   + f_2 * ig_s_331[k];

        t_332[k] = -f_10 * id_s_134[k]
                   + f_2 * ig_s_332[k]
                   + f_6 * id_134[k]
                   + pb_x[k] * if__222[k];

        t_333[k] = pa_z[k] * hg_228[k]
                   + f_2 * ig_s_333[k];
    }

#pragma omp simd aligned(t_334, t_335, t_336, pa_z, pb_x, hf_151, hg_229, id_s_137, ig_s_334, \
                         ig_s_335, ig_s_336, id_137, if__225, if__226 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_334[k] = f_5 * hf_151[k]
                   + pa_z[k] * hg_229[k]
                   + f_2 * ig_s_334[k];

        t_335[k] = -f_4 * id_s_137[k]
                   + f_2 * ig_s_335[k]
                   + f_5 * id_137[k]
                   + pb_x[k] * if__225[k];

        t_336[k] = f_2 * ig_s_336[k]
                   + pb_x[k] * if__226[k];
    }

#pragma omp simd aligned(t_337, t_338, t_339, t_340, pa_z, pb_x, hg_235, ig_s_337, ig_s_338, \
                         ig_s_339, ig_s_340, if__227, if__228, \
                         if__229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_337[k] = f_2 * ig_s_337[k]
                   + pb_x[k] * if__227[k];

        t_338[k] = f_2 * ig_s_338[k]
                   + pb_x[k] * if__228[k];

        t_339[k] = f_2 * ig_s_339[k]
                   + pb_x[k] * if__229[k];

        t_340[k] = pa_z[k] * hg_235[k]
                   + f_2 * ig_s_340[k];
    }

#pragma omp simd aligned(t_341, t_342, t_343, pa_z, pb_y, pb_z, hf_156, hf_157, hf_169, \
                         hg_237, ig_s_341, ig_s_342, ig_s_343, if__226, \
                         if__229 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_341[k] = f_5 * hf_156[k]
                   + f_2 * ig_s_341[k]
                   + pb_z[k] * if__226[k];

        t_342[k] = f_6 * hf_157[k]
                   + pa_z[k] * hg_237[k]
                   + f_2 * ig_s_342[k];

        t_343[k] = f_7 * hf_169[k]
                   + f_2 * ig_s_343[k]
                   + pb_y[k] * if__229[k];
    }

#pragma omp simd aligned(t_344, t_345, pa_y, pb_x, gg_s_179, gg_179, hg_254, id_s_138, \
                         ig_s_344, ig_s_345, id_138, if__230 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_344[k] = -f_8 * gg_s_179[k]
                   + f_9 * gg_179[k]
                   + pa_y[k] * hg_254[k]
                   + f_2 * ig_s_344[k];

        t_345[k] = -f_1 * id_s_138[k]
                   + f_2 * ig_s_345[k]
                   + f_3 * id_138[k]
                   + pb_x[k] * if__230[k];
    }

#pragma omp simd aligned(t_346, t_347, t_348, pb_x, id_s_139, id_s_140, id_s_141, ig_s_346, \
                         ig_s_347, ig_s_348, id_139, id_140, id_141, if__231, if__232, \
                         if__233 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_346[k] = -f_10 * id_s_139[k]
                   + f_2 * ig_s_346[k]
                   + f_6 * id_139[k]
                   + pb_x[k] * if__231[k];

        t_347[k] = -f_10 * id_s_140[k]
                   + f_2 * ig_s_347[k]
                   + f_6 * id_140[k]
                   + pb_x[k] * if__232[k];

        t_348[k] = -f_4 * id_s_141[k]
                   + f_2 * ig_s_348[k]
                   + f_5 * id_141[k]
                   + pb_x[k] * if__233[k];
    }

#pragma omp simd aligned(t_349, t_350, t_351, pb_x, id_s_142, id_s_143, ig_s_349, ig_s_350, \
                         ig_s_351, id_142, id_143, if__234, if__235, \
                         if__236 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_349[k] = -f_4 * id_s_142[k]
                   + f_2 * ig_s_349[k]
                   + f_5 * id_142[k]
                   + pb_x[k] * if__234[k];

        t_350[k] = -f_4 * id_s_143[k]
                   + f_2 * ig_s_350[k]
                   + f_5 * id_143[k]
                   + pb_x[k] * if__235[k];

        t_351[k] = f_2 * ig_s_351[k]
                   + pb_x[k] * if__236[k];
    }

#pragma omp simd aligned(t_352, t_353, t_354, t_355, pa_z, pb_x, gg_s_160, gg_160, hg_250, \
                         ig_s_352, ig_s_353, ig_s_354, ig_s_355, if__237, if__238, \
                         if__239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_352[k] = f_2 * ig_s_352[k]
                   + pb_x[k] * if__237[k];

        t_353[k] = f_2 * ig_s_353[k]
                   + pb_x[k] * if__238[k];

        t_354[k] = f_2 * ig_s_354[k]
                   + pb_x[k] * if__239[k];

        t_355[k] = -f_11 * gg_s_160[k]
                   + f_5 * gg_160[k]
                   + pa_z[k] * hg_250[k]
                   + f_2 * ig_s_355[k];
    }

#pragma omp simd aligned(t_356, t_357, t_358, pb_y, pb_z, hf_166, hf_178, hf_179, id_s_143, \
                         ig_s_356, ig_s_357, ig_s_358, id_143, if__236, if__238, \
                         if__239 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_356[k] = f_6 * hf_166[k]
                   + f_2 * ig_s_356[k]
                   + pb_z[k] * if__236[k];

        t_357[k] = f_9 * hf_178[k]
                   - f_4 * id_s_143[k]
                   + f_2 * ig_s_357[k]
                   + f_5 * id_143[k]
                   + pb_y[k] * if__238[k];

        t_358[k] = f_9 * hf_179[k]
                   + f_2 * ig_s_358[k]
                   + pb_y[k] * if__239[k];
    }

#pragma omp simd aligned(t_359, t_360, pa_y, pb_x, gg_s_194, gg_194, hg_269, id_s_144, \
                         ig_s_359, ig_s_360, id_144, if__240 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_359[k] = -f_12 * gg_s_194[k]
                   + f_3 * gg_194[k]
                   + pa_y[k] * hg_269[k]
                   + f_2 * ig_s_359[k];

        t_360[k] = -f_1 * id_s_144[k]
                   + f_2 * ig_s_360[k]
                   + f_3 * id_144[k]
                   + pb_x[k] * if__240[k];
    }

#pragma omp simd aligned(t_361, t_362, t_363, pb_x, id_s_145, id_s_146, id_s_147, ig_s_361, \
                         ig_s_362, ig_s_363, id_145, id_146, id_147, if__241, if__242, \
                         if__243 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_361[k] = -f_10 * id_s_145[k]
                   + f_2 * ig_s_361[k]
                   + f_6 * id_145[k]
                   + pb_x[k] * if__241[k];

        t_362[k] = -f_10 * id_s_146[k]
                   + f_2 * ig_s_362[k]
                   + f_6 * id_146[k]
                   + pb_x[k] * if__242[k];

        t_363[k] = -f_4 * id_s_147[k]
                   + f_2 * ig_s_363[k]
                   + f_5 * id_147[k]
                   + pb_x[k] * if__243[k];
    }

#pragma omp simd aligned(t_364, t_365, t_366, pb_x, id_s_148, id_s_149, ig_s_364, ig_s_365, \
                         ig_s_366, id_148, id_149, if__244, if__245, \
                         if__246 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_364[k] = -f_4 * id_s_148[k]
                   + f_2 * ig_s_364[k]
                   + f_5 * id_148[k]
                   + pb_x[k] * if__244[k];

        t_365[k] = -f_4 * id_s_149[k]
                   + f_2 * ig_s_365[k]
                   + f_5 * id_149[k]
                   + pb_x[k] * if__245[k];

        t_366[k] = f_2 * ig_s_366[k]
                   + pb_x[k] * if__246[k];
    }

#pragma omp simd aligned(t_367, t_368, t_369, t_370, pa_z, pb_x, gg_s_175, gg_175, hg_265, \
                         ig_s_367, ig_s_368, ig_s_369, ig_s_370, if__247, if__248, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_367[k] = f_2 * ig_s_367[k]
                   + pb_x[k] * if__247[k];

        t_368[k] = f_2 * ig_s_368[k]
                   + pb_x[k] * if__248[k];

        t_369[k] = f_2 * ig_s_369[k]
                   + pb_x[k] * if__249[k];

        t_370[k] = -f_13 * gg_s_175[k]
                   + f_6 * gg_175[k]
                   + pa_z[k] * hg_265[k]
                   + f_2 * ig_s_370[k];
    }

#pragma omp simd aligned(t_371, t_372, t_373, pb_y, pb_z, hf_176, hf_188, hf_189, id_s_149, \
                         ig_s_371, ig_s_372, ig_s_373, id_149, if__246, if__248, \
                         if__249 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_371[k] = f_3 * hf_176[k]
                   + f_2 * ig_s_371[k]
                   + pb_z[k] * if__246[k];

        t_372[k] = f_3 * hf_188[k]
                   - f_4 * id_s_149[k]
                   + f_2 * ig_s_372[k]
                   + f_5 * id_149[k]
                   + pb_y[k] * if__248[k];

        t_373[k] = f_3 * hf_189[k]
                   + f_2 * ig_s_373[k]
                   + pb_y[k] * if__249[k];
    }

#pragma omp simd aligned(t_374, t_375, pa_y, pb_x, gg_s_209, gg_209, hg_284, id_s_150, \
                         ig_s_374, ig_s_375, id_150, if__250 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_374[k] = -f_13 * gg_s_209[k]
                   + f_6 * gg_209[k]
                   + pa_y[k] * hg_284[k]
                   + f_2 * ig_s_374[k];

        t_375[k] = -f_1 * id_s_150[k]
                   + f_2 * ig_s_375[k]
                   + f_3 * id_150[k]
                   + pb_x[k] * if__250[k];
    }

#pragma omp simd aligned(t_376, t_377, t_378, pb_x, id_s_151, id_s_152, id_s_153, ig_s_376, \
                         ig_s_377, ig_s_378, id_151, id_152, id_153, if__251, if__252, \
                         if__253 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_376[k] = -f_10 * id_s_151[k]
                   + f_2 * ig_s_376[k]
                   + f_6 * id_151[k]
                   + pb_x[k] * if__251[k];

        t_377[k] = -f_10 * id_s_152[k]
                   + f_2 * ig_s_377[k]
                   + f_6 * id_152[k]
                   + pb_x[k] * if__252[k];

        t_378[k] = -f_4 * id_s_153[k]
                   + f_2 * ig_s_378[k]
                   + f_5 * id_153[k]
                   + pb_x[k] * if__253[k];
    }

#pragma omp simd aligned(t_379, t_380, t_381, pb_x, id_s_154, id_s_155, ig_s_379, ig_s_380, \
                         ig_s_381, id_154, id_155, if__254, if__255, \
                         if__256 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_379[k] = -f_4 * id_s_154[k]
                   + f_2 * ig_s_379[k]
                   + f_5 * id_154[k]
                   + pb_x[k] * if__254[k];

        t_380[k] = -f_4 * id_s_155[k]
                   + f_2 * ig_s_380[k]
                   + f_5 * id_155[k]
                   + pb_x[k] * if__255[k];

        t_381[k] = f_2 * ig_s_381[k]
                   + pb_x[k] * if__256[k];
    }

#pragma omp simd aligned(t_382, t_383, t_384, t_385, pa_z, pb_x, gg_s_190, gg_190, hg_280, \
                         ig_s_382, ig_s_383, ig_s_384, ig_s_385, if__257, if__258, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_382[k] = f_2 * ig_s_382[k]
                   + pb_x[k] * if__257[k];

        t_383[k] = f_2 * ig_s_383[k]
                   + pb_x[k] * if__258[k];

        t_384[k] = f_2 * ig_s_384[k]
                   + pb_x[k] * if__259[k];

        t_385[k] = -f_12 * gg_s_190[k]
                   + f_3 * gg_190[k]
                   + pa_z[k] * hg_280[k]
                   + f_2 * ig_s_385[k];
    }

#pragma omp simd aligned(t_386, t_387, t_388, pb_y, pb_z, hf_186, hf_198, hf_199, id_s_155, \
                         ig_s_386, ig_s_387, ig_s_388, id_155, if__256, if__258, \
                         if__259 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_386[k] = f_9 * hf_186[k]
                   + f_2 * ig_s_386[k]
                   + pb_z[k] * if__256[k];

        t_387[k] = f_6 * hf_198[k]
                   - f_4 * id_s_155[k]
                   + f_2 * ig_s_387[k]
                   + f_5 * id_155[k]
                   + pb_y[k] * if__258[k];

        t_388[k] = f_6 * hf_199[k]
                   + f_2 * ig_s_388[k]
                   + pb_y[k] * if__259[k];
    }

#pragma omp simd aligned(t_389, t_390, t_391, t_392, pa_y, gg_s_224, gg_224, hf_200, hg_299, \
                         hg_300, hg_301, hg_302, ig_s_389, ig_s_390, ig_s_391, \
                         ig_s_392 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_389[k] = -f_11 * gg_s_224[k]
                   + f_5 * gg_224[k]
                   + pa_y[k] * hg_299[k]
                   + f_2 * ig_s_389[k];

        t_390[k] = pa_y[k] * hg_300[k]
                   + f_2 * ig_s_390[k];

        t_391[k] = f_5 * hf_200[k]
                   + pa_y[k] * hg_301[k]
                   + f_2 * ig_s_391[k];

        t_392[k] = pa_y[k] * hg_302[k]
                   + f_2 * ig_s_392[k];
    }

#pragma omp simd aligned(t_393, t_394, t_395, t_396, pa_y, pb_x, hf_201, hf_202, hg_303, \
                         hg_304, hg_305, ig_s_393, ig_s_394, ig_s_395, ig_s_396, \
                         if__266 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_393[k] = f_6 * hf_201[k]
                   + pa_y[k] * hg_303[k]
                   + f_2 * ig_s_393[k];

        t_394[k] = f_5 * hf_202[k]
                   + pa_y[k] * hg_304[k]
                   + f_2 * ig_s_394[k];

        t_395[k] = pa_y[k] * hg_305[k]
                   + f_2 * ig_s_395[k];

        t_396[k] = f_2 * ig_s_396[k]
                   + pb_x[k] * if__266[k];
    }

#pragma omp simd aligned(t_397, t_398, t_399, t_400, pa_y, pb_x, hf_206, hg_310, ig_s_397, \
                         ig_s_398, ig_s_399, ig_s_400, if__267, if__268, \
                         if__269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_397[k] = f_2 * ig_s_397[k]
                   + pb_x[k] * if__267[k];

        t_398[k] = f_2 * ig_s_398[k]
                   + pb_x[k] * if__268[k];

        t_399[k] = f_2 * ig_s_399[k]
                   + pb_x[k] * if__269[k];

        t_400[k] = f_9 * hf_206[k]
                   + pa_y[k] * hg_310[k]
                   + f_2 * ig_s_400[k];
    }

#pragma omp simd aligned(t_401, t_402, t_403, pa_y, pb_y, pb_z, hf_196, hf_208, hf_209, \
                         hg_312, ig_s_401, ig_s_402, ig_s_403, if__266, \
                         if__269 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_401[k] = f_7 * hf_196[k]
                   + f_2 * ig_s_401[k]
                   + pb_z[k] * if__266[k];

        t_402[k] = f_6 * hf_208[k]
                   + pa_y[k] * hg_312[k]
                   + f_2 * ig_s_402[k];

        t_403[k] = f_5 * hf_209[k]
                   + f_2 * ig_s_403[k]
                   + pb_y[k] * if__269[k];
    }

#pragma omp simd aligned(t_404, t_405, t_406, pa_y, pb_x, pb_y, hg_314, id_s_162, ig_s_404, \
                         ig_s_405, ig_s_406, id_162, if__270 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_404[k] = pa_y[k] * hg_314[k]
                   + f_2 * ig_s_404[k];

        t_405[k] = -f_1 * id_s_162[k]
                   + f_2 * ig_s_405[k]
                   + f_3 * id_162[k]
                   + pb_x[k] * if__270[k];

        t_406[k] = f_2 * ig_s_406[k]
                   + pb_y[k] * if__270[k];
    }

#pragma omp simd aligned(t_407, t_408, t_409, pb_x, pb_y, id_s_164, id_s_165, ig_s_407, \
                         ig_s_408, ig_s_409, id_164, id_165, if__272, \
                         if__273 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_407[k] = -f_10 * id_s_164[k]
                   + f_2 * ig_s_407[k]
                   + f_6 * id_164[k]
                   + pb_x[k] * if__272[k];

        t_408[k] = -f_4 * id_s_165[k]
                   + f_2 * ig_s_408[k]
                   + f_5 * id_165[k]
                   + pb_x[k] * if__273[k];

        t_409[k] = f_2 * ig_s_409[k]
                   + pb_y[k] * if__272[k];
    }

#pragma omp simd aligned(t_410, t_411, t_412, t_413, pb_x, id_s_167, ig_s_410, ig_s_411, \
                         ig_s_412, ig_s_413, id_167, if__275, if__276, if__277, \
                         if__278 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_410[k] = -f_4 * id_s_167[k]
                   + f_2 * ig_s_410[k]
                   + f_5 * id_167[k]
                   + pb_x[k] * if__275[k];

        t_411[k] = f_2 * ig_s_411[k]
                   + pb_x[k] * if__276[k];

        t_412[k] = f_2 * ig_s_412[k]
                   + pb_x[k] * if__277[k];

        t_413[k] = f_2 * ig_s_413[k]
                   + pb_x[k] * if__278[k];
    }

#pragma omp simd aligned(t_414, t_415, t_416, pb_x, pb_y, id_s_165, id_s_166, ig_s_414, \
                         ig_s_415, ig_s_416, id_165, id_166, if__276, if__277, \
                         if__279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_414[k] = f_2 * ig_s_414[k]
                   + pb_x[k] * if__279[k];

        t_415[k] = -f_1 * id_s_165[k]
                   + f_2 * ig_s_415[k]
                   + f_3 * id_165[k]
                   + pb_y[k] * if__276[k];

        t_416[k] = -f_10 * id_s_166[k]
                   + f_2 * ig_s_416[k]
                   + f_6 * id_166[k]
                   + pb_y[k] * if__277[k];
    }
}

static auto
compute_prim_ig_kinetic_energy_0_piece4(CSimdMatrix &buffer, const size_t target,
                                        const size_t pb, const size_t hf, const size_t id_s,
                                        const size_t ig_s, const size_t id, const size_t if_,
                                        const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 3.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.5 / p;
    const auto f_4 = alpha / p;
    const auto f_5 = 0.5 / p;

    auto *t_417 = buffer.data(target + 417);
    auto *t_418 = buffer.data(target + 418);
    auto *t_419 = buffer.data(target + 419);

    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *hf_209 = buffer.data(hf + 209);

    const auto *id_s_167 = buffer.data(id_s + 167);

    const auto *ig_s_417 = buffer.data(ig_s + 417);
    const auto *ig_s_418 = buffer.data(ig_s + 418);
    const auto *ig_s_419 = buffer.data(ig_s + 419);

    const auto *id_167 = buffer.data(id + 167);

    const auto *if__278 = buffer.data(if_ + 278);
    const auto *if__279 = buffer.data(if_ + 279);

#pragma omp simd aligned(t_417, t_418, t_419, pb_y, pb_z, hf_209, id_s_167, ig_s_417, \
                         ig_s_418, ig_s_419, id_167, if__278, if__279 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_417[k] = -f_4 * id_s_167[k]
                   + f_2 * ig_s_417[k]
                   + f_5 * id_167[k]
                   + pb_y[k] * if__278[k];

        t_418[k] = f_2 * ig_s_418[k]
                   + pb_y[k] * if__279[k];

        t_419[k] = f_0 * hf_209[k]
                   - f_1 * id_s_167[k]
                   + f_2 * ig_s_419[k]
                   + f_3 * id_167[k]
                   + pb_z[k] * if__279[k];
    }
}

auto
compute_prim_ig_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gg_s, const size_t gg,
                                 const size_t hf, const size_t hg, const size_t id_s,
                                 const size_t ig_s, const size_t id, const size_t if_,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_ig_kinetic_energy_0_piece0(buffer, target, pa, pb, gg_s, gg, hf, hg, id_s,
                                            ig_s, id, if_, ncols, alpha, beta, p);

    compute_prim_ig_kinetic_energy_0_piece1(buffer, target, pa, pb, gg_s, gg, hf, hg, id_s,
                                            ig_s, id, if_, ncols, alpha, beta, p);

    compute_prim_ig_kinetic_energy_0_piece2(buffer, target, pa, pb, gg_s, gg, hf, hg, id_s,
                                            ig_s, id, if_, ncols, alpha, beta, p);

    compute_prim_ig_kinetic_energy_0_piece3(buffer, target, pa, pb, gg_s, gg, hf, hg, id_s,
                                            ig_s, id, if_, ncols, alpha, beta, p);

    compute_prim_ig_kinetic_energy_0_piece4(buffer, target, pb, hf, id_s, ig_s, id, if_, ncols,
                                            alpha, beta, p);
}

}  // namespace simdkin
