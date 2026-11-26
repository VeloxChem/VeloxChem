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


#include "SimdKineticEnergyVrrRecIF.hpp"

#include "SimdAlign.hpp"

namespace simdkin {  // simdkin namespace

static auto
compute_prim_if_kinetic_energy_0_piece0(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gf_s,
                                        const size_t gf, const size_t hd, const size_t hf,
                                        const size_t ip_s, const size_t if_s, const size_t ip,
                                        const size_t id, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.5 / p;
    const auto f_6 = 4.0 * beta / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = alpha / p;
    const auto f_9 = beta / p;
    const auto f_10 = 3.0 * beta / p;
    const auto f_11 = 1.5 / p;
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

    const auto *pa_x = buffer.data(pa + 0);
    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf_s_0 = buffer.data(gf_s + 0);
    const auto *gf_s_10 = buffer.data(gf_s + 10);
    const auto *gf_s_16 = buffer.data(gf_s + 16);
    const auto *gf_s_20 = buffer.data(gf_s + 20);
    const auto *gf_s_29 = buffer.data(gf_s + 29);
    const auto *gf_s_30 = buffer.data(gf_s + 30);
    const auto *gf_s_36 = buffer.data(gf_s + 36);
    const auto *gf_s_59 = buffer.data(gf_s + 59);
    const auto *gf_s_66 = buffer.data(gf_s + 66);
    const auto *gf_s_79 = buffer.data(gf_s + 79);
    const auto *gf_s_86 = buffer.data(gf_s + 86);
    const auto *gf_s_99 = buffer.data(gf_s + 99);

    const auto *gf_0 = buffer.data(gf + 0);
    const auto *gf_10 = buffer.data(gf + 10);
    const auto *gf_16 = buffer.data(gf + 16);
    const auto *gf_20 = buffer.data(gf + 20);
    const auto *gf_29 = buffer.data(gf + 29);
    const auto *gf_30 = buffer.data(gf + 30);
    const auto *gf_36 = buffer.data(gf + 36);
    const auto *gf_59 = buffer.data(gf + 59);
    const auto *gf_66 = buffer.data(gf + 66);
    const auto *gf_79 = buffer.data(gf + 79);
    const auto *gf_86 = buffer.data(gf + 86);
    const auto *gf_99 = buffer.data(gf + 99);

    const auto *hd_0 = buffer.data(hd + 0);
    const auto *hd_3 = buffer.data(hd + 3);
    const auto *hd_5 = buffer.data(hd + 5);
    const auto *hd_6 = buffer.data(hd + 6);
    const auto *hd_9 = buffer.data(hd + 9);
    const auto *hd_11 = buffer.data(hd + 11);
    const auto *hd_12 = buffer.data(hd + 12);
    const auto *hd_17 = buffer.data(hd + 17);
    const auto *hd_18 = buffer.data(hd + 18);
    const auto *hd_21 = buffer.data(hd + 21);
    const auto *hd_23 = buffer.data(hd + 23);
    const auto *hd_27 = buffer.data(hd + 27);
    const auto *hd_28 = buffer.data(hd + 28);
    const auto *hd_29 = buffer.data(hd + 29);
    const auto *hd_30 = buffer.data(hd + 30);
    const auto *hd_33 = buffer.data(hd + 33);
    const auto *hd_35 = buffer.data(hd + 35);
    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_46 = buffer.data(hd + 46);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_52 = buffer.data(hd + 52);
    const auto *hd_57 = buffer.data(hd + 57);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_63 = buffer.data(hd + 63);

    const auto *hf_0 = buffer.data(hf + 0);
    const auto *hf_3 = buffer.data(hf + 3);
    const auto *hf_5 = buffer.data(hf + 5);
    const auto *hf_6 = buffer.data(hf + 6);
    const auto *hf_9 = buffer.data(hf + 9);
    const auto *hf_10 = buffer.data(hf + 10);
    const auto *hf_11 = buffer.data(hf + 11);
    const auto *hf_13 = buffer.data(hf + 13);
    const auto *hf_16 = buffer.data(hf + 16);
    const auto *hf_20 = buffer.data(hf + 20);
    const auto *hf_22 = buffer.data(hf + 22);
    const auto *hf_25 = buffer.data(hf + 25);
    const auto *hf_29 = buffer.data(hf + 29);
    const auto *hf_30 = buffer.data(hf + 30);
    const auto *hf_31 = buffer.data(hf + 31);
    const auto *hf_33 = buffer.data(hf + 33);
    const auto *hf_36 = buffer.data(hf + 36);
    const auto *hf_50 = buffer.data(hf + 50);
    const auto *hf_52 = buffer.data(hf + 52);
    const auto *hf_55 = buffer.data(hf + 55);
    const auto *hf_59 = buffer.data(hf + 59);
    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_79 = buffer.data(hf + 79);
    const auto *hf_86 = buffer.data(hf + 86);
    const auto *hf_99 = buffer.data(hf + 99);

    const auto *ip_s_0 = buffer.data(ip_s + 0);
    const auto *ip_s_1 = buffer.data(ip_s + 1);
    const auto *ip_s_2 = buffer.data(ip_s + 2);
    const auto *ip_s_8 = buffer.data(ip_s + 8);
    const auto *ip_s_11 = buffer.data(ip_s + 11);
    const auto *ip_s_16 = buffer.data(ip_s + 16);
    const auto *ip_s_17 = buffer.data(ip_s + 17);
    const auto *ip_s_20 = buffer.data(ip_s + 20);
    const auto *ip_s_28 = buffer.data(ip_s + 28);
    const auto *ip_s_29 = buffer.data(ip_s + 29);

    const auto *if_s_0 = buffer.data(if_s + 0);
    const auto *if_s_1 = buffer.data(if_s + 1);
    const auto *if_s_2 = buffer.data(if_s + 2);
    const auto *if_s_3 = buffer.data(if_s + 3);
    const auto *if_s_4 = buffer.data(if_s + 4);
    const auto *if_s_5 = buffer.data(if_s + 5);
    const auto *if_s_6 = buffer.data(if_s + 6);
    const auto *if_s_7 = buffer.data(if_s + 7);
    const auto *if_s_8 = buffer.data(if_s + 8);
    const auto *if_s_9 = buffer.data(if_s + 9);
    const auto *if_s_10 = buffer.data(if_s + 10);
    const auto *if_s_11 = buffer.data(if_s + 11);
    const auto *if_s_12 = buffer.data(if_s + 12);
    const auto *if_s_13 = buffer.data(if_s + 13);
    const auto *if_s_14 = buffer.data(if_s + 14);
    const auto *if_s_15 = buffer.data(if_s + 15);
    const auto *if_s_16 = buffer.data(if_s + 16);
    const auto *if_s_17 = buffer.data(if_s + 17);
    const auto *if_s_18 = buffer.data(if_s + 18);
    const auto *if_s_19 = buffer.data(if_s + 19);
    const auto *if_s_20 = buffer.data(if_s + 20);
    const auto *if_s_21 = buffer.data(if_s + 21);
    const auto *if_s_22 = buffer.data(if_s + 22);
    const auto *if_s_23 = buffer.data(if_s + 23);
    const auto *if_s_24 = buffer.data(if_s + 24);
    const auto *if_s_25 = buffer.data(if_s + 25);
    const auto *if_s_26 = buffer.data(if_s + 26);
    const auto *if_s_27 = buffer.data(if_s + 27);
    const auto *if_s_28 = buffer.data(if_s + 28);
    const auto *if_s_29 = buffer.data(if_s + 29);
    const auto *if_s_30 = buffer.data(if_s + 30);
    const auto *if_s_31 = buffer.data(if_s + 31);
    const auto *if_s_32 = buffer.data(if_s + 32);
    const auto *if_s_33 = buffer.data(if_s + 33);
    const auto *if_s_34 = buffer.data(if_s + 34);
    const auto *if_s_35 = buffer.data(if_s + 35);
    const auto *if_s_36 = buffer.data(if_s + 36);
    const auto *if_s_37 = buffer.data(if_s + 37);
    const auto *if_s_38 = buffer.data(if_s + 38);
    const auto *if_s_39 = buffer.data(if_s + 39);
    const auto *if_s_40 = buffer.data(if_s + 40);
    const auto *if_s_41 = buffer.data(if_s + 41);
    const auto *if_s_42 = buffer.data(if_s + 42);
    const auto *if_s_43 = buffer.data(if_s + 43);
    const auto *if_s_44 = buffer.data(if_s + 44);
    const auto *if_s_45 = buffer.data(if_s + 45);
    const auto *if_s_46 = buffer.data(if_s + 46);
    const auto *if_s_47 = buffer.data(if_s + 47);
    const auto *if_s_48 = buffer.data(if_s + 48);
    const auto *if_s_49 = buffer.data(if_s + 49);
    const auto *if_s_50 = buffer.data(if_s + 50);
    const auto *if_s_51 = buffer.data(if_s + 51);
    const auto *if_s_52 = buffer.data(if_s + 52);
    const auto *if_s_53 = buffer.data(if_s + 53);
    const auto *if_s_54 = buffer.data(if_s + 54);
    const auto *if_s_55 = buffer.data(if_s + 55);
    const auto *if_s_56 = buffer.data(if_s + 56);
    const auto *if_s_57 = buffer.data(if_s + 57);
    const auto *if_s_58 = buffer.data(if_s + 58);
    const auto *if_s_59 = buffer.data(if_s + 59);
    const auto *if_s_60 = buffer.data(if_s + 60);
    const auto *if_s_61 = buffer.data(if_s + 61);
    const auto *if_s_62 = buffer.data(if_s + 62);
    const auto *if_s_63 = buffer.data(if_s + 63);
    const auto *if_s_64 = buffer.data(if_s + 64);
    const auto *if_s_65 = buffer.data(if_s + 65);
    const auto *if_s_66 = buffer.data(if_s + 66);
    const auto *if_s_67 = buffer.data(if_s + 67);
    const auto *if_s_68 = buffer.data(if_s + 68);
    const auto *if_s_69 = buffer.data(if_s + 69);
    const auto *if_s_70 = buffer.data(if_s + 70);
    const auto *if_s_71 = buffer.data(if_s + 71);
    const auto *if_s_72 = buffer.data(if_s + 72);
    const auto *if_s_73 = buffer.data(if_s + 73);
    const auto *if_s_74 = buffer.data(if_s + 74);
    const auto *if_s_75 = buffer.data(if_s + 75);
    const auto *if_s_76 = buffer.data(if_s + 76);
    const auto *if_s_77 = buffer.data(if_s + 77);
    const auto *if_s_78 = buffer.data(if_s + 78);
    const auto *if_s_79 = buffer.data(if_s + 79);
    const auto *if_s_80 = buffer.data(if_s + 80);
    const auto *if_s_81 = buffer.data(if_s + 81);
    const auto *if_s_82 = buffer.data(if_s + 82);
    const auto *if_s_83 = buffer.data(if_s + 83);
    const auto *if_s_84 = buffer.data(if_s + 84);
    const auto *if_s_85 = buffer.data(if_s + 85);
    const auto *if_s_86 = buffer.data(if_s + 86);
    const auto *if_s_87 = buffer.data(if_s + 87);
    const auto *if_s_88 = buffer.data(if_s + 88);
    const auto *if_s_89 = buffer.data(if_s + 89);
    const auto *if_s_90 = buffer.data(if_s + 90);
    const auto *if_s_91 = buffer.data(if_s + 91);
    const auto *if_s_92 = buffer.data(if_s + 92);
    const auto *if_s_93 = buffer.data(if_s + 93);
    const auto *if_s_94 = buffer.data(if_s + 94);
    const auto *if_s_95 = buffer.data(if_s + 95);
    const auto *if_s_96 = buffer.data(if_s + 96);
    const auto *if_s_97 = buffer.data(if_s + 97);
    const auto *if_s_98 = buffer.data(if_s + 98);
    const auto *if_s_99 = buffer.data(if_s + 99);
    const auto *if_s_100 = buffer.data(if_s + 100);
    const auto *if_s_101 = buffer.data(if_s + 101);
    const auto *if_s_102 = buffer.data(if_s + 102);
    const auto *if_s_103 = buffer.data(if_s + 103);
    const auto *if_s_104 = buffer.data(if_s + 104);

    const auto *ip_0 = buffer.data(ip + 0);
    const auto *ip_1 = buffer.data(ip + 1);
    const auto *ip_2 = buffer.data(ip + 2);
    const auto *ip_8 = buffer.data(ip + 8);
    const auto *ip_11 = buffer.data(ip + 11);
    const auto *ip_16 = buffer.data(ip + 16);
    const auto *ip_17 = buffer.data(ip + 17);
    const auto *ip_20 = buffer.data(ip + 20);
    const auto *ip_28 = buffer.data(ip + 28);
    const auto *ip_29 = buffer.data(ip + 29);

    const auto *id_0 = buffer.data(id + 0);
    const auto *id_2 = buffer.data(id + 2);
    const auto *id_3 = buffer.data(id + 3);
    const auto *id_5 = buffer.data(id + 5);
    const auto *id_6 = buffer.data(id + 6);
    const auto *id_7 = buffer.data(id + 7);
    const auto *id_9 = buffer.data(id + 9);
    const auto *id_11 = buffer.data(id + 11);
    const auto *id_12 = buffer.data(id + 12);
    const auto *id_14 = buffer.data(id + 14);
    const auto *id_16 = buffer.data(id + 16);
    const auto *id_17 = buffer.data(id + 17);
    const auto *id_18 = buffer.data(id + 18);
    const auto *id_19 = buffer.data(id + 19);
    const auto *id_21 = buffer.data(id + 21);
    const auto *id_23 = buffer.data(id + 23);
    const auto *id_27 = buffer.data(id + 27);
    const auto *id_28 = buffer.data(id + 28);
    const auto *id_29 = buffer.data(id + 29);
    const auto *id_30 = buffer.data(id + 30);
    const auto *id_32 = buffer.data(id + 32);
    const auto *id_33 = buffer.data(id + 33);
    const auto *id_34 = buffer.data(id + 34);
    const auto *id_35 = buffer.data(id + 35);
    const auto *id_36 = buffer.data(id + 36);
    const auto *id_37 = buffer.data(id + 37);
    const auto *id_39 = buffer.data(id + 39);
    const auto *id_41 = buffer.data(id + 41);
    const auto *id_42 = buffer.data(id + 42);
    const auto *id_45 = buffer.data(id + 45);
    const auto *id_46 = buffer.data(id + 46);
    const auto *id_47 = buffer.data(id + 47);
    const auto *id_48 = buffer.data(id + 48);
    const auto *id_51 = buffer.data(id + 51);
    const auto *id_52 = buffer.data(id + 52);
    const auto *id_53 = buffer.data(id + 53);
    const auto *id_54 = buffer.data(id + 54);
    const auto *id_56 = buffer.data(id + 56);
    const auto *id_57 = buffer.data(id + 57);
    const auto *id_58 = buffer.data(id + 58);
    const auto *id_59 = buffer.data(id + 59);
    const auto *id_60 = buffer.data(id + 60);
    const auto *id_61 = buffer.data(id + 61);
    const auto *id_63 = buffer.data(id + 63);

#pragma omp simd aligned(t_0, t_1, t_2, pb_x, pb_y, pb_z, hd_0, ip_s_0, if_s_0, if_s_1, \
                         if_s_2, ip_0, id_0 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_0[k] = f_0 * hd_0[k]
                 - f_1 * ip_s_0[k]
                 + f_2 * if_s_0[k]
                 + f_3 * ip_0[k]
                 + pb_x[k] * id_0[k];

        t_1[k] = f_2 * if_s_1[k]
                 + pb_y[k] * id_0[k];

        t_2[k] = f_2 * if_s_2[k]
                 + pb_z[k] * id_0[k];
    }

#pragma omp simd aligned(t_3, t_4, t_5, pb_x, pb_y, hd_3, hd_5, if_s_3, if_s_4, if_s_5, id_2, \
                         id_3, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_3[k] = f_0 * hd_3[k]
                 + f_2 * if_s_3[k]
                 + pb_x[k] * id_3[k];

        t_4[k] = f_2 * if_s_4[k]
                 + pb_y[k] * id_2[k];

        t_5[k] = f_0 * hd_5[k]
                 + f_2 * if_s_5[k]
                 + pb_x[k] * id_5[k];
    }

#pragma omp simd aligned(t_6, t_7, t_8, t_9, pb_y, pb_z, ip_s_1, ip_s_2, if_s_6, if_s_7, \
                         if_s_8, if_s_9, ip_1, ip_2, id_3, id_5 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_6[k] = -f_1 * ip_s_1[k]
                 + f_2 * if_s_6[k]
                 + f_3 * ip_1[k]
                 + pb_y[k] * id_3[k];

        t_7[k] = f_2 * if_s_7[k]
                 + pb_z[k] * id_3[k];

        t_8[k] = f_2 * if_s_8[k]
                 + pb_y[k] * id_5[k];

        t_9[k] = -f_1 * ip_s_2[k]
                 + f_2 * if_s_9[k]
                 + f_3 * ip_2[k]
                 + pb_z[k] * id_5[k];
    }

#pragma omp simd aligned(t_10, t_11, t_12, pa_y, pb_y, pb_z, hd_0, hf_0, if_s_10, if_s_11, \
                         if_s_12, id_6 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_10[k] = pa_y[k] * hf_0[k]
                  + f_2 * if_s_10[k];

        t_11[k] = f_4 * hd_0[k]
                  + f_2 * if_s_11[k]
                  + pb_y[k] * id_6[k];

        t_12[k] = f_2 * if_s_12[k]
                  + pb_z[k] * id_6[k];
    }

#pragma omp simd aligned(t_13, t_14, t_15, pa_y, pb_x, pb_z, hd_9, hf_5, if_s_13, if_s_14, \
                         if_s_15, id_7, id_9 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_13[k] = f_5 * hd_9[k]
                  + f_2 * if_s_13[k]
                  + pb_x[k] * id_9[k];

        t_14[k] = f_2 * if_s_14[k]
                  + pb_z[k] * id_7[k];

        t_15[k] = pa_y[k] * hf_5[k]
                  + f_2 * if_s_15[k];
    }

#pragma omp simd aligned(t_16, t_17, t_18, pa_x, pb_y, pb_z, gf_s_16, gf_16, hd_5, hf_16, \
                         if_s_16, if_s_17, if_s_18, id_9, id_11 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_16[k] = -f_6 * gf_s_16[k]
                  + f_7 * gf_16[k]
                  + pa_x[k] * hf_16[k]
                  + f_2 * if_s_16[k];

        t_17[k] = f_2 * if_s_17[k]
                  + pb_z[k] * id_9[k];

        t_18[k] = f_4 * hd_5[k]
                  + f_2 * if_s_18[k]
                  + pb_y[k] * id_11[k];
    }

#pragma omp simd aligned(t_19, t_20, t_21, t_22, pa_y, pa_z, pb_y, pb_z, hd_0, hf_0, hf_9, \
                         if_s_19, if_s_20, if_s_21, if_s_22, id_12 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_19[k] = pa_y[k] * hf_9[k]
                  + f_2 * if_s_19[k];

        t_20[k] = pa_z[k] * hf_0[k]
                  + f_2 * if_s_20[k];

        t_21[k] = f_2 * if_s_21[k]
                  + pb_y[k] * id_12[k];

        t_22[k] = f_4 * hd_0[k]
                  + f_2 * if_s_22[k]
                  + pb_z[k] * id_12[k];
    }

#pragma omp simd aligned(t_23, t_24, t_25, t_26, pa_z, pb_x, pb_y, hd_17, hf_3, hf_6, if_s_23, \
                         if_s_24, if_s_25, if_s_26, id_14, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_23[k] = pa_z[k] * hf_3[k]
                  + f_2 * if_s_23[k];

        t_24[k] = f_2 * if_s_24[k]
                  + pb_y[k] * id_14[k];

        t_25[k] = f_5 * hd_17[k]
                  + f_2 * if_s_25[k]
                  + pb_x[k] * id_17[k];

        t_26[k] = pa_z[k] * hf_6[k]
                  + f_2 * if_s_26[k];
    }

#pragma omp simd aligned(t_27, t_28, t_29, pa_x, pb_y, gf_s_29, gf_29, hf_29, ip_s_8, if_s_27, \
                         if_s_28, if_s_29, ip_8, id_16, id_17 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_27[k] = -f_8 * ip_s_8[k]
                  + f_2 * if_s_27[k]
                  + f_4 * ip_8[k]
                  + pb_y[k] * id_16[k];

        t_28[k] = f_2 * if_s_28[k]
                  + pb_y[k] * id_17[k];

        t_29[k] = -f_6 * gf_s_29[k]
                  + f_7 * gf_29[k]
                  + pa_x[k] * hf_29[k]
                  + f_2 * if_s_29[k];
    }

#pragma omp simd aligned(t_30, t_31, t_32, pa_y, pb_y, pb_z, gf_s_0, gf_0, hd_6, hf_10, \
                         if_s_30, if_s_31, if_s_32, id_18 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_30[k] = -f_9 * gf_s_0[k]
                  + f_4 * gf_0[k]
                  + pa_y[k] * hf_10[k]
                  + f_2 * if_s_30[k];

        t_31[k] = f_3 * hd_6[k]
                  + f_2 * if_s_31[k]
                  + pb_y[k] * id_18[k];

        t_32[k] = f_2 * if_s_32[k]
                  + pb_z[k] * id_18[k];
    }

#pragma omp simd aligned(t_33, t_34, t_35, pb_x, pb_z, hd_21, hd_23, if_s_33, if_s_34, \
                         if_s_35, id_19, id_21, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_33[k] = f_7 * hd_21[k]
                  + f_2 * if_s_33[k]
                  + pb_x[k] * id_21[k];

        t_34[k] = f_2 * if_s_34[k]
                  + pb_z[k] * id_19[k];

        t_35[k] = f_7 * hd_23[k]
                  + f_2 * if_s_35[k]
                  + pb_x[k] * id_23[k];
    }

#pragma omp simd aligned(t_36, t_37, t_38, pa_x, pb_y, pb_z, gf_s_36, gf_36, hd_11, hf_36, \
                         if_s_36, if_s_37, if_s_38, id_21, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_36[k] = -f_10 * gf_s_36[k]
                  + f_11 * gf_36[k]
                  + pa_x[k] * hf_36[k]
                  + f_2 * if_s_36[k];

        t_37[k] = f_2 * if_s_37[k]
                  + pb_z[k] * id_21[k];

        t_38[k] = f_3 * hd_11[k]
                  + f_2 * if_s_38[k]
                  + pb_y[k] * id_23[k];
    }

#pragma omp simd aligned(t_39, t_40, t_41, pa_y, pa_z, pb_z, hf_11, hf_20, ip_s_11, if_s_39, \
                         if_s_40, if_s_41, ip_11, id_23 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_39[k] = -f_1 * ip_s_11[k]
                  + f_2 * if_s_39[k]
                  + f_3 * ip_11[k]
                  + pb_z[k] * id_23[k];

        t_40[k] = pa_y[k] * hf_20[k]
                  + f_2 * if_s_40[k];

        t_41[k] = pa_z[k] * hf_11[k]
                  + f_2 * if_s_41[k];
    }

#pragma omp simd aligned(t_42, t_43, t_44, t_45, pa_y, pa_z, pb_x, hd_28, hf_13, hf_22, hf_25, \
                         if_s_42, if_s_43, if_s_44, if_s_45, id_28 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_42[k] = pa_y[k] * hf_22[k]
                  + f_2 * if_s_42[k];

        t_43[k] = pa_z[k] * hf_13[k]
                  + f_2 * if_s_43[k];

        t_44[k] = f_7 * hd_28[k]
                  + f_2 * if_s_44[k]
                  + pb_x[k] * id_28[k];

        t_45[k] = pa_y[k] * hf_25[k]
                  + f_2 * if_s_45[k];
    }

#pragma omp simd aligned(t_46, t_47, t_48, pa_z, pb_y, pb_z, hd_9, hd_17, hf_16, if_s_46, \
                         if_s_47, if_s_48, id_27, id_29 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_46[k] = pa_z[k] * hf_16[k]
                  + f_2 * if_s_46[k];

        t_47[k] = f_4 * hd_9[k]
                  + f_2 * if_s_47[k]
                  + pb_z[k] * id_27[k];

        t_48[k] = f_4 * hd_17[k]
                  + f_2 * if_s_48[k]
                  + pb_y[k] * id_29[k];
    }

#pragma omp simd aligned(t_49, t_50, t_51, pa_y, pa_z, pb_y, gf_s_0, gf_0, hf_20, hf_29, \
                         if_s_49, if_s_50, if_s_51, id_30 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_49[k] = pa_y[k] * hf_29[k]
                  + f_2 * if_s_49[k];

        t_50[k] = -f_9 * gf_s_0[k]
                  + f_4 * gf_0[k]
                  + pa_z[k] * hf_20[k]
                  + f_2 * if_s_50[k];

        t_51[k] = f_2 * if_s_51[k]
                  + pb_y[k] * id_30[k];
    }

#pragma omp simd aligned(t_52, t_53, t_54, pb_x, pb_y, pb_z, hd_12, hd_33, if_s_52, if_s_53, \
                         if_s_54, id_30, id_32, id_33 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_52[k] = f_3 * hd_12[k]
                  + f_2 * if_s_52[k]
                  + pb_z[k] * id_30[k];

        t_53[k] = f_7 * hd_33[k]
                  + f_2 * if_s_53[k]
                  + pb_x[k] * id_33[k];

        t_54[k] = f_2 * if_s_54[k]
                  + pb_y[k] * id_32[k];
    }

#pragma omp simd aligned(t_55, t_56, t_57, pb_x, pb_y, hd_35, ip_s_16, ip_s_17, if_s_55, \
                         if_s_56, if_s_57, ip_16, ip_17, id_33, id_34, \
                         id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_55[k] = f_7 * hd_35[k]
                  + f_2 * if_s_55[k]
                  + pb_x[k] * id_35[k];

        t_56[k] = -f_1 * ip_s_16[k]
                  + f_2 * if_s_56[k]
                  + f_3 * ip_16[k]
                  + pb_y[k] * id_33[k];

        t_57[k] = -f_8 * ip_s_17[k]
                  + f_2 * if_s_57[k]
                  + f_4 * ip_17[k]
                  + pb_y[k] * id_34[k];
    }

#pragma omp simd aligned(t_58, t_59, t_60, pa_x, pa_y, pb_y, gf_s_10, gf_s_59, gf_10, gf_59, \
                         hf_30, hf_59, if_s_58, if_s_59, if_s_60, \
                         id_35 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_58[k] = f_2 * if_s_58[k]
                  + pb_y[k] * id_35[k];

        t_59[k] = -f_10 * gf_s_59[k]
                  + f_11 * gf_59[k]
                  + pa_x[k] * hf_59[k]
                  + f_2 * if_s_59[k];

        t_60[k] = -f_12 * gf_s_10[k]
                  + f_3 * gf_10[k]
                  + pa_y[k] * hf_30[k]
                  + f_2 * if_s_60[k];
    }

#pragma omp simd aligned(t_61, t_62, t_63, t_64, pb_x, pb_y, pb_z, hd_18, hd_39, if_s_61, \
                         if_s_62, if_s_63, if_s_64, id_36, id_37, \
                         id_39 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_61[k] = f_11 * hd_18[k]
                  + f_2 * if_s_61[k]
                  + pb_y[k] * id_36[k];

        t_62[k] = f_2 * if_s_62[k]
                  + pb_z[k] * id_36[k];

        t_63[k] = f_11 * hd_39[k]
                  + f_2 * if_s_63[k]
                  + pb_x[k] * id_39[k];

        t_64[k] = f_2 * if_s_64[k]
                  + pb_z[k] * id_37[k];
    }

#pragma omp simd aligned(t_65, t_66, t_67, pa_x, pb_x, pb_z, gf_s_66, gf_66, hd_41, hf_66, \
                         if_s_65, if_s_66, if_s_67, id_39, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_65[k] = f_11 * hd_41[k]
                  + f_2 * if_s_65[k]
                  + pb_x[k] * id_41[k];

        t_66[k] = -f_12 * gf_s_66[k]
                  + f_3 * gf_66[k]
                  + pa_x[k] * hf_66[k]
                  + f_2 * if_s_66[k];

        t_67[k] = f_2 * if_s_67[k]
                  + pb_z[k] * id_39[k];
    }

#pragma omp simd aligned(t_68, t_69, t_70, pa_z, pb_y, pb_z, hd_23, hf_30, ip_s_20, if_s_68, \
                         if_s_69, if_s_70, ip_20, id_41 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_68[k] = f_11 * hd_23[k]
                  + f_2 * if_s_68[k]
                  + pb_y[k] * id_41[k];

        t_69[k] = -f_1 * ip_s_20[k]
                  + f_2 * if_s_69[k]
                  + f_3 * ip_20[k]
                  + pb_z[k] * id_41[k];

        t_70[k] = pa_z[k] * hf_30[k]
                  + f_2 * if_s_70[k];
    }

#pragma omp simd aligned(t_71, t_72, t_73, pa_z, pb_z, hd_18, hf_31, hf_33, if_s_71, if_s_72, \
                         if_s_73, id_42 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_71[k] = pa_z[k] * hf_31[k]
                  + f_2 * if_s_71[k];

        t_72[k] = f_4 * hd_18[k]
                  + f_2 * if_s_72[k]
                  + pb_z[k] * id_42[k];

        t_73[k] = pa_z[k] * hf_33[k]
                  + f_2 * if_s_73[k];
    }

#pragma omp simd aligned(t_74, t_75, t_76, pa_z, pb_x, hd_46, hd_47, hf_36, if_s_74, if_s_75, \
                         if_s_76, id_46, id_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_74[k] = f_11 * hd_46[k]
                  + f_2 * if_s_74[k]
                  + pb_x[k] * id_46[k];

        t_75[k] = f_11 * hd_47[k]
                  + f_2 * if_s_75[k]
                  + pb_x[k] * id_47[k];

        t_76[k] = pa_z[k] * hf_36[k]
                  + f_2 * if_s_76[k];
    }

#pragma omp simd aligned(t_77, t_78, t_79, pa_x, pb_y, pb_z, gf_s_79, gf_79, hd_21, hd_29, \
                         hf_79, if_s_77, if_s_78, if_s_79, id_45, \
                         id_47 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_77[k] = f_4 * hd_21[k]
                  + f_2 * if_s_77[k]
                  + pb_z[k] * id_45[k];

        t_78[k] = f_3 * hd_29[k]
                  + f_2 * if_s_78[k]
                  + pb_y[k] * id_47[k];

        t_79[k] = -f_12 * gf_s_79[k]
                  + f_3 * gf_79[k]
                  + pa_x[k] * hf_79[k]
                  + f_2 * if_s_79[k];
    }

#pragma omp simd aligned(t_80, t_81, t_82, pa_y, pb_y, hd_30, hf_50, hf_52, if_s_80, if_s_81, \
                         if_s_82, id_48 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_80[k] = pa_y[k] * hf_50[k]
                  + f_2 * if_s_80[k];

        t_81[k] = f_4 * hd_30[k]
                  + f_2 * if_s_81[k]
                  + pb_y[k] * id_48[k];

        t_82[k] = pa_y[k] * hf_52[k]
                  + f_2 * if_s_82[k];
    }

#pragma omp simd aligned(t_83, t_84, t_85, pa_y, pb_x, hd_51, hd_52, hf_55, if_s_83, if_s_84, \
                         if_s_85, id_51, id_52 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_83[k] = f_11 * hd_51[k]
                  + f_2 * if_s_83[k]
                  + pb_x[k] * id_51[k];

        t_84[k] = f_11 * hd_52[k]
                  + f_2 * if_s_84[k]
                  + pb_x[k] * id_52[k];

        t_85[k] = pa_y[k] * hf_55[k]
                  + f_2 * if_s_85[k];
    }

#pragma omp simd aligned(t_86, t_87, t_88, pa_x, pb_y, pb_z, gf_s_86, gf_86, hd_27, hd_35, \
                         hf_86, if_s_86, if_s_87, if_s_88, id_51, \
                         id_53 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_86[k] = -f_12 * gf_s_86[k]
                  + f_3 * gf_86[k]
                  + pa_x[k] * hf_86[k]
                  + f_2 * if_s_86[k];

        t_87[k] = f_3 * hd_27[k]
                  + f_2 * if_s_87[k]
                  + pb_z[k] * id_51[k];

        t_88[k] = f_4 * hd_35[k]
                  + f_2 * if_s_88[k]
                  + pb_y[k] * id_53[k];
    }

#pragma omp simd aligned(t_89, t_90, t_91, pa_y, pa_z, pb_y, gf_s_20, gf_20, hf_50, hf_59, \
                         if_s_89, if_s_90, if_s_91, id_54 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_89[k] = pa_y[k] * hf_59[k]
                  + f_2 * if_s_89[k];

        t_90[k] = -f_12 * gf_s_20[k]
                  + f_3 * gf_20[k]
                  + pa_z[k] * hf_50[k]
                  + f_2 * if_s_90[k];

        t_91[k] = f_2 * if_s_91[k]
                  + pb_y[k] * id_54[k];
    }

#pragma omp simd aligned(t_92, t_93, t_94, pb_x, pb_y, pb_z, hd_30, hd_57, if_s_92, if_s_93, \
                         if_s_94, id_54, id_56, id_57 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_92[k] = f_11 * hd_30[k]
                  + f_2 * if_s_92[k]
                  + pb_z[k] * id_54[k];

        t_93[k] = f_11 * hd_57[k]
                  + f_2 * if_s_93[k]
                  + pb_x[k] * id_57[k];

        t_94[k] = f_2 * if_s_94[k]
                  + pb_y[k] * id_56[k];
    }

#pragma omp simd aligned(t_95, t_96, t_97, pb_x, pb_y, hd_59, ip_s_28, ip_s_29, if_s_95, \
                         if_s_96, if_s_97, ip_28, ip_29, id_57, id_58, \
                         id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_95[k] = f_11 * hd_59[k]
                  + f_2 * if_s_95[k]
                  + pb_x[k] * id_59[k];

        t_96[k] = -f_1 * ip_s_28[k]
                  + f_2 * if_s_96[k]
                  + f_3 * ip_28[k]
                  + pb_y[k] * id_57[k];

        t_97[k] = -f_8 * ip_s_29[k]
                  + f_2 * if_s_97[k]
                  + f_4 * ip_29[k]
                  + pb_y[k] * id_58[k];
    }

#pragma omp simd aligned(t_98, t_99, t_100, pa_x, pa_y, pb_y, gf_s_30, gf_s_99, gf_30, gf_99, \
                         hf_60, hf_99, if_s_98, if_s_99, if_s_100, \
                         id_59 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_98[k] = f_2 * if_s_98[k]
                  + pb_y[k] * id_59[k];

        t_99[k] = -f_12 * gf_s_99[k]
                  + f_3 * gf_99[k]
                  + pa_x[k] * hf_99[k]
                  + f_2 * if_s_99[k];

        t_100[k] = -f_10 * gf_s_30[k]
                   + f_11 * gf_30[k]
                   + pa_y[k] * hf_60[k]
                   + f_2 * if_s_100[k];
    }

#pragma omp simd aligned(t_101, t_102, t_103, t_104, pb_x, pb_y, pb_z, hd_36, hd_63, if_s_101, \
                         if_s_102, if_s_103, if_s_104, id_60, id_61, \
                         id_63 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_101[k] = f_7 * hd_36[k]
                   + f_2 * if_s_101[k]
                   + pb_y[k] * id_60[k];

        t_102[k] = f_2 * if_s_102[k]
                   + pb_z[k] * id_60[k];

        t_103[k] = f_3 * hd_63[k]
                   + f_2 * if_s_103[k]
                   + pb_x[k] * id_63[k];

        t_104[k] = f_2 * if_s_104[k]
                   + pb_z[k] * id_61[k];
    }
}

static auto
compute_prim_if_kinetic_energy_0_piece1(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gf_s,
                                        const size_t gf, const size_t hd, const size_t hf,
                                        const size_t ip_s, const size_t if_s, const size_t ip,
                                        const size_t id, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.5 / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = alpha / p;
    const auto f_9 = beta / p;
    const auto f_10 = 3.0 * beta / p;
    const auto f_11 = 1.5 / p;

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

    const auto *gf_s_50 = buffer.data(gf_s + 50);
    const auto *gf_s_106 = buffer.data(gf_s + 106);
    const auto *gf_s_119 = buffer.data(gf_s + 119);
    const auto *gf_s_126 = buffer.data(gf_s + 126);
    const auto *gf_s_129 = buffer.data(gf_s + 129);
    const auto *gf_s_136 = buffer.data(gf_s + 136);
    const auto *gf_s_149 = buffer.data(gf_s + 149);

    const auto *gf_50 = buffer.data(gf + 50);
    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_136 = buffer.data(gf + 136);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *hd_36 = buffer.data(hd + 36);
    const auto *hd_39 = buffer.data(hd + 39);
    const auto *hd_41 = buffer.data(hd + 41);
    const auto *hd_42 = buffer.data(hd + 42);
    const auto *hd_45 = buffer.data(hd + 45);
    const auto *hd_47 = buffer.data(hd + 47);
    const auto *hd_48 = buffer.data(hd + 48);
    const auto *hd_51 = buffer.data(hd + 51);
    const auto *hd_53 = buffer.data(hd + 53);
    const auto *hd_54 = buffer.data(hd + 54);
    const auto *hd_59 = buffer.data(hd + 59);
    const auto *hd_60 = buffer.data(hd + 60);
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
    const auto *hd_87 = buffer.data(hd + 87);
    const auto *hd_89 = buffer.data(hd + 89);
    const auto *hd_90 = buffer.data(hd + 90);
    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_100 = buffer.data(hd + 100);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_102 = buffer.data(hd + 102);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_106 = buffer.data(hd + 106);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_108 = buffer.data(hd + 108);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_112 = buffer.data(hd + 112);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_118 = buffer.data(hd + 118);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *hf_60 = buffer.data(hf + 60);
    const auto *hf_61 = buffer.data(hf + 61);
    const auto *hf_63 = buffer.data(hf + 63);
    const auto *hf_66 = buffer.data(hf + 66);
    const auto *hf_80 = buffer.data(hf + 80);
    const auto *hf_90 = buffer.data(hf + 90);
    const auto *hf_92 = buffer.data(hf + 92);
    const auto *hf_95 = buffer.data(hf + 95);
    const auto *hf_99 = buffer.data(hf + 99);
    const auto *hf_100 = buffer.data(hf + 100);
    const auto *hf_101 = buffer.data(hf + 101);
    const auto *hf_103 = buffer.data(hf + 103);
    const auto *hf_106 = buffer.data(hf + 106);
    const auto *hf_119 = buffer.data(hf + 119);
    const auto *hf_126 = buffer.data(hf + 126);
    const auto *hf_129 = buffer.data(hf + 129);
    const auto *hf_136 = buffer.data(hf + 136);
    const auto *hf_140 = buffer.data(hf + 140);
    const auto *hf_142 = buffer.data(hf + 142);
    const auto *hf_145 = buffer.data(hf + 145);
    const auto *hf_149 = buffer.data(hf + 149);
    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_158 = buffer.data(hf + 158);
    const auto *hf_159 = buffer.data(hf + 159);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_167 = buffer.data(hf + 167);
    const auto *hf_168 = buffer.data(hf + 168);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_170 = buffer.data(hf + 170);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_177 = buffer.data(hf + 177);
    const auto *hf_178 = buffer.data(hf + 178);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_180 = buffer.data(hf + 180);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_187 = buffer.data(hf + 187);
    const auto *hf_188 = buffer.data(hf + 188);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_196 = buffer.data(hf + 196);
    const auto *hf_197 = buffer.data(hf + 197);
    const auto *hf_198 = buffer.data(hf + 198);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_207 = buffer.data(hf + 207);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *ip_s_32 = buffer.data(ip_s + 32);
    const auto *ip_s_43 = buffer.data(ip_s + 43);
    const auto *ip_s_44 = buffer.data(ip_s + 44);
    const auto *ip_s_63 = buffer.data(ip_s + 63);
    const auto *ip_s_64 = buffer.data(ip_s + 64);

    const auto *if_s_105 = buffer.data(if_s + 105);
    const auto *if_s_106 = buffer.data(if_s + 106);
    const auto *if_s_107 = buffer.data(if_s + 107);
    const auto *if_s_108 = buffer.data(if_s + 108);
    const auto *if_s_109 = buffer.data(if_s + 109);
    const auto *if_s_110 = buffer.data(if_s + 110);
    const auto *if_s_111 = buffer.data(if_s + 111);
    const auto *if_s_112 = buffer.data(if_s + 112);
    const auto *if_s_113 = buffer.data(if_s + 113);
    const auto *if_s_114 = buffer.data(if_s + 114);
    const auto *if_s_115 = buffer.data(if_s + 115);
    const auto *if_s_116 = buffer.data(if_s + 116);
    const auto *if_s_117 = buffer.data(if_s + 117);
    const auto *if_s_118 = buffer.data(if_s + 118);
    const auto *if_s_119 = buffer.data(if_s + 119);
    const auto *if_s_120 = buffer.data(if_s + 120);
    const auto *if_s_121 = buffer.data(if_s + 121);
    const auto *if_s_122 = buffer.data(if_s + 122);
    const auto *if_s_123 = buffer.data(if_s + 123);
    const auto *if_s_124 = buffer.data(if_s + 124);
    const auto *if_s_125 = buffer.data(if_s + 125);
    const auto *if_s_126 = buffer.data(if_s + 126);
    const auto *if_s_127 = buffer.data(if_s + 127);
    const auto *if_s_128 = buffer.data(if_s + 128);
    const auto *if_s_129 = buffer.data(if_s + 129);
    const auto *if_s_130 = buffer.data(if_s + 130);
    const auto *if_s_131 = buffer.data(if_s + 131);
    const auto *if_s_132 = buffer.data(if_s + 132);
    const auto *if_s_133 = buffer.data(if_s + 133);
    const auto *if_s_134 = buffer.data(if_s + 134);
    const auto *if_s_135 = buffer.data(if_s + 135);
    const auto *if_s_136 = buffer.data(if_s + 136);
    const auto *if_s_137 = buffer.data(if_s + 137);
    const auto *if_s_138 = buffer.data(if_s + 138);
    const auto *if_s_139 = buffer.data(if_s + 139);
    const auto *if_s_140 = buffer.data(if_s + 140);
    const auto *if_s_141 = buffer.data(if_s + 141);
    const auto *if_s_142 = buffer.data(if_s + 142);
    const auto *if_s_143 = buffer.data(if_s + 143);
    const auto *if_s_144 = buffer.data(if_s + 144);
    const auto *if_s_145 = buffer.data(if_s + 145);
    const auto *if_s_146 = buffer.data(if_s + 146);
    const auto *if_s_147 = buffer.data(if_s + 147);
    const auto *if_s_148 = buffer.data(if_s + 148);
    const auto *if_s_149 = buffer.data(if_s + 149);
    const auto *if_s_150 = buffer.data(if_s + 150);
    const auto *if_s_151 = buffer.data(if_s + 151);
    const auto *if_s_152 = buffer.data(if_s + 152);
    const auto *if_s_153 = buffer.data(if_s + 153);
    const auto *if_s_154 = buffer.data(if_s + 154);
    const auto *if_s_155 = buffer.data(if_s + 155);
    const auto *if_s_156 = buffer.data(if_s + 156);
    const auto *if_s_157 = buffer.data(if_s + 157);
    const auto *if_s_158 = buffer.data(if_s + 158);
    const auto *if_s_159 = buffer.data(if_s + 159);
    const auto *if_s_160 = buffer.data(if_s + 160);
    const auto *if_s_161 = buffer.data(if_s + 161);
    const auto *if_s_162 = buffer.data(if_s + 162);
    const auto *if_s_163 = buffer.data(if_s + 163);
    const auto *if_s_164 = buffer.data(if_s + 164);
    const auto *if_s_165 = buffer.data(if_s + 165);
    const auto *if_s_166 = buffer.data(if_s + 166);
    const auto *if_s_167 = buffer.data(if_s + 167);
    const auto *if_s_168 = buffer.data(if_s + 168);
    const auto *if_s_169 = buffer.data(if_s + 169);
    const auto *if_s_170 = buffer.data(if_s + 170);
    const auto *if_s_171 = buffer.data(if_s + 171);
    const auto *if_s_172 = buffer.data(if_s + 172);
    const auto *if_s_173 = buffer.data(if_s + 173);
    const auto *if_s_174 = buffer.data(if_s + 174);
    const auto *if_s_175 = buffer.data(if_s + 175);
    const auto *if_s_176 = buffer.data(if_s + 176);
    const auto *if_s_177 = buffer.data(if_s + 177);
    const auto *if_s_178 = buffer.data(if_s + 178);
    const auto *if_s_179 = buffer.data(if_s + 179);
    const auto *if_s_180 = buffer.data(if_s + 180);
    const auto *if_s_181 = buffer.data(if_s + 181);
    const auto *if_s_182 = buffer.data(if_s + 182);
    const auto *if_s_183 = buffer.data(if_s + 183);
    const auto *if_s_184 = buffer.data(if_s + 184);
    const auto *if_s_185 = buffer.data(if_s + 185);
    const auto *if_s_186 = buffer.data(if_s + 186);
    const auto *if_s_187 = buffer.data(if_s + 187);
    const auto *if_s_188 = buffer.data(if_s + 188);
    const auto *if_s_189 = buffer.data(if_s + 189);
    const auto *if_s_190 = buffer.data(if_s + 190);
    const auto *if_s_191 = buffer.data(if_s + 191);
    const auto *if_s_192 = buffer.data(if_s + 192);
    const auto *if_s_193 = buffer.data(if_s + 193);
    const auto *if_s_194 = buffer.data(if_s + 194);
    const auto *if_s_195 = buffer.data(if_s + 195);
    const auto *if_s_196 = buffer.data(if_s + 196);
    const auto *if_s_197 = buffer.data(if_s + 197);
    const auto *if_s_198 = buffer.data(if_s + 198);
    const auto *if_s_199 = buffer.data(if_s + 199);
    const auto *if_s_200 = buffer.data(if_s + 200);
    const auto *if_s_201 = buffer.data(if_s + 201);
    const auto *if_s_202 = buffer.data(if_s + 202);
    const auto *if_s_203 = buffer.data(if_s + 203);
    const auto *if_s_204 = buffer.data(if_s + 204);
    const auto *if_s_205 = buffer.data(if_s + 205);
    const auto *if_s_206 = buffer.data(if_s + 206);
    const auto *if_s_207 = buffer.data(if_s + 207);
    const auto *if_s_208 = buffer.data(if_s + 208);
    const auto *if_s_209 = buffer.data(if_s + 209);
    const auto *if_s_210 = buffer.data(if_s + 210);
    const auto *if_s_211 = buffer.data(if_s + 211);
    const auto *if_s_212 = buffer.data(if_s + 212);
    const auto *if_s_213 = buffer.data(if_s + 213);
    const auto *if_s_214 = buffer.data(if_s + 214);
    const auto *if_s_215 = buffer.data(if_s + 215);
    const auto *if_s_216 = buffer.data(if_s + 216);

    const auto *ip_32 = buffer.data(ip + 32);
    const auto *ip_43 = buffer.data(ip + 43);
    const auto *ip_44 = buffer.data(ip + 44);
    const auto *ip_63 = buffer.data(ip + 63);
    const auto *ip_64 = buffer.data(ip + 64);

    const auto *id_63 = buffer.data(id + 63);
    const auto *id_65 = buffer.data(id + 65);
    const auto *id_66 = buffer.data(id + 66);
    const auto *id_69 = buffer.data(id + 69);
    const auto *id_70 = buffer.data(id + 70);
    const auto *id_71 = buffer.data(id + 71);
    const auto *id_72 = buffer.data(id + 72);
    const auto *id_75 = buffer.data(id + 75);
    const auto *id_76 = buffer.data(id + 76);
    const auto *id_77 = buffer.data(id + 77);
    const auto *id_78 = buffer.data(id + 78);
    const auto *id_81 = buffer.data(id + 81);
    const auto *id_82 = buffer.data(id + 82);
    const auto *id_83 = buffer.data(id + 83);
    const auto *id_84 = buffer.data(id + 84);
    const auto *id_86 = buffer.data(id + 86);
    const auto *id_87 = buffer.data(id + 87);
    const auto *id_88 = buffer.data(id + 88);
    const auto *id_89 = buffer.data(id + 89);
    const auto *id_90 = buffer.data(id + 90);
    const auto *id_91 = buffer.data(id + 91);
    const auto *id_93 = buffer.data(id + 93);
    const auto *id_95 = buffer.data(id + 95);
    const auto *id_96 = buffer.data(id + 96);
    const auto *id_100 = buffer.data(id + 100);
    const auto *id_101 = buffer.data(id + 101);
    const auto *id_102 = buffer.data(id + 102);
    const auto *id_105 = buffer.data(id + 105);
    const auto *id_106 = buffer.data(id + 106);
    const auto *id_107 = buffer.data(id + 107);
    const auto *id_108 = buffer.data(id + 108);
    const auto *id_111 = buffer.data(id + 111);
    const auto *id_112 = buffer.data(id + 112);
    const auto *id_113 = buffer.data(id + 113);
    const auto *id_114 = buffer.data(id + 114);
    const auto *id_117 = buffer.data(id + 117);
    const auto *id_118 = buffer.data(id + 118);
    const auto *id_120 = buffer.data(id + 120);
    const auto *id_122 = buffer.data(id + 122);
    const auto *id_123 = buffer.data(id + 123);
    const auto *id_125 = buffer.data(id + 125);
    const auto *id_126 = buffer.data(id + 126);
    const auto *id_127 = buffer.data(id + 127);
    const auto *id_129 = buffer.data(id + 129);
    const auto *id_130 = buffer.data(id + 130);
    const auto *id_131 = buffer.data(id + 131);

#pragma omp simd aligned(t_105, t_106, t_107, pa_x, pb_x, pb_z, gf_s_106, gf_106, hd_65, \
                         hf_106, if_s_105, if_s_106, if_s_107, id_63, \
                         id_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_105[k] = f_3 * hd_65[k]
                   + f_2 * if_s_105[k]
                   + pb_x[k] * id_65[k];

        t_106[k] = -f_9 * gf_s_106[k]
                   + f_4 * gf_106[k]
                   + pa_x[k] * hf_106[k]
                   + f_2 * if_s_106[k];

        t_107[k] = f_2 * if_s_107[k]
                   + pb_z[k] * id_63[k];
    }

#pragma omp simd aligned(t_108, t_109, t_110, pa_z, pb_y, pb_z, hd_41, hf_60, ip_s_32, \
                         if_s_108, if_s_109, if_s_110, ip_32, id_65 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_108[k] = f_7 * hd_41[k]
                   + f_2 * if_s_108[k]
                   + pb_y[k] * id_65[k];

        t_109[k] = -f_1 * ip_s_32[k]
                   + f_2 * if_s_109[k]
                   + f_3 * ip_32[k]
                   + pb_z[k] * id_65[k];

        t_110[k] = pa_z[k] * hf_60[k]
                   + f_2 * if_s_110[k];
    }

#pragma omp simd aligned(t_111, t_112, t_113, pa_z, pb_z, hd_36, hf_61, hf_63, if_s_111, \
                         if_s_112, if_s_113, id_66 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_111[k] = pa_z[k] * hf_61[k]
                   + f_2 * if_s_111[k];

        t_112[k] = f_4 * hd_36[k]
                   + f_2 * if_s_112[k]
                   + pb_z[k] * id_66[k];

        t_113[k] = pa_z[k] * hf_63[k]
                   + f_2 * if_s_113[k];
    }

#pragma omp simd aligned(t_114, t_115, t_116, pa_z, pb_x, hd_70, hd_71, hf_66, if_s_114, \
                         if_s_115, if_s_116, id_70, id_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_114[k] = f_3 * hd_70[k]
                   + f_2 * if_s_114[k]
                   + pb_x[k] * id_70[k];

        t_115[k] = f_3 * hd_71[k]
                   + f_2 * if_s_115[k]
                   + pb_x[k] * id_71[k];

        t_116[k] = pa_z[k] * hf_66[k]
                   + f_2 * if_s_116[k];
    }

#pragma omp simd aligned(t_117, t_118, t_119, pa_x, pb_y, pb_z, gf_s_119, gf_119, hd_39, \
                         hd_47, hf_119, if_s_117, if_s_118, if_s_119, id_69, \
                         id_71 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_117[k] = f_4 * hd_39[k]
                   + f_2 * if_s_117[k]
                   + pb_z[k] * id_69[k];

        t_118[k] = f_11 * hd_47[k]
                   + f_2 * if_s_118[k]
                   + pb_y[k] * id_71[k];

        t_119[k] = -f_9 * gf_s_119[k]
                   + f_4 * gf_119[k]
                   + pa_x[k] * hf_119[k]
                   + f_2 * if_s_119[k];
    }

#pragma omp simd aligned(t_120, t_121, t_122, pa_y, pb_y, pb_z, gf_s_50, gf_50, hd_42, hd_48, \
                         hf_80, if_s_120, if_s_121, if_s_122, id_72 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_120[k] = -f_9 * gf_s_50[k]
                   + f_4 * gf_50[k]
                   + pa_y[k] * hf_80[k]
                   + f_2 * if_s_120[k];

        t_121[k] = f_3 * hd_48[k]
                   + f_2 * if_s_121[k]
                   + pb_y[k] * id_72[k];

        t_122[k] = f_3 * hd_42[k]
                   + f_2 * if_s_122[k]
                   + pb_z[k] * id_72[k];
    }

#pragma omp simd aligned(t_123, t_124, t_125, pb_x, hd_75, hd_76, hd_77, if_s_123, if_s_124, \
                         if_s_125, id_75, id_76, id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_123[k] = f_3 * hd_75[k]
                   + f_2 * if_s_123[k]
                   + pb_x[k] * id_75[k];

        t_124[k] = f_3 * hd_76[k]
                   + f_2 * if_s_124[k]
                   + pb_x[k] * id_76[k];

        t_125[k] = f_3 * hd_77[k]
                   + f_2 * if_s_125[k]
                   + pb_x[k] * id_77[k];
    }

#pragma omp simd aligned(t_126, t_127, t_128, pa_x, pb_y, pb_z, gf_s_126, gf_126, hd_45, \
                         hd_53, hf_126, if_s_126, if_s_127, if_s_128, id_75, \
                         id_77 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_126[k] = -f_9 * gf_s_126[k]
                   + f_4 * gf_126[k]
                   + pa_x[k] * hf_126[k]
                   + f_2 * if_s_126[k];

        t_127[k] = f_3 * hd_45[k]
                   + f_2 * if_s_127[k]
                   + pb_z[k] * id_75[k];

        t_128[k] = f_3 * hd_53[k]
                   + f_2 * if_s_128[k]
                   + pb_y[k] * id_77[k];
    }

#pragma omp simd aligned(t_129, t_130, t_131, pa_x, pa_y, pb_y, gf_s_129, gf_129, hd_54, \
                         hf_90, hf_129, if_s_129, if_s_130, if_s_131, \
                         id_78 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_129[k] = -f_9 * gf_s_129[k]
                   + f_4 * gf_129[k]
                   + pa_x[k] * hf_129[k]
                   + f_2 * if_s_129[k];

        t_130[k] = pa_y[k] * hf_90[k]
                   + f_2 * if_s_130[k];

        t_131[k] = f_4 * hd_54[k]
                   + f_2 * if_s_131[k]
                   + pb_y[k] * id_78[k];
    }

#pragma omp simd aligned(t_132, t_133, t_134, t_135, pa_y, pb_x, hd_81, hd_82, hf_92, hf_95, \
                         if_s_132, if_s_133, if_s_134, if_s_135, id_81, \
                         id_82 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_132[k] = pa_y[k] * hf_92[k]
                   + f_2 * if_s_132[k];

        t_133[k] = f_3 * hd_81[k]
                   + f_2 * if_s_133[k]
                   + pb_x[k] * id_81[k];

        t_134[k] = f_3 * hd_82[k]
                   + f_2 * if_s_134[k]
                   + pb_x[k] * id_82[k];

        t_135[k] = pa_y[k] * hf_95[k]
                   + f_2 * if_s_135[k];
    }

#pragma omp simd aligned(t_136, t_137, t_138, pa_x, pb_y, pb_z, gf_s_136, gf_136, hd_51, \
                         hd_59, hf_136, if_s_136, if_s_137, if_s_138, id_81, \
                         id_83 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_136[k] = -f_9 * gf_s_136[k]
                   + f_4 * gf_136[k]
                   + pa_x[k] * hf_136[k]
                   + f_2 * if_s_136[k];

        t_137[k] = f_11 * hd_51[k]
                   + f_2 * if_s_137[k]
                   + pb_z[k] * id_81[k];

        t_138[k] = f_4 * hd_59[k]
                   + f_2 * if_s_138[k]
                   + pb_y[k] * id_83[k];
    }

#pragma omp simd aligned(t_139, t_140, t_141, pa_y, pa_z, pb_y, gf_s_50, gf_50, hf_90, hf_99, \
                         if_s_139, if_s_140, if_s_141, id_84 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_139[k] = pa_y[k] * hf_99[k]
                   + f_2 * if_s_139[k];

        t_140[k] = -f_10 * gf_s_50[k]
                   + f_11 * gf_50[k]
                   + pa_z[k] * hf_90[k]
                   + f_2 * if_s_140[k];

        t_141[k] = f_2 * if_s_141[k]
                   + pb_y[k] * id_84[k];
    }

#pragma omp simd aligned(t_142, t_143, t_144, pb_x, pb_y, pb_z, hd_54, hd_87, if_s_142, \
                         if_s_143, if_s_144, id_84, id_86, id_87 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_142[k] = f_7 * hd_54[k]
                   + f_2 * if_s_142[k]
                   + pb_z[k] * id_84[k];

        t_143[k] = f_3 * hd_87[k]
                   + f_2 * if_s_143[k]
                   + pb_x[k] * id_87[k];

        t_144[k] = f_2 * if_s_144[k]
                   + pb_y[k] * id_86[k];
    }

#pragma omp simd aligned(t_145, t_146, t_147, pb_x, pb_y, hd_89, ip_s_43, ip_s_44, if_s_145, \
                         if_s_146, if_s_147, ip_43, ip_44, id_87, id_88, \
                         id_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_145[k] = f_3 * hd_89[k]
                   + f_2 * if_s_145[k]
                   + pb_x[k] * id_89[k];

        t_146[k] = -f_1 * ip_s_43[k]
                   + f_2 * if_s_146[k]
                   + f_3 * ip_43[k]
                   + pb_y[k] * id_87[k];

        t_147[k] = -f_8 * ip_s_44[k]
                   + f_2 * if_s_147[k]
                   + f_4 * ip_44[k]
                   + pb_y[k] * id_88[k];
    }

#pragma omp simd aligned(t_148, t_149, t_150, pa_x, pb_y, gf_s_149, gf_149, hd_90, hf_149, \
                         hf_150, if_s_148, if_s_149, if_s_150, id_89 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_148[k] = f_2 * if_s_148[k]
                   + pb_y[k] * id_89[k];

        t_149[k] = -f_9 * gf_s_149[k]
                   + f_4 * gf_149[k]
                   + pa_x[k] * hf_149[k]
                   + f_2 * if_s_149[k];

        t_150[k] = f_11 * hd_90[k]
                   + pa_x[k] * hf_150[k]
                   + f_2 * if_s_150[k];
    }

#pragma omp simd aligned(t_151, t_152, t_153, t_154, pb_x, pb_y, pb_z, hd_60, hd_93, if_s_151, \
                         if_s_152, if_s_153, if_s_154, id_90, id_91, \
                         id_93 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_151[k] = f_5 * hd_60[k]
                   + f_2 * if_s_151[k]
                   + pb_y[k] * id_90[k];

        t_152[k] = f_2 * if_s_152[k]
                   + pb_z[k] * id_90[k];

        t_153[k] = f_4 * hd_93[k]
                   + f_2 * if_s_153[k]
                   + pb_x[k] * id_93[k];

        t_154[k] = f_2 * if_s_154[k]
                   + pb_z[k] * id_91[k];
    }

#pragma omp simd aligned(t_155, t_156, t_157, t_158, pa_x, pb_x, pb_z, hd_95, hf_156, hf_158, \
                         if_s_155, if_s_156, if_s_157, if_s_158, id_93, \
                         id_95 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_155[k] = f_4 * hd_95[k]
                   + f_2 * if_s_155[k]
                   + pb_x[k] * id_95[k];

        t_156[k] = pa_x[k] * hf_156[k]
                   + f_2 * if_s_156[k];

        t_157[k] = f_2 * if_s_157[k]
                   + pb_z[k] * id_93[k];

        t_158[k] = pa_x[k] * hf_158[k]
                   + f_2 * if_s_158[k];
    }

#pragma omp simd aligned(t_159, t_160, t_161, t_162, pa_x, pa_z, pb_z, hd_60, hf_100, hf_101, \
                         hf_159, if_s_159, if_s_160, if_s_161, if_s_162, \
                         id_96 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_159[k] = pa_x[k] * hf_159[k]
                   + f_2 * if_s_159[k];

        t_160[k] = pa_z[k] * hf_100[k]
                   + f_2 * if_s_160[k];

        t_161[k] = pa_z[k] * hf_101[k]
                   + f_2 * if_s_161[k];

        t_162[k] = f_4 * hd_60[k]
                   + f_2 * if_s_162[k]
                   + pb_z[k] * id_96[k];
    }

#pragma omp simd aligned(t_163, t_164, t_165, pa_z, pb_x, hd_100, hd_101, hf_103, if_s_163, \
                         if_s_164, if_s_165, id_100, id_101 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_163[k] = pa_z[k] * hf_103[k]
                   + f_2 * if_s_163[k];

        t_164[k] = f_4 * hd_100[k]
                   + f_2 * if_s_164[k]
                   + pb_x[k] * id_100[k];

        t_165[k] = f_4 * hd_101[k]
                   + f_2 * if_s_165[k]
                   + pb_x[k] * id_101[k];
    }

#pragma omp simd aligned(t_166, t_167, t_168, t_169, pa_x, hf_166, hf_167, hf_168, hf_169, \
                         if_s_166, if_s_167, if_s_168, if_s_169 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_166[k] = pa_x[k] * hf_166[k]
                   + f_2 * if_s_166[k];

        t_167[k] = pa_x[k] * hf_167[k]
                   + f_2 * if_s_167[k];

        t_168[k] = pa_x[k] * hf_168[k]
                   + f_2 * if_s_168[k];

        t_169[k] = pa_x[k] * hf_169[k]
                   + f_2 * if_s_169[k];
    }

#pragma omp simd aligned(t_170, t_171, t_172, pa_x, pb_y, pb_z, hd_66, hd_72, hd_102, hf_170, \
                         if_s_170, if_s_171, if_s_172, id_102 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_170[k] = f_11 * hd_102[k]
                   + pa_x[k] * hf_170[k]
                   + f_2 * if_s_170[k];

        t_171[k] = f_11 * hd_72[k]
                   + f_2 * if_s_171[k]
                   + pb_y[k] * id_102[k];

        t_172[k] = f_3 * hd_66[k]
                   + f_2 * if_s_172[k]
                   + pb_z[k] * id_102[k];
    }

#pragma omp simd aligned(t_173, t_174, t_175, pb_x, hd_105, hd_106, hd_107, if_s_173, \
                         if_s_174, if_s_175, id_105, id_106, id_107 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_173[k] = f_4 * hd_105[k]
                   + f_2 * if_s_173[k]
                   + pb_x[k] * id_105[k];

        t_174[k] = f_4 * hd_106[k]
                   + f_2 * if_s_174[k]
                   + pb_x[k] * id_106[k];

        t_175[k] = f_4 * hd_107[k]
                   + f_2 * if_s_175[k]
                   + pb_x[k] * id_107[k];
    }

#pragma omp simd aligned(t_176, t_177, t_178, t_179, pa_x, hf_176, hf_177, hf_178, hf_179, \
                         if_s_176, if_s_177, if_s_178, if_s_179 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_176[k] = pa_x[k] * hf_176[k]
                   + f_2 * if_s_176[k];

        t_177[k] = pa_x[k] * hf_177[k]
                   + f_2 * if_s_177[k];

        t_178[k] = pa_x[k] * hf_178[k]
                   + f_2 * if_s_178[k];

        t_179[k] = pa_x[k] * hf_179[k]
                   + f_2 * if_s_179[k];
    }

#pragma omp simd aligned(t_180, t_181, t_182, pa_x, pb_y, pb_z, hd_72, hd_78, hd_108, hf_180, \
                         if_s_180, if_s_181, if_s_182, id_108 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_180[k] = f_11 * hd_108[k]
                   + pa_x[k] * hf_180[k]
                   + f_2 * if_s_180[k];

        t_181[k] = f_3 * hd_78[k]
                   + f_2 * if_s_181[k]
                   + pb_y[k] * id_108[k];

        t_182[k] = f_11 * hd_72[k]
                   + f_2 * if_s_182[k]
                   + pb_z[k] * id_108[k];
    }

#pragma omp simd aligned(t_183, t_184, t_185, pb_x, hd_111, hd_112, hd_113, if_s_183, \
                         if_s_184, if_s_185, id_111, id_112, id_113 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_183[k] = f_4 * hd_111[k]
                   + f_2 * if_s_183[k]
                   + pb_x[k] * id_111[k];

        t_184[k] = f_4 * hd_112[k]
                   + f_2 * if_s_184[k]
                   + pb_x[k] * id_112[k];

        t_185[k] = f_4 * hd_113[k]
                   + f_2 * if_s_185[k]
                   + pb_x[k] * id_113[k];
    }

#pragma omp simd aligned(t_186, t_187, t_188, t_189, pa_x, hf_186, hf_187, hf_188, hf_189, \
                         if_s_186, if_s_187, if_s_188, if_s_189 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_186[k] = pa_x[k] * hf_186[k]
                   + f_2 * if_s_186[k];

        t_187[k] = pa_x[k] * hf_187[k]
                   + f_2 * if_s_187[k];

        t_188[k] = pa_x[k] * hf_188[k]
                   + f_2 * if_s_188[k];

        t_189[k] = pa_x[k] * hf_189[k]
                   + f_2 * if_s_189[k];
    }

#pragma omp simd aligned(t_190, t_191, t_192, pa_y, pb_y, hd_84, hf_140, hf_142, if_s_190, \
                         if_s_191, if_s_192, id_114 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_190[k] = pa_y[k] * hf_140[k]
                   + f_2 * if_s_190[k];

        t_191[k] = f_4 * hd_84[k]
                   + f_2 * if_s_191[k]
                   + pb_y[k] * id_114[k];

        t_192[k] = pa_y[k] * hf_142[k]
                   + f_2 * if_s_192[k];
    }

#pragma omp simd aligned(t_193, t_194, t_195, pa_y, pb_x, hd_117, hd_118, hf_145, if_s_193, \
                         if_s_194, if_s_195, id_117, id_118 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_193[k] = f_4 * hd_117[k]
                   + f_2 * if_s_193[k]
                   + pb_x[k] * id_117[k];

        t_194[k] = f_4 * hd_118[k]
                   + f_2 * if_s_194[k]
                   + pb_x[k] * id_118[k];

        t_195[k] = pa_y[k] * hf_145[k]
                   + f_2 * if_s_195[k];
    }

#pragma omp simd aligned(t_196, t_197, t_198, t_199, pa_x, hf_196, hf_197, hf_198, hf_199, \
                         if_s_196, if_s_197, if_s_198, if_s_199 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_196[k] = pa_x[k] * hf_196[k]
                   + f_2 * if_s_196[k];

        t_197[k] = pa_x[k] * hf_197[k]
                   + f_2 * if_s_197[k];

        t_198[k] = pa_x[k] * hf_198[k]
                   + f_2 * if_s_198[k];

        t_199[k] = pa_x[k] * hf_199[k]
                   + f_2 * if_s_199[k];
    }

#pragma omp simd aligned(t_200, t_201, t_202, pa_x, pb_y, pb_z, hd_84, hd_120, hf_200, \
                         if_s_200, if_s_201, if_s_202, id_120 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_200[k] = f_11 * hd_120[k]
                   + pa_x[k] * hf_200[k]
                   + f_2 * if_s_200[k];

        t_201[k] = f_2 * if_s_201[k]
                   + pb_y[k] * id_120[k];

        t_202[k] = f_5 * hd_84[k]
                   + f_2 * if_s_202[k]
                   + pb_z[k] * id_120[k];
    }

#pragma omp simd aligned(t_203, t_204, t_205, pb_x, pb_y, hd_123, hd_125, if_s_203, if_s_204, \
                         if_s_205, id_122, id_123, id_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_203[k] = f_4 * hd_123[k]
                   + f_2 * if_s_203[k]
                   + pb_x[k] * id_123[k];

        t_204[k] = f_2 * if_s_204[k]
                   + pb_y[k] * id_122[k];

        t_205[k] = f_4 * hd_125[k]
                   + f_2 * if_s_205[k]
                   + pb_x[k] * id_125[k];
    }

#pragma omp simd aligned(t_206, t_207, t_208, t_209, pa_x, pb_y, hf_206, hf_207, hf_209, \
                         if_s_206, if_s_207, if_s_208, if_s_209, \
                         id_125 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_206[k] = pa_x[k] * hf_206[k]
                   + f_2 * if_s_206[k];

        t_207[k] = pa_x[k] * hf_207[k]
                   + f_2 * if_s_207[k];

        t_208[k] = f_2 * if_s_208[k]
                   + pb_y[k] * id_125[k];

        t_209[k] = pa_x[k] * hf_209[k]
                   + f_2 * if_s_209[k];
    }

#pragma omp simd aligned(t_210, t_211, t_212, pb_x, pb_z, ip_s_63, ip_s_64, if_s_210, \
                         if_s_211, if_s_212, ip_63, ip_64, id_126, \
                         id_127 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_210[k] = -f_1 * ip_s_63[k]
                   + f_2 * if_s_210[k]
                   + f_3 * ip_63[k]
                   + pb_x[k] * id_126[k];

        t_211[k] = -f_8 * ip_s_64[k]
                   + f_2 * if_s_211[k]
                   + f_4 * ip_64[k]
                   + pb_x[k] * id_127[k];

        t_212[k] = f_2 * if_s_212[k]
                   + pb_z[k] * id_126[k];
    }

#pragma omp simd aligned(t_213, t_214, t_215, t_216, pb_x, pb_y, hd_93, ip_s_64, if_s_213, \
                         if_s_214, if_s_215, if_s_216, ip_64, id_129, id_130, \
                         id_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_213[k] = f_2 * if_s_213[k]
                   + pb_x[k] * id_129[k];

        t_214[k] = f_2 * if_s_214[k]
                   + pb_x[k] * id_130[k];

        t_215[k] = f_2 * if_s_215[k]
                   + pb_x[k] * id_131[k];

        t_216[k] = f_0 * hd_93[k]
                   - f_1 * ip_s_64[k]
                   + f_2 * if_s_216[k]
                   + f_3 * ip_64[k]
                   + pb_y[k] * id_129[k];
    }
}

static auto
compute_prim_if_kinetic_energy_0_piece2(CSimdMatrix &buffer, const size_t target,
                                        const size_t pa, const size_t pb, const size_t gf_s,
                                        const size_t gf, const size_t hd, const size_t hf,
                                        const size_t ip_s, const size_t if_s, const size_t ip,
                                        const size_t id, const size_t ncols, const double alpha,
                                        const double beta, const double p) -> void
{
    // NOTE: the factors are fixed by the pair of primitives, so they are formed
    // once rather than for every atom pair the pair reaches.

    const auto f_0 = 3.0 / p;
    const auto f_1 = 2.0 * alpha / p;
    const auto f_2 = 2.0 * alpha * beta / p;
    const auto f_3 = 1.0 / p;
    const auto f_4 = 0.5 / p;
    const auto f_5 = 2.5 / p;
    const auto f_6 = 4.0 * beta / p;
    const auto f_7 = 2.0 / p;
    const auto f_8 = alpha / p;
    const auto f_9 = beta / p;
    const auto f_10 = 3.0 * beta / p;
    const auto f_11 = 1.5 / p;
    const auto f_12 = 2.0 * beta / p;

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

    const auto *pa_y = buffer.data(pa + 1);
    const auto *pa_z = buffer.data(pa + 2);

    const auto *pb_x = buffer.data(pb + 0);
    const auto *pb_y = buffer.data(pb + 1);
    const auto *pb_z = buffer.data(pb + 2);

    const auto *gf_s_106 = buffer.data(gf_s + 106);
    const auto *gf_s_116 = buffer.data(gf_s + 116);
    const auto *gf_s_119 = buffer.data(gf_s + 119);
    const auto *gf_s_126 = buffer.data(gf_s + 126);
    const auto *gf_s_129 = buffer.data(gf_s + 129);
    const auto *gf_s_139 = buffer.data(gf_s + 139);
    const auto *gf_s_149 = buffer.data(gf_s + 149);

    const auto *gf_106 = buffer.data(gf + 106);
    const auto *gf_116 = buffer.data(gf + 116);
    const auto *gf_119 = buffer.data(gf + 119);
    const auto *gf_126 = buffer.data(gf + 126);
    const auto *gf_129 = buffer.data(gf + 129);
    const auto *gf_139 = buffer.data(gf + 139);
    const auto *gf_149 = buffer.data(gf + 149);

    const auto *hd_93 = buffer.data(hd + 93);
    const auto *hd_95 = buffer.data(hd + 95);
    const auto *hd_99 = buffer.data(hd + 99);
    const auto *hd_101 = buffer.data(hd + 101);
    const auto *hd_105 = buffer.data(hd + 105);
    const auto *hd_107 = buffer.data(hd + 107);
    const auto *hd_111 = buffer.data(hd + 111);
    const auto *hd_113 = buffer.data(hd + 113);
    const auto *hd_117 = buffer.data(hd + 117);
    const auto *hd_119 = buffer.data(hd + 119);
    const auto *hd_120 = buffer.data(hd + 120);
    const auto *hd_123 = buffer.data(hd + 123);
    const auto *hd_125 = buffer.data(hd + 125);

    const auto *hf_150 = buffer.data(hf + 150);
    const auto *hf_151 = buffer.data(hf + 151);
    const auto *hf_156 = buffer.data(hf + 156);
    const auto *hf_166 = buffer.data(hf + 166);
    const auto *hf_169 = buffer.data(hf + 169);
    const auto *hf_176 = buffer.data(hf + 176);
    const auto *hf_179 = buffer.data(hf + 179);
    const auto *hf_186 = buffer.data(hf + 186);
    const auto *hf_189 = buffer.data(hf + 189);
    const auto *hf_199 = buffer.data(hf + 199);
    const auto *hf_200 = buffer.data(hf + 200);
    const auto *hf_201 = buffer.data(hf + 201);
    const auto *hf_202 = buffer.data(hf + 202);
    const auto *hf_206 = buffer.data(hf + 206);
    const auto *hf_209 = buffer.data(hf + 209);

    const auto *ip_s_65 = buffer.data(ip_s + 65);
    const auto *ip_s_68 = buffer.data(ip_s + 68);
    const auto *ip_s_69 = buffer.data(ip_s + 69);
    const auto *ip_s_70 = buffer.data(ip_s + 70);
    const auto *ip_s_71 = buffer.data(ip_s + 71);
    const auto *ip_s_72 = buffer.data(ip_s + 72);
    const auto *ip_s_73 = buffer.data(ip_s + 73);
    const auto *ip_s_74 = buffer.data(ip_s + 74);
    const auto *ip_s_75 = buffer.data(ip_s + 75);
    const auto *ip_s_76 = buffer.data(ip_s + 76);
    const auto *ip_s_77 = buffer.data(ip_s + 77);
    const auto *ip_s_81 = buffer.data(ip_s + 81);
    const auto *ip_s_82 = buffer.data(ip_s + 82);
    const auto *ip_s_83 = buffer.data(ip_s + 83);

    const auto *if_s_217 = buffer.data(if_s + 217);
    const auto *if_s_218 = buffer.data(if_s + 218);
    const auto *if_s_219 = buffer.data(if_s + 219);
    const auto *if_s_220 = buffer.data(if_s + 220);
    const auto *if_s_221 = buffer.data(if_s + 221);
    const auto *if_s_222 = buffer.data(if_s + 222);
    const auto *if_s_223 = buffer.data(if_s + 223);
    const auto *if_s_224 = buffer.data(if_s + 224);
    const auto *if_s_225 = buffer.data(if_s + 225);
    const auto *if_s_226 = buffer.data(if_s + 226);
    const auto *if_s_227 = buffer.data(if_s + 227);
    const auto *if_s_228 = buffer.data(if_s + 228);
    const auto *if_s_229 = buffer.data(if_s + 229);
    const auto *if_s_230 = buffer.data(if_s + 230);
    const auto *if_s_231 = buffer.data(if_s + 231);
    const auto *if_s_232 = buffer.data(if_s + 232);
    const auto *if_s_233 = buffer.data(if_s + 233);
    const auto *if_s_234 = buffer.data(if_s + 234);
    const auto *if_s_235 = buffer.data(if_s + 235);
    const auto *if_s_236 = buffer.data(if_s + 236);
    const auto *if_s_237 = buffer.data(if_s + 237);
    const auto *if_s_238 = buffer.data(if_s + 238);
    const auto *if_s_239 = buffer.data(if_s + 239);
    const auto *if_s_240 = buffer.data(if_s + 240);
    const auto *if_s_241 = buffer.data(if_s + 241);
    const auto *if_s_242 = buffer.data(if_s + 242);
    const auto *if_s_243 = buffer.data(if_s + 243);
    const auto *if_s_244 = buffer.data(if_s + 244);
    const auto *if_s_245 = buffer.data(if_s + 245);
    const auto *if_s_246 = buffer.data(if_s + 246);
    const auto *if_s_247 = buffer.data(if_s + 247);
    const auto *if_s_248 = buffer.data(if_s + 248);
    const auto *if_s_249 = buffer.data(if_s + 249);
    const auto *if_s_250 = buffer.data(if_s + 250);
    const auto *if_s_251 = buffer.data(if_s + 251);
    const auto *if_s_252 = buffer.data(if_s + 252);
    const auto *if_s_253 = buffer.data(if_s + 253);
    const auto *if_s_254 = buffer.data(if_s + 254);
    const auto *if_s_255 = buffer.data(if_s + 255);
    const auto *if_s_256 = buffer.data(if_s + 256);
    const auto *if_s_257 = buffer.data(if_s + 257);
    const auto *if_s_258 = buffer.data(if_s + 258);
    const auto *if_s_259 = buffer.data(if_s + 259);
    const auto *if_s_260 = buffer.data(if_s + 260);
    const auto *if_s_261 = buffer.data(if_s + 261);
    const auto *if_s_262 = buffer.data(if_s + 262);
    const auto *if_s_263 = buffer.data(if_s + 263);
    const auto *if_s_264 = buffer.data(if_s + 264);
    const auto *if_s_265 = buffer.data(if_s + 265);
    const auto *if_s_266 = buffer.data(if_s + 266);
    const auto *if_s_267 = buffer.data(if_s + 267);
    const auto *if_s_268 = buffer.data(if_s + 268);
    const auto *if_s_269 = buffer.data(if_s + 269);
    const auto *if_s_270 = buffer.data(if_s + 270);
    const auto *if_s_271 = buffer.data(if_s + 271);
    const auto *if_s_272 = buffer.data(if_s + 272);
    const auto *if_s_273 = buffer.data(if_s + 273);
    const auto *if_s_274 = buffer.data(if_s + 274);
    const auto *if_s_275 = buffer.data(if_s + 275);
    const auto *if_s_276 = buffer.data(if_s + 276);
    const auto *if_s_277 = buffer.data(if_s + 277);
    const auto *if_s_278 = buffer.data(if_s + 278);
    const auto *if_s_279 = buffer.data(if_s + 279);

    const auto *ip_65 = buffer.data(ip + 65);
    const auto *ip_68 = buffer.data(ip + 68);
    const auto *ip_69 = buffer.data(ip + 69);
    const auto *ip_70 = buffer.data(ip + 70);
    const auto *ip_71 = buffer.data(ip + 71);
    const auto *ip_72 = buffer.data(ip + 72);
    const auto *ip_73 = buffer.data(ip + 73);
    const auto *ip_74 = buffer.data(ip + 74);
    const auto *ip_75 = buffer.data(ip + 75);
    const auto *ip_76 = buffer.data(ip + 76);
    const auto *ip_77 = buffer.data(ip + 77);
    const auto *ip_81 = buffer.data(ip + 81);
    const auto *ip_82 = buffer.data(ip + 82);
    const auto *ip_83 = buffer.data(ip + 83);

    const auto *id_129 = buffer.data(id + 129);
    const auto *id_131 = buffer.data(id + 131);
    const auto *id_134 = buffer.data(id + 134);
    const auto *id_135 = buffer.data(id + 135);
    const auto *id_136 = buffer.data(id + 136);
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
    const auto *id_159 = buffer.data(id + 159);
    const auto *id_160 = buffer.data(id + 160);
    const auto *id_161 = buffer.data(id + 161);
    const auto *id_162 = buffer.data(id + 162);
    const auto *id_164 = buffer.data(id + 164);
    const auto *id_165 = buffer.data(id + 165);
    const auto *id_166 = buffer.data(id + 166);
    const auto *id_167 = buffer.data(id + 167);

#pragma omp simd aligned(t_217, t_218, t_219, pb_y, pb_z, hd_95, ip_s_65, if_s_217, if_s_218, \
                         if_s_219, ip_65, id_129, id_131 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_217[k] = f_2 * if_s_217[k]
                   + pb_z[k] * id_129[k];

        t_218[k] = f_0 * hd_95[k]
                   + f_2 * if_s_218[k]
                   + pb_y[k] * id_131[k];

        t_219[k] = -f_1 * ip_s_65[k]
                   + f_2 * if_s_219[k]
                   + f_3 * ip_65[k]
                   + pb_z[k] * id_131[k];
    }

#pragma omp simd aligned(t_220, t_221, t_222, t_223, pa_z, pb_x, hf_150, hf_151, ip_s_68, \
                         if_s_220, if_s_221, if_s_222, if_s_223, ip_68, id_134, \
                         id_135 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_220[k] = pa_z[k] * hf_150[k]
                   + f_2 * if_s_220[k];

        t_221[k] = pa_z[k] * hf_151[k]
                   + f_2 * if_s_221[k];

        t_222[k] = -f_8 * ip_s_68[k]
                   + f_2 * if_s_222[k]
                   + f_4 * ip_68[k]
                   + pb_x[k] * id_134[k];

        t_223[k] = f_2 * if_s_223[k]
                   + pb_x[k] * id_135[k];
    }

#pragma omp simd aligned(t_224, t_225, t_226, t_227, pa_z, pb_x, pb_z, hd_93, hf_156, \
                         if_s_224, if_s_225, if_s_226, if_s_227, id_135, id_136, \
                         id_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_224[k] = f_2 * if_s_224[k]
                   + pb_x[k] * id_136[k];

        t_225[k] = f_2 * if_s_225[k]
                   + pb_x[k] * id_137[k];

        t_226[k] = pa_z[k] * hf_156[k]
                   + f_2 * if_s_226[k];

        t_227[k] = f_4 * hd_93[k]
                   + f_2 * if_s_227[k]
                   + pb_z[k] * id_135[k];
    }

#pragma omp simd aligned(t_228, t_229, pa_y, pb_y, gf_s_119, gf_119, hd_101, hf_169, if_s_228, \
                         if_s_229, id_137 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_228[k] = f_5 * hd_101[k]
                   + f_2 * if_s_228[k]
                   + pb_y[k] * id_137[k];

        t_229[k] = -f_6 * gf_s_119[k]
                   + f_7 * gf_119[k]
                   + pa_y[k] * hf_169[k]
                   + f_2 * if_s_229[k];
    }

#pragma omp simd aligned(t_230, t_231, t_232, pb_x, ip_s_69, ip_s_70, ip_s_71, if_s_230, \
                         if_s_231, if_s_232, ip_69, ip_70, ip_71, id_138, id_139, \
                         id_140 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_230[k] = -f_1 * ip_s_69[k]
                   + f_2 * if_s_230[k]
                   + f_3 * ip_69[k]
                   + pb_x[k] * id_138[k];

        t_231[k] = -f_8 * ip_s_70[k]
                   + f_2 * if_s_231[k]
                   + f_4 * ip_70[k]
                   + pb_x[k] * id_139[k];

        t_232[k] = -f_8 * ip_s_71[k]
                   + f_2 * if_s_232[k]
                   + f_4 * ip_71[k]
                   + pb_x[k] * id_140[k];
    }

#pragma omp simd aligned(t_233, t_234, t_235, t_236, pa_z, pb_x, gf_s_106, gf_106, hf_166, \
                         if_s_233, if_s_234, if_s_235, if_s_236, id_141, id_142, \
                         id_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_233[k] = f_2 * if_s_233[k]
                   + pb_x[k] * id_141[k];

        t_234[k] = f_2 * if_s_234[k]
                   + pb_x[k] * id_142[k];

        t_235[k] = f_2 * if_s_235[k]
                   + pb_x[k] * id_143[k];

        t_236[k] = -f_9 * gf_s_106[k]
                   + f_4 * gf_106[k]
                   + pa_z[k] * hf_166[k]
                   + f_2 * if_s_236[k];
    }

#pragma omp simd aligned(t_237, t_238, t_239, pa_y, pb_y, pb_z, gf_s_129, gf_129, hd_99, \
                         hd_107, hf_179, if_s_237, if_s_238, if_s_239, id_141, \
                         id_143 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_237[k] = f_3 * hd_99[k]
                   + f_2 * if_s_237[k]
                   + pb_z[k] * id_141[k];

        t_238[k] = f_7 * hd_107[k]
                   + f_2 * if_s_238[k]
                   + pb_y[k] * id_143[k];

        t_239[k] = -f_10 * gf_s_129[k]
                   + f_11 * gf_129[k]
                   + pa_y[k] * hf_179[k]
                   + f_2 * if_s_239[k];
    }

#pragma omp simd aligned(t_240, t_241, t_242, pb_x, ip_s_72, ip_s_73, ip_s_74, if_s_240, \
                         if_s_241, if_s_242, ip_72, ip_73, ip_74, id_144, id_145, \
                         id_146 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_240[k] = -f_1 * ip_s_72[k]
                   + f_2 * if_s_240[k]
                   + f_3 * ip_72[k]
                   + pb_x[k] * id_144[k];

        t_241[k] = -f_8 * ip_s_73[k]
                   + f_2 * if_s_241[k]
                   + f_4 * ip_73[k]
                   + pb_x[k] * id_145[k];

        t_242[k] = -f_8 * ip_s_74[k]
                   + f_2 * if_s_242[k]
                   + f_4 * ip_74[k]
                   + pb_x[k] * id_146[k];
    }

#pragma omp simd aligned(t_243, t_244, t_245, t_246, pa_z, pb_x, gf_s_116, gf_116, hf_176, \
                         if_s_243, if_s_244, if_s_245, if_s_246, id_147, id_148, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_243[k] = f_2 * if_s_243[k]
                   + pb_x[k] * id_147[k];

        t_244[k] = f_2 * if_s_244[k]
                   + pb_x[k] * id_148[k];

        t_245[k] = f_2 * if_s_245[k]
                   + pb_x[k] * id_149[k];

        t_246[k] = -f_12 * gf_s_116[k]
                   + f_3 * gf_116[k]
                   + pa_z[k] * hf_176[k]
                   + f_2 * if_s_246[k];
    }

#pragma omp simd aligned(t_247, t_248, t_249, pa_y, pb_y, pb_z, gf_s_139, gf_139, hd_105, \
                         hd_113, hf_189, if_s_247, if_s_248, if_s_249, id_147, \
                         id_149 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_247[k] = f_11 * hd_105[k]
                   + f_2 * if_s_247[k]
                   + pb_z[k] * id_147[k];

        t_248[k] = f_11 * hd_113[k]
                   + f_2 * if_s_248[k]
                   + pb_y[k] * id_149[k];

        t_249[k] = -f_12 * gf_s_139[k]
                   + f_3 * gf_139[k]
                   + pa_y[k] * hf_189[k]
                   + f_2 * if_s_249[k];
    }

#pragma omp simd aligned(t_250, t_251, t_252, pb_x, ip_s_75, ip_s_76, ip_s_77, if_s_250, \
                         if_s_251, if_s_252, ip_75, ip_76, ip_77, id_150, id_151, \
                         id_152 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_250[k] = -f_1 * ip_s_75[k]
                   + f_2 * if_s_250[k]
                   + f_3 * ip_75[k]
                   + pb_x[k] * id_150[k];

        t_251[k] = -f_8 * ip_s_76[k]
                   + f_2 * if_s_251[k]
                   + f_4 * ip_76[k]
                   + pb_x[k] * id_151[k];

        t_252[k] = -f_8 * ip_s_77[k]
                   + f_2 * if_s_252[k]
                   + f_4 * ip_77[k]
                   + pb_x[k] * id_152[k];
    }

#pragma omp simd aligned(t_253, t_254, t_255, t_256, pa_z, pb_x, gf_s_126, gf_126, hf_186, \
                         if_s_253, if_s_254, if_s_255, if_s_256, id_153, id_154, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_253[k] = f_2 * if_s_253[k]
                   + pb_x[k] * id_153[k];

        t_254[k] = f_2 * if_s_254[k]
                   + pb_x[k] * id_154[k];

        t_255[k] = f_2 * if_s_255[k]
                   + pb_x[k] * id_155[k];

        t_256[k] = -f_10 * gf_s_126[k]
                   + f_11 * gf_126[k]
                   + pa_z[k] * hf_186[k]
                   + f_2 * if_s_256[k];
    }

#pragma omp simd aligned(t_257, t_258, t_259, pa_y, pb_y, pb_z, gf_s_149, gf_149, hd_111, \
                         hd_119, hf_199, if_s_257, if_s_258, if_s_259, id_153, \
                         id_155 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_257[k] = f_7 * hd_111[k]
                   + f_2 * if_s_257[k]
                   + pb_z[k] * id_153[k];

        t_258[k] = f_3 * hd_119[k]
                   + f_2 * if_s_258[k]
                   + pb_y[k] * id_155[k];

        t_259[k] = -f_9 * gf_s_149[k]
                   + f_4 * gf_149[k]
                   + pa_y[k] * hf_199[k]
                   + f_2 * if_s_259[k];
    }

#pragma omp simd aligned(t_260, t_261, t_262, t_263, pa_y, pb_x, hd_120, hf_200, hf_201, \
                         hf_202, if_s_260, if_s_261, if_s_262, if_s_263, \
                         id_159 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_260[k] = pa_y[k] * hf_200[k]
                   + f_2 * if_s_260[k];

        t_261[k] = f_4 * hd_120[k]
                   + pa_y[k] * hf_201[k]
                   + f_2 * if_s_261[k];

        t_262[k] = pa_y[k] * hf_202[k]
                   + f_2 * if_s_262[k];

        t_263[k] = f_2 * if_s_263[k]
                   + pb_x[k] * id_159[k];
    }

#pragma omp simd aligned(t_264, t_265, t_266, pa_y, pb_x, hd_123, hf_206, if_s_264, if_s_265, \
                         if_s_266, id_160, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_264[k] = f_2 * if_s_264[k]
                   + pb_x[k] * id_160[k];

        t_265[k] = f_2 * if_s_265[k]
                   + pb_x[k] * id_161[k];

        t_266[k] = f_11 * hd_123[k]
                   + pa_y[k] * hf_206[k]
                   + f_2 * if_s_266[k];
    }

#pragma omp simd aligned(t_267, t_268, t_269, pa_y, pb_y, pb_z, hd_117, hd_125, hf_209, \
                         if_s_267, if_s_268, if_s_269, id_159, id_161 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_267[k] = f_5 * hd_117[k]
                   + f_2 * if_s_267[k]
                   + pb_z[k] * id_159[k];

        t_268[k] = f_4 * hd_125[k]
                   + f_2 * if_s_268[k]
                   + pb_y[k] * id_161[k];

        t_269[k] = pa_y[k] * hf_209[k]
                   + f_2 * if_s_269[k];
    }

#pragma omp simd aligned(t_270, t_271, t_272, pb_x, pb_y, ip_s_81, ip_s_83, if_s_270, \
                         if_s_271, if_s_272, ip_81, ip_83, id_162, \
                         id_164 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_270[k] = -f_1 * ip_s_81[k]
                   + f_2 * if_s_270[k]
                   + f_3 * ip_81[k]
                   + pb_x[k] * id_162[k];

        t_271[k] = f_2 * if_s_271[k]
                   + pb_y[k] * id_162[k];

        t_272[k] = -f_8 * ip_s_83[k]
                   + f_2 * if_s_272[k]
                   + f_4 * ip_83[k]
                   + pb_x[k] * id_164[k];
    }

#pragma omp simd aligned(t_273, t_274, t_275, t_276, pb_x, pb_y, ip_s_82, if_s_273, if_s_274, \
                         if_s_275, if_s_276, ip_82, id_165, id_166, \
                         id_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_273[k] = f_2 * if_s_273[k]
                   + pb_x[k] * id_165[k];

        t_274[k] = f_2 * if_s_274[k]
                   + pb_x[k] * id_166[k];

        t_275[k] = f_2 * if_s_275[k]
                   + pb_x[k] * id_167[k];

        t_276[k] = -f_1 * ip_s_82[k]
                   + f_2 * if_s_276[k]
                   + f_3 * ip_82[k]
                   + pb_y[k] * id_165[k];
    }

#pragma omp simd aligned(t_277, t_278, t_279, pb_y, pb_z, hd_125, ip_s_83, if_s_277, if_s_278, \
                         if_s_279, ip_83, id_166, id_167 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        t_277[k] = -f_8 * ip_s_83[k]
                   + f_2 * if_s_277[k]
                   + f_4 * ip_83[k]
                   + pb_y[k] * id_166[k];

        t_278[k] = f_2 * if_s_278[k]
                   + pb_y[k] * id_167[k];

        t_279[k] = f_0 * hd_125[k]
                   - f_1 * ip_s_83[k]
                   + f_2 * if_s_279[k]
                   + f_3 * ip_83[k]
                   + pb_z[k] * id_167[k];
    }
}

auto
compute_prim_if_kinetic_energy_0(CSimdMatrix &buffer, const size_t target, const size_t pa,
                                 const size_t pb, const size_t gf_s, const size_t gf,
                                 const size_t hd, const size_t hf, const size_t ip_s,
                                 const size_t if_s, const size_t ip, const size_t id,
                                 const size_t ncols, const double alpha, const double beta,
                                 const double p) -> void
{
    compute_prim_if_kinetic_energy_0_piece0(buffer, target, pa, pb, gf_s, gf, hd, hf, ip_s,
                                            if_s, ip, id, ncols, alpha, beta, p);

    compute_prim_if_kinetic_energy_0_piece1(buffer, target, pa, pb, gf_s, gf, hd, hf, ip_s,
                                            if_s, ip, id, ncols, alpha, beta, p);

    compute_prim_if_kinetic_energy_0_piece2(buffer, target, pa, pb, gf_s, gf, hd, hf, ip_s,
                                            if_s, ip, id, ncols, alpha, beta, p);
}

}  // namespace simdkin
