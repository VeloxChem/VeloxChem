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

#include "NuclearPotentialPrimRecDH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_dh(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_dh,
                               const size_t              idx_npot_0_sh,
                               const size_t              idx_npot_1_sh,
                               const size_t              idx_npot_0_pg,
                               const size_t              idx_npot_1_pg,
                               const size_t              idx_npot_0_ph,
                               const size_t              idx_npot_1_ph,
                               const CSimdArray<double>& factors,
                               const size_t              idx_rpa,
                               const size_t              idx_rpc,
                               const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_0 = pbuffer.data(idx_npot_0_sh);

    auto ta_0_xxxxy_0 = pbuffer.data(idx_npot_0_sh + 1);

    auto ta_0_xxxxz_0 = pbuffer.data(idx_npot_0_sh + 2);

    auto ta_0_xxxyy_0 = pbuffer.data(idx_npot_0_sh + 3);

    auto ta_0_xxxyz_0 = pbuffer.data(idx_npot_0_sh + 4);

    auto ta_0_xxxzz_0 = pbuffer.data(idx_npot_0_sh + 5);

    auto ta_0_xxyyy_0 = pbuffer.data(idx_npot_0_sh + 6);

    auto ta_0_xxyyz_0 = pbuffer.data(idx_npot_0_sh + 7);

    auto ta_0_xxyzz_0 = pbuffer.data(idx_npot_0_sh + 8);

    auto ta_0_xxzzz_0 = pbuffer.data(idx_npot_0_sh + 9);

    auto ta_0_xyyyy_0 = pbuffer.data(idx_npot_0_sh + 10);

    auto ta_0_xyyyz_0 = pbuffer.data(idx_npot_0_sh + 11);

    auto ta_0_xyyzz_0 = pbuffer.data(idx_npot_0_sh + 12);

    auto ta_0_xyzzz_0 = pbuffer.data(idx_npot_0_sh + 13);

    auto ta_0_xzzzz_0 = pbuffer.data(idx_npot_0_sh + 14);

    auto ta_0_yyyyy_0 = pbuffer.data(idx_npot_0_sh + 15);

    auto ta_0_yyyyz_0 = pbuffer.data(idx_npot_0_sh + 16);

    auto ta_0_yyyzz_0 = pbuffer.data(idx_npot_0_sh + 17);

    auto ta_0_yyzzz_0 = pbuffer.data(idx_npot_0_sh + 18);

    auto ta_0_yzzzz_0 = pbuffer.data(idx_npot_0_sh + 19);

    auto ta_0_zzzzz_0 = pbuffer.data(idx_npot_0_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_1 = pbuffer.data(idx_npot_1_sh);

    auto ta_0_xxxxy_1 = pbuffer.data(idx_npot_1_sh + 1);

    auto ta_0_xxxxz_1 = pbuffer.data(idx_npot_1_sh + 2);

    auto ta_0_xxxyy_1 = pbuffer.data(idx_npot_1_sh + 3);

    auto ta_0_xxxyz_1 = pbuffer.data(idx_npot_1_sh + 4);

    auto ta_0_xxxzz_1 = pbuffer.data(idx_npot_1_sh + 5);

    auto ta_0_xxyyy_1 = pbuffer.data(idx_npot_1_sh + 6);

    auto ta_0_xxyyz_1 = pbuffer.data(idx_npot_1_sh + 7);

    auto ta_0_xxyzz_1 = pbuffer.data(idx_npot_1_sh + 8);

    auto ta_0_xxzzz_1 = pbuffer.data(idx_npot_1_sh + 9);

    auto ta_0_xyyyy_1 = pbuffer.data(idx_npot_1_sh + 10);

    auto ta_0_xyyyz_1 = pbuffer.data(idx_npot_1_sh + 11);

    auto ta_0_xyyzz_1 = pbuffer.data(idx_npot_1_sh + 12);

    auto ta_0_xyzzz_1 = pbuffer.data(idx_npot_1_sh + 13);

    auto ta_0_xzzzz_1 = pbuffer.data(idx_npot_1_sh + 14);

    auto ta_0_yyyyy_1 = pbuffer.data(idx_npot_1_sh + 15);

    auto ta_0_yyyyz_1 = pbuffer.data(idx_npot_1_sh + 16);

    auto ta_0_yyyzz_1 = pbuffer.data(idx_npot_1_sh + 17);

    auto ta_0_yyzzz_1 = pbuffer.data(idx_npot_1_sh + 18);

    auto ta_0_yzzzz_1 = pbuffer.data(idx_npot_1_sh + 19);

    auto ta_0_zzzzz_1 = pbuffer.data(idx_npot_1_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_0 = pbuffer.data(idx_npot_0_pg);

    auto ta_x_xxxy_0 = pbuffer.data(idx_npot_0_pg + 1);

    auto ta_x_xxxz_0 = pbuffer.data(idx_npot_0_pg + 2);

    auto ta_x_xxyy_0 = pbuffer.data(idx_npot_0_pg + 3);

    auto ta_x_xxyz_0 = pbuffer.data(idx_npot_0_pg + 4);

    auto ta_x_xxzz_0 = pbuffer.data(idx_npot_0_pg + 5);

    auto ta_x_xyyy_0 = pbuffer.data(idx_npot_0_pg + 6);

    auto ta_x_xyyz_0 = pbuffer.data(idx_npot_0_pg + 7);

    auto ta_x_xyzz_0 = pbuffer.data(idx_npot_0_pg + 8);

    auto ta_x_xzzz_0 = pbuffer.data(idx_npot_0_pg + 9);

    auto ta_x_yyyy_0 = pbuffer.data(idx_npot_0_pg + 10);

    auto ta_x_yyyz_0 = pbuffer.data(idx_npot_0_pg + 11);

    auto ta_x_yyzz_0 = pbuffer.data(idx_npot_0_pg + 12);

    auto ta_x_yzzz_0 = pbuffer.data(idx_npot_0_pg + 13);

    auto ta_x_zzzz_0 = pbuffer.data(idx_npot_0_pg + 14);

    auto ta_y_xxxx_0 = pbuffer.data(idx_npot_0_pg + 15);

    auto ta_y_xxxy_0 = pbuffer.data(idx_npot_0_pg + 16);

    auto ta_y_xxxz_0 = pbuffer.data(idx_npot_0_pg + 17);

    auto ta_y_xxyy_0 = pbuffer.data(idx_npot_0_pg + 18);

    auto ta_y_xxyz_0 = pbuffer.data(idx_npot_0_pg + 19);

    auto ta_y_xxzz_0 = pbuffer.data(idx_npot_0_pg + 20);

    auto ta_y_xyyy_0 = pbuffer.data(idx_npot_0_pg + 21);

    auto ta_y_xyyz_0 = pbuffer.data(idx_npot_0_pg + 22);

    auto ta_y_xyzz_0 = pbuffer.data(idx_npot_0_pg + 23);

    auto ta_y_xzzz_0 = pbuffer.data(idx_npot_0_pg + 24);

    auto ta_y_yyyy_0 = pbuffer.data(idx_npot_0_pg + 25);

    auto ta_y_yyyz_0 = pbuffer.data(idx_npot_0_pg + 26);

    auto ta_y_yyzz_0 = pbuffer.data(idx_npot_0_pg + 27);

    auto ta_y_yzzz_0 = pbuffer.data(idx_npot_0_pg + 28);

    auto ta_y_zzzz_0 = pbuffer.data(idx_npot_0_pg + 29);

    auto ta_z_xxxx_0 = pbuffer.data(idx_npot_0_pg + 30);

    auto ta_z_xxxy_0 = pbuffer.data(idx_npot_0_pg + 31);

    auto ta_z_xxxz_0 = pbuffer.data(idx_npot_0_pg + 32);

    auto ta_z_xxyy_0 = pbuffer.data(idx_npot_0_pg + 33);

    auto ta_z_xxyz_0 = pbuffer.data(idx_npot_0_pg + 34);

    auto ta_z_xxzz_0 = pbuffer.data(idx_npot_0_pg + 35);

    auto ta_z_xyyy_0 = pbuffer.data(idx_npot_0_pg + 36);

    auto ta_z_xyyz_0 = pbuffer.data(idx_npot_0_pg + 37);

    auto ta_z_xyzz_0 = pbuffer.data(idx_npot_0_pg + 38);

    auto ta_z_xzzz_0 = pbuffer.data(idx_npot_0_pg + 39);

    auto ta_z_yyyy_0 = pbuffer.data(idx_npot_0_pg + 40);

    auto ta_z_yyyz_0 = pbuffer.data(idx_npot_0_pg + 41);

    auto ta_z_yyzz_0 = pbuffer.data(idx_npot_0_pg + 42);

    auto ta_z_yzzz_0 = pbuffer.data(idx_npot_0_pg + 43);

    auto ta_z_zzzz_0 = pbuffer.data(idx_npot_0_pg + 44);

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_1 = pbuffer.data(idx_npot_1_pg);

    auto ta_x_xxxy_1 = pbuffer.data(idx_npot_1_pg + 1);

    auto ta_x_xxxz_1 = pbuffer.data(idx_npot_1_pg + 2);

    auto ta_x_xxyy_1 = pbuffer.data(idx_npot_1_pg + 3);

    auto ta_x_xxyz_1 = pbuffer.data(idx_npot_1_pg + 4);

    auto ta_x_xxzz_1 = pbuffer.data(idx_npot_1_pg + 5);

    auto ta_x_xyyy_1 = pbuffer.data(idx_npot_1_pg + 6);

    auto ta_x_xyyz_1 = pbuffer.data(idx_npot_1_pg + 7);

    auto ta_x_xyzz_1 = pbuffer.data(idx_npot_1_pg + 8);

    auto ta_x_xzzz_1 = pbuffer.data(idx_npot_1_pg + 9);

    auto ta_x_yyyy_1 = pbuffer.data(idx_npot_1_pg + 10);

    auto ta_x_yyyz_1 = pbuffer.data(idx_npot_1_pg + 11);

    auto ta_x_yyzz_1 = pbuffer.data(idx_npot_1_pg + 12);

    auto ta_x_yzzz_1 = pbuffer.data(idx_npot_1_pg + 13);

    auto ta_x_zzzz_1 = pbuffer.data(idx_npot_1_pg + 14);

    auto ta_y_xxxx_1 = pbuffer.data(idx_npot_1_pg + 15);

    auto ta_y_xxxy_1 = pbuffer.data(idx_npot_1_pg + 16);

    auto ta_y_xxxz_1 = pbuffer.data(idx_npot_1_pg + 17);

    auto ta_y_xxyy_1 = pbuffer.data(idx_npot_1_pg + 18);

    auto ta_y_xxyz_1 = pbuffer.data(idx_npot_1_pg + 19);

    auto ta_y_xxzz_1 = pbuffer.data(idx_npot_1_pg + 20);

    auto ta_y_xyyy_1 = pbuffer.data(idx_npot_1_pg + 21);

    auto ta_y_xyyz_1 = pbuffer.data(idx_npot_1_pg + 22);

    auto ta_y_xyzz_1 = pbuffer.data(idx_npot_1_pg + 23);

    auto ta_y_xzzz_1 = pbuffer.data(idx_npot_1_pg + 24);

    auto ta_y_yyyy_1 = pbuffer.data(idx_npot_1_pg + 25);

    auto ta_y_yyyz_1 = pbuffer.data(idx_npot_1_pg + 26);

    auto ta_y_yyzz_1 = pbuffer.data(idx_npot_1_pg + 27);

    auto ta_y_yzzz_1 = pbuffer.data(idx_npot_1_pg + 28);

    auto ta_y_zzzz_1 = pbuffer.data(idx_npot_1_pg + 29);

    auto ta_z_xxxx_1 = pbuffer.data(idx_npot_1_pg + 30);

    auto ta_z_xxxy_1 = pbuffer.data(idx_npot_1_pg + 31);

    auto ta_z_xxxz_1 = pbuffer.data(idx_npot_1_pg + 32);

    auto ta_z_xxyy_1 = pbuffer.data(idx_npot_1_pg + 33);

    auto ta_z_xxyz_1 = pbuffer.data(idx_npot_1_pg + 34);

    auto ta_z_xxzz_1 = pbuffer.data(idx_npot_1_pg + 35);

    auto ta_z_xyyy_1 = pbuffer.data(idx_npot_1_pg + 36);

    auto ta_z_xyyz_1 = pbuffer.data(idx_npot_1_pg + 37);

    auto ta_z_xyzz_1 = pbuffer.data(idx_npot_1_pg + 38);

    auto ta_z_xzzz_1 = pbuffer.data(idx_npot_1_pg + 39);

    auto ta_z_yyyy_1 = pbuffer.data(idx_npot_1_pg + 40);

    auto ta_z_yyyz_1 = pbuffer.data(idx_npot_1_pg + 41);

    auto ta_z_yyzz_1 = pbuffer.data(idx_npot_1_pg + 42);

    auto ta_z_yzzz_1 = pbuffer.data(idx_npot_1_pg + 43);

    auto ta_z_zzzz_1 = pbuffer.data(idx_npot_1_pg + 44);

    // Set up components of auxiliary buffer : PH

    auto ta_x_xxxxx_0 = pbuffer.data(idx_npot_0_ph);

    auto ta_x_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 1);

    auto ta_x_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 2);

    auto ta_x_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 3);

    auto ta_x_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 4);

    auto ta_x_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 5);

    auto ta_x_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 6);

    auto ta_x_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 7);

    auto ta_x_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 8);

    auto ta_x_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 9);

    auto ta_x_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 10);

    auto ta_x_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 11);

    auto ta_x_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 12);

    auto ta_x_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 13);

    auto ta_x_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 14);

    auto ta_x_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 15);

    auto ta_x_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 16);

    auto ta_x_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 17);

    auto ta_x_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 18);

    auto ta_x_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 19);

    auto ta_x_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 20);

    auto ta_y_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 21);

    auto ta_y_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 22);

    auto ta_y_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 23);

    auto ta_y_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 24);

    auto ta_y_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 25);

    auto ta_y_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 26);

    auto ta_y_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 27);

    auto ta_y_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 28);

    auto ta_y_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 29);

    auto ta_y_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 30);

    auto ta_y_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 31);

    auto ta_y_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 32);

    auto ta_y_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 33);

    auto ta_y_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 34);

    auto ta_y_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 35);

    auto ta_y_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 36);

    auto ta_y_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 37);

    auto ta_y_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 38);

    auto ta_y_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 39);

    auto ta_y_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 40);

    auto ta_y_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 41);

    auto ta_z_xxxxx_0 = pbuffer.data(idx_npot_0_ph + 42);

    auto ta_z_xxxxy_0 = pbuffer.data(idx_npot_0_ph + 43);

    auto ta_z_xxxxz_0 = pbuffer.data(idx_npot_0_ph + 44);

    auto ta_z_xxxyy_0 = pbuffer.data(idx_npot_0_ph + 45);

    auto ta_z_xxxyz_0 = pbuffer.data(idx_npot_0_ph + 46);

    auto ta_z_xxxzz_0 = pbuffer.data(idx_npot_0_ph + 47);

    auto ta_z_xxyyy_0 = pbuffer.data(idx_npot_0_ph + 48);

    auto ta_z_xxyyz_0 = pbuffer.data(idx_npot_0_ph + 49);

    auto ta_z_xxyzz_0 = pbuffer.data(idx_npot_0_ph + 50);

    auto ta_z_xxzzz_0 = pbuffer.data(idx_npot_0_ph + 51);

    auto ta_z_xyyyy_0 = pbuffer.data(idx_npot_0_ph + 52);

    auto ta_z_xyyyz_0 = pbuffer.data(idx_npot_0_ph + 53);

    auto ta_z_xyyzz_0 = pbuffer.data(idx_npot_0_ph + 54);

    auto ta_z_xyzzz_0 = pbuffer.data(idx_npot_0_ph + 55);

    auto ta_z_xzzzz_0 = pbuffer.data(idx_npot_0_ph + 56);

    auto ta_z_yyyyy_0 = pbuffer.data(idx_npot_0_ph + 57);

    auto ta_z_yyyyz_0 = pbuffer.data(idx_npot_0_ph + 58);

    auto ta_z_yyyzz_0 = pbuffer.data(idx_npot_0_ph + 59);

    auto ta_z_yyzzz_0 = pbuffer.data(idx_npot_0_ph + 60);

    auto ta_z_yzzzz_0 = pbuffer.data(idx_npot_0_ph + 61);

    auto ta_z_zzzzz_0 = pbuffer.data(idx_npot_0_ph + 62);

    // Set up components of auxiliary buffer : PH

    auto ta_x_xxxxx_1 = pbuffer.data(idx_npot_1_ph);

    auto ta_x_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 1);

    auto ta_x_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 2);

    auto ta_x_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 3);

    auto ta_x_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 4);

    auto ta_x_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 5);

    auto ta_x_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 6);

    auto ta_x_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 7);

    auto ta_x_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 8);

    auto ta_x_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 9);

    auto ta_x_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 10);

    auto ta_x_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 11);

    auto ta_x_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 12);

    auto ta_x_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 13);

    auto ta_x_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 14);

    auto ta_x_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 15);

    auto ta_x_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 16);

    auto ta_x_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 17);

    auto ta_x_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 18);

    auto ta_x_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 19);

    auto ta_x_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 20);

    auto ta_y_xxxxx_1 = pbuffer.data(idx_npot_1_ph + 21);

    auto ta_y_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 22);

    auto ta_y_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 23);

    auto ta_y_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 24);

    auto ta_y_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 25);

    auto ta_y_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 26);

    auto ta_y_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 27);

    auto ta_y_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 28);

    auto ta_y_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 29);

    auto ta_y_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 30);

    auto ta_y_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 31);

    auto ta_y_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 32);

    auto ta_y_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 33);

    auto ta_y_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 34);

    auto ta_y_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 35);

    auto ta_y_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 36);

    auto ta_y_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 37);

    auto ta_y_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 38);

    auto ta_y_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 39);

    auto ta_y_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 40);

    auto ta_y_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 41);

    auto ta_z_xxxxx_1 = pbuffer.data(idx_npot_1_ph + 42);

    auto ta_z_xxxxy_1 = pbuffer.data(idx_npot_1_ph + 43);

    auto ta_z_xxxxz_1 = pbuffer.data(idx_npot_1_ph + 44);

    auto ta_z_xxxyy_1 = pbuffer.data(idx_npot_1_ph + 45);

    auto ta_z_xxxyz_1 = pbuffer.data(idx_npot_1_ph + 46);

    auto ta_z_xxxzz_1 = pbuffer.data(idx_npot_1_ph + 47);

    auto ta_z_xxyyy_1 = pbuffer.data(idx_npot_1_ph + 48);

    auto ta_z_xxyyz_1 = pbuffer.data(idx_npot_1_ph + 49);

    auto ta_z_xxyzz_1 = pbuffer.data(idx_npot_1_ph + 50);

    auto ta_z_xxzzz_1 = pbuffer.data(idx_npot_1_ph + 51);

    auto ta_z_xyyyy_1 = pbuffer.data(idx_npot_1_ph + 52);

    auto ta_z_xyyyz_1 = pbuffer.data(idx_npot_1_ph + 53);

    auto ta_z_xyyzz_1 = pbuffer.data(idx_npot_1_ph + 54);

    auto ta_z_xyzzz_1 = pbuffer.data(idx_npot_1_ph + 55);

    auto ta_z_xzzzz_1 = pbuffer.data(idx_npot_1_ph + 56);

    auto ta_z_yyyyy_1 = pbuffer.data(idx_npot_1_ph + 57);

    auto ta_z_yyyyz_1 = pbuffer.data(idx_npot_1_ph + 58);

    auto ta_z_yyyzz_1 = pbuffer.data(idx_npot_1_ph + 59);

    auto ta_z_yyzzz_1 = pbuffer.data(idx_npot_1_ph + 60);

    auto ta_z_yzzzz_1 = pbuffer.data(idx_npot_1_ph + 61);

    auto ta_z_zzzzz_1 = pbuffer.data(idx_npot_1_ph + 62);

    // Set up 0-21 components of targeted buffer : DH

    auto ta_xx_xxxxx_0 = pbuffer.data(idx_npot_0_dh);

    auto ta_xx_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 1);

    auto ta_xx_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 2);

    auto ta_xx_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 3);

    auto ta_xx_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 4);

    auto ta_xx_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 5);

    auto ta_xx_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 6);

    auto ta_xx_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 7);

    auto ta_xx_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 8);

    auto ta_xx_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 9);

    auto ta_xx_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 10);

    auto ta_xx_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 11);

    auto ta_xx_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 12);

    auto ta_xx_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 13);

    auto ta_xx_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 14);

    auto ta_xx_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 15);

    auto ta_xx_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 16);

    auto ta_xx_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 17);

    auto ta_xx_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 18);

    auto ta_xx_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 19);

    auto ta_xx_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 20);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_0_xxxxx_0,  \
                             ta_0_xxxxx_1,  \
                             ta_0_xxxxy_0,  \
                             ta_0_xxxxy_1,  \
                             ta_0_xxxxz_0,  \
                             ta_0_xxxxz_1,  \
                             ta_0_xxxyy_0,  \
                             ta_0_xxxyy_1,  \
                             ta_0_xxxyz_0,  \
                             ta_0_xxxyz_1,  \
                             ta_0_xxxzz_0,  \
                             ta_0_xxxzz_1,  \
                             ta_0_xxyyy_0,  \
                             ta_0_xxyyy_1,  \
                             ta_0_xxyyz_0,  \
                             ta_0_xxyyz_1,  \
                             ta_0_xxyzz_0,  \
                             ta_0_xxyzz_1,  \
                             ta_0_xxzzz_0,  \
                             ta_0_xxzzz_1,  \
                             ta_0_xyyyy_0,  \
                             ta_0_xyyyy_1,  \
                             ta_0_xyyyz_0,  \
                             ta_0_xyyyz_1,  \
                             ta_0_xyyzz_0,  \
                             ta_0_xyyzz_1,  \
                             ta_0_xyzzz_0,  \
                             ta_0_xyzzz_1,  \
                             ta_0_xzzzz_0,  \
                             ta_0_xzzzz_1,  \
                             ta_0_yyyyy_0,  \
                             ta_0_yyyyy_1,  \
                             ta_0_yyyyz_0,  \
                             ta_0_yyyyz_1,  \
                             ta_0_yyyzz_0,  \
                             ta_0_yyyzz_1,  \
                             ta_0_yyzzz_0,  \
                             ta_0_yyzzz_1,  \
                             ta_0_yzzzz_0,  \
                             ta_0_yzzzz_1,  \
                             ta_0_zzzzz_0,  \
                             ta_0_zzzzz_1,  \
                             ta_x_xxxx_0,   \
                             ta_x_xxxx_1,   \
                             ta_x_xxxxx_0,  \
                             ta_x_xxxxx_1,  \
                             ta_x_xxxxy_0,  \
                             ta_x_xxxxy_1,  \
                             ta_x_xxxxz_0,  \
                             ta_x_xxxxz_1,  \
                             ta_x_xxxy_0,   \
                             ta_x_xxxy_1,   \
                             ta_x_xxxyy_0,  \
                             ta_x_xxxyy_1,  \
                             ta_x_xxxyz_0,  \
                             ta_x_xxxyz_1,  \
                             ta_x_xxxz_0,   \
                             ta_x_xxxz_1,   \
                             ta_x_xxxzz_0,  \
                             ta_x_xxxzz_1,  \
                             ta_x_xxyy_0,   \
                             ta_x_xxyy_1,   \
                             ta_x_xxyyy_0,  \
                             ta_x_xxyyy_1,  \
                             ta_x_xxyyz_0,  \
                             ta_x_xxyyz_1,  \
                             ta_x_xxyz_0,   \
                             ta_x_xxyz_1,   \
                             ta_x_xxyzz_0,  \
                             ta_x_xxyzz_1,  \
                             ta_x_xxzz_0,   \
                             ta_x_xxzz_1,   \
                             ta_x_xxzzz_0,  \
                             ta_x_xxzzz_1,  \
                             ta_x_xyyy_0,   \
                             ta_x_xyyy_1,   \
                             ta_x_xyyyy_0,  \
                             ta_x_xyyyy_1,  \
                             ta_x_xyyyz_0,  \
                             ta_x_xyyyz_1,  \
                             ta_x_xyyz_0,   \
                             ta_x_xyyz_1,   \
                             ta_x_xyyzz_0,  \
                             ta_x_xyyzz_1,  \
                             ta_x_xyzz_0,   \
                             ta_x_xyzz_1,   \
                             ta_x_xyzzz_0,  \
                             ta_x_xyzzz_1,  \
                             ta_x_xzzz_0,   \
                             ta_x_xzzz_1,   \
                             ta_x_xzzzz_0,  \
                             ta_x_xzzzz_1,  \
                             ta_x_yyyy_0,   \
                             ta_x_yyyy_1,   \
                             ta_x_yyyyy_0,  \
                             ta_x_yyyyy_1,  \
                             ta_x_yyyyz_0,  \
                             ta_x_yyyyz_1,  \
                             ta_x_yyyz_0,   \
                             ta_x_yyyz_1,   \
                             ta_x_yyyzz_0,  \
                             ta_x_yyyzz_1,  \
                             ta_x_yyzz_0,   \
                             ta_x_yyzz_1,   \
                             ta_x_yyzzz_0,  \
                             ta_x_yyzzz_1,  \
                             ta_x_yzzz_0,   \
                             ta_x_yzzz_1,   \
                             ta_x_yzzzz_0,  \
                             ta_x_yzzzz_1,  \
                             ta_x_zzzz_0,   \
                             ta_x_zzzz_1,   \
                             ta_x_zzzzz_0,  \
                             ta_x_zzzzz_1,  \
                             ta_xx_xxxxx_0, \
                             ta_xx_xxxxy_0, \
                             ta_xx_xxxxz_0, \
                             ta_xx_xxxyy_0, \
                             ta_xx_xxxyz_0, \
                             ta_xx_xxxzz_0, \
                             ta_xx_xxyyy_0, \
                             ta_xx_xxyyz_0, \
                             ta_xx_xxyzz_0, \
                             ta_xx_xxzzz_0, \
                             ta_xx_xyyyy_0, \
                             ta_xx_xyyyz_0, \
                             ta_xx_xyyzz_0, \
                             ta_xx_xyzzz_0, \
                             ta_xx_xzzzz_0, \
                             ta_xx_yyyyy_0, \
                             ta_xx_yyyyz_0, \
                             ta_xx_yyyzz_0, \
                             ta_xx_yyzzz_0, \
                             ta_xx_yzzzz_0, \
                             ta_xx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xx_xxxxx_0[i] = ta_0_xxxxx_0[i] * fe_0 - ta_0_xxxxx_1[i] * fe_0 + 5.0 * ta_x_xxxx_0[i] * fe_0 - 5.0 * ta_x_xxxx_1[i] * fe_0 +
                           ta_x_xxxxx_0[i] * pa_x[i] - ta_x_xxxxx_1[i] * pc_x[i];

        ta_xx_xxxxy_0[i] = ta_0_xxxxy_0[i] * fe_0 - ta_0_xxxxy_1[i] * fe_0 + 4.0 * ta_x_xxxy_0[i] * fe_0 - 4.0 * ta_x_xxxy_1[i] * fe_0 +
                           ta_x_xxxxy_0[i] * pa_x[i] - ta_x_xxxxy_1[i] * pc_x[i];

        ta_xx_xxxxz_0[i] = ta_0_xxxxz_0[i] * fe_0 - ta_0_xxxxz_1[i] * fe_0 + 4.0 * ta_x_xxxz_0[i] * fe_0 - 4.0 * ta_x_xxxz_1[i] * fe_0 +
                           ta_x_xxxxz_0[i] * pa_x[i] - ta_x_xxxxz_1[i] * pc_x[i];

        ta_xx_xxxyy_0[i] = ta_0_xxxyy_0[i] * fe_0 - ta_0_xxxyy_1[i] * fe_0 + 3.0 * ta_x_xxyy_0[i] * fe_0 - 3.0 * ta_x_xxyy_1[i] * fe_0 +
                           ta_x_xxxyy_0[i] * pa_x[i] - ta_x_xxxyy_1[i] * pc_x[i];

        ta_xx_xxxyz_0[i] = ta_0_xxxyz_0[i] * fe_0 - ta_0_xxxyz_1[i] * fe_0 + 3.0 * ta_x_xxyz_0[i] * fe_0 - 3.0 * ta_x_xxyz_1[i] * fe_0 +
                           ta_x_xxxyz_0[i] * pa_x[i] - ta_x_xxxyz_1[i] * pc_x[i];

        ta_xx_xxxzz_0[i] = ta_0_xxxzz_0[i] * fe_0 - ta_0_xxxzz_1[i] * fe_0 + 3.0 * ta_x_xxzz_0[i] * fe_0 - 3.0 * ta_x_xxzz_1[i] * fe_0 +
                           ta_x_xxxzz_0[i] * pa_x[i] - ta_x_xxxzz_1[i] * pc_x[i];

        ta_xx_xxyyy_0[i] = ta_0_xxyyy_0[i] * fe_0 - ta_0_xxyyy_1[i] * fe_0 + 2.0 * ta_x_xyyy_0[i] * fe_0 - 2.0 * ta_x_xyyy_1[i] * fe_0 +
                           ta_x_xxyyy_0[i] * pa_x[i] - ta_x_xxyyy_1[i] * pc_x[i];

        ta_xx_xxyyz_0[i] = ta_0_xxyyz_0[i] * fe_0 - ta_0_xxyyz_1[i] * fe_0 + 2.0 * ta_x_xyyz_0[i] * fe_0 - 2.0 * ta_x_xyyz_1[i] * fe_0 +
                           ta_x_xxyyz_0[i] * pa_x[i] - ta_x_xxyyz_1[i] * pc_x[i];

        ta_xx_xxyzz_0[i] = ta_0_xxyzz_0[i] * fe_0 - ta_0_xxyzz_1[i] * fe_0 + 2.0 * ta_x_xyzz_0[i] * fe_0 - 2.0 * ta_x_xyzz_1[i] * fe_0 +
                           ta_x_xxyzz_0[i] * pa_x[i] - ta_x_xxyzz_1[i] * pc_x[i];

        ta_xx_xxzzz_0[i] = ta_0_xxzzz_0[i] * fe_0 - ta_0_xxzzz_1[i] * fe_0 + 2.0 * ta_x_xzzz_0[i] * fe_0 - 2.0 * ta_x_xzzz_1[i] * fe_0 +
                           ta_x_xxzzz_0[i] * pa_x[i] - ta_x_xxzzz_1[i] * pc_x[i];

        ta_xx_xyyyy_0[i] = ta_0_xyyyy_0[i] * fe_0 - ta_0_xyyyy_1[i] * fe_0 + ta_x_yyyy_0[i] * fe_0 - ta_x_yyyy_1[i] * fe_0 +
                           ta_x_xyyyy_0[i] * pa_x[i] - ta_x_xyyyy_1[i] * pc_x[i];

        ta_xx_xyyyz_0[i] = ta_0_xyyyz_0[i] * fe_0 - ta_0_xyyyz_1[i] * fe_0 + ta_x_yyyz_0[i] * fe_0 - ta_x_yyyz_1[i] * fe_0 +
                           ta_x_xyyyz_0[i] * pa_x[i] - ta_x_xyyyz_1[i] * pc_x[i];

        ta_xx_xyyzz_0[i] = ta_0_xyyzz_0[i] * fe_0 - ta_0_xyyzz_1[i] * fe_0 + ta_x_yyzz_0[i] * fe_0 - ta_x_yyzz_1[i] * fe_0 +
                           ta_x_xyyzz_0[i] * pa_x[i] - ta_x_xyyzz_1[i] * pc_x[i];

        ta_xx_xyzzz_0[i] = ta_0_xyzzz_0[i] * fe_0 - ta_0_xyzzz_1[i] * fe_0 + ta_x_yzzz_0[i] * fe_0 - ta_x_yzzz_1[i] * fe_0 +
                           ta_x_xyzzz_0[i] * pa_x[i] - ta_x_xyzzz_1[i] * pc_x[i];

        ta_xx_xzzzz_0[i] = ta_0_xzzzz_0[i] * fe_0 - ta_0_xzzzz_1[i] * fe_0 + ta_x_zzzz_0[i] * fe_0 - ta_x_zzzz_1[i] * fe_0 +
                           ta_x_xzzzz_0[i] * pa_x[i] - ta_x_xzzzz_1[i] * pc_x[i];

        ta_xx_yyyyy_0[i] = ta_0_yyyyy_0[i] * fe_0 - ta_0_yyyyy_1[i] * fe_0 + ta_x_yyyyy_0[i] * pa_x[i] - ta_x_yyyyy_1[i] * pc_x[i];

        ta_xx_yyyyz_0[i] = ta_0_yyyyz_0[i] * fe_0 - ta_0_yyyyz_1[i] * fe_0 + ta_x_yyyyz_0[i] * pa_x[i] - ta_x_yyyyz_1[i] * pc_x[i];

        ta_xx_yyyzz_0[i] = ta_0_yyyzz_0[i] * fe_0 - ta_0_yyyzz_1[i] * fe_0 + ta_x_yyyzz_0[i] * pa_x[i] - ta_x_yyyzz_1[i] * pc_x[i];

        ta_xx_yyzzz_0[i] = ta_0_yyzzz_0[i] * fe_0 - ta_0_yyzzz_1[i] * fe_0 + ta_x_yyzzz_0[i] * pa_x[i] - ta_x_yyzzz_1[i] * pc_x[i];

        ta_xx_yzzzz_0[i] = ta_0_yzzzz_0[i] * fe_0 - ta_0_yzzzz_1[i] * fe_0 + ta_x_yzzzz_0[i] * pa_x[i] - ta_x_yzzzz_1[i] * pc_x[i];

        ta_xx_zzzzz_0[i] = ta_0_zzzzz_0[i] * fe_0 - ta_0_zzzzz_1[i] * fe_0 + ta_x_zzzzz_0[i] * pa_x[i] - ta_x_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : DH

    auto ta_xy_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 21);

    auto ta_xy_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 22);

    auto ta_xy_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 23);

    auto ta_xy_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 24);

    auto ta_xy_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 25);

    auto ta_xy_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 26);

    auto ta_xy_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 27);

    auto ta_xy_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 28);

    auto ta_xy_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 29);

    auto ta_xy_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 30);

    auto ta_xy_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 31);

    auto ta_xy_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 32);

    auto ta_xy_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 33);

    auto ta_xy_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 34);

    auto ta_xy_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 35);

    auto ta_xy_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 36);

    auto ta_xy_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 37);

    auto ta_xy_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 38);

    auto ta_xy_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 39);

    auto ta_xy_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 40);

    auto ta_xy_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 41);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_x_xxxxx_0,  \
                             ta_x_xxxxx_1,  \
                             ta_x_xxxxz_0,  \
                             ta_x_xxxxz_1,  \
                             ta_x_xxxzz_0,  \
                             ta_x_xxxzz_1,  \
                             ta_x_xxzzz_0,  \
                             ta_x_xxzzz_1,  \
                             ta_x_xzzzz_0,  \
                             ta_x_xzzzz_1,  \
                             ta_xy_xxxxx_0, \
                             ta_xy_xxxxy_0, \
                             ta_xy_xxxxz_0, \
                             ta_xy_xxxyy_0, \
                             ta_xy_xxxyz_0, \
                             ta_xy_xxxzz_0, \
                             ta_xy_xxyyy_0, \
                             ta_xy_xxyyz_0, \
                             ta_xy_xxyzz_0, \
                             ta_xy_xxzzz_0, \
                             ta_xy_xyyyy_0, \
                             ta_xy_xyyyz_0, \
                             ta_xy_xyyzz_0, \
                             ta_xy_xyzzz_0, \
                             ta_xy_xzzzz_0, \
                             ta_xy_yyyyy_0, \
                             ta_xy_yyyyz_0, \
                             ta_xy_yyyzz_0, \
                             ta_xy_yyzzz_0, \
                             ta_xy_yzzzz_0, \
                             ta_xy_zzzzz_0, \
                             ta_y_xxxxy_0,  \
                             ta_y_xxxxy_1,  \
                             ta_y_xxxy_0,   \
                             ta_y_xxxy_1,   \
                             ta_y_xxxyy_0,  \
                             ta_y_xxxyy_1,  \
                             ta_y_xxxyz_0,  \
                             ta_y_xxxyz_1,  \
                             ta_y_xxyy_0,   \
                             ta_y_xxyy_1,   \
                             ta_y_xxyyy_0,  \
                             ta_y_xxyyy_1,  \
                             ta_y_xxyyz_0,  \
                             ta_y_xxyyz_1,  \
                             ta_y_xxyz_0,   \
                             ta_y_xxyz_1,   \
                             ta_y_xxyzz_0,  \
                             ta_y_xxyzz_1,  \
                             ta_y_xyyy_0,   \
                             ta_y_xyyy_1,   \
                             ta_y_xyyyy_0,  \
                             ta_y_xyyyy_1,  \
                             ta_y_xyyyz_0,  \
                             ta_y_xyyyz_1,  \
                             ta_y_xyyz_0,   \
                             ta_y_xyyz_1,   \
                             ta_y_xyyzz_0,  \
                             ta_y_xyyzz_1,  \
                             ta_y_xyzz_0,   \
                             ta_y_xyzz_1,   \
                             ta_y_xyzzz_0,  \
                             ta_y_xyzzz_1,  \
                             ta_y_yyyy_0,   \
                             ta_y_yyyy_1,   \
                             ta_y_yyyyy_0,  \
                             ta_y_yyyyy_1,  \
                             ta_y_yyyyz_0,  \
                             ta_y_yyyyz_1,  \
                             ta_y_yyyz_0,   \
                             ta_y_yyyz_1,   \
                             ta_y_yyyzz_0,  \
                             ta_y_yyyzz_1,  \
                             ta_y_yyzz_0,   \
                             ta_y_yyzz_1,   \
                             ta_y_yyzzz_0,  \
                             ta_y_yyzzz_1,  \
                             ta_y_yzzz_0,   \
                             ta_y_yzzz_1,   \
                             ta_y_yzzzz_0,  \
                             ta_y_yzzzz_1,  \
                             ta_y_zzzzz_0,  \
                             ta_y_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xy_xxxxx_0[i] = ta_x_xxxxx_0[i] * pa_y[i] - ta_x_xxxxx_1[i] * pc_y[i];

        ta_xy_xxxxy_0[i] = 4.0 * ta_y_xxxy_0[i] * fe_0 - 4.0 * ta_y_xxxy_1[i] * fe_0 + ta_y_xxxxy_0[i] * pa_x[i] - ta_y_xxxxy_1[i] * pc_x[i];

        ta_xy_xxxxz_0[i] = ta_x_xxxxz_0[i] * pa_y[i] - ta_x_xxxxz_1[i] * pc_y[i];

        ta_xy_xxxyy_0[i] = 3.0 * ta_y_xxyy_0[i] * fe_0 - 3.0 * ta_y_xxyy_1[i] * fe_0 + ta_y_xxxyy_0[i] * pa_x[i] - ta_y_xxxyy_1[i] * pc_x[i];

        ta_xy_xxxyz_0[i] = 3.0 * ta_y_xxyz_0[i] * fe_0 - 3.0 * ta_y_xxyz_1[i] * fe_0 + ta_y_xxxyz_0[i] * pa_x[i] - ta_y_xxxyz_1[i] * pc_x[i];

        ta_xy_xxxzz_0[i] = ta_x_xxxzz_0[i] * pa_y[i] - ta_x_xxxzz_1[i] * pc_y[i];

        ta_xy_xxyyy_0[i] = 2.0 * ta_y_xyyy_0[i] * fe_0 - 2.0 * ta_y_xyyy_1[i] * fe_0 + ta_y_xxyyy_0[i] * pa_x[i] - ta_y_xxyyy_1[i] * pc_x[i];

        ta_xy_xxyyz_0[i] = 2.0 * ta_y_xyyz_0[i] * fe_0 - 2.0 * ta_y_xyyz_1[i] * fe_0 + ta_y_xxyyz_0[i] * pa_x[i] - ta_y_xxyyz_1[i] * pc_x[i];

        ta_xy_xxyzz_0[i] = 2.0 * ta_y_xyzz_0[i] * fe_0 - 2.0 * ta_y_xyzz_1[i] * fe_0 + ta_y_xxyzz_0[i] * pa_x[i] - ta_y_xxyzz_1[i] * pc_x[i];

        ta_xy_xxzzz_0[i] = ta_x_xxzzz_0[i] * pa_y[i] - ta_x_xxzzz_1[i] * pc_y[i];

        ta_xy_xyyyy_0[i] = ta_y_yyyy_0[i] * fe_0 - ta_y_yyyy_1[i] * fe_0 + ta_y_xyyyy_0[i] * pa_x[i] - ta_y_xyyyy_1[i] * pc_x[i];

        ta_xy_xyyyz_0[i] = ta_y_yyyz_0[i] * fe_0 - ta_y_yyyz_1[i] * fe_0 + ta_y_xyyyz_0[i] * pa_x[i] - ta_y_xyyyz_1[i] * pc_x[i];

        ta_xy_xyyzz_0[i] = ta_y_yyzz_0[i] * fe_0 - ta_y_yyzz_1[i] * fe_0 + ta_y_xyyzz_0[i] * pa_x[i] - ta_y_xyyzz_1[i] * pc_x[i];

        ta_xy_xyzzz_0[i] = ta_y_yzzz_0[i] * fe_0 - ta_y_yzzz_1[i] * fe_0 + ta_y_xyzzz_0[i] * pa_x[i] - ta_y_xyzzz_1[i] * pc_x[i];

        ta_xy_xzzzz_0[i] = ta_x_xzzzz_0[i] * pa_y[i] - ta_x_xzzzz_1[i] * pc_y[i];

        ta_xy_yyyyy_0[i] = ta_y_yyyyy_0[i] * pa_x[i] - ta_y_yyyyy_1[i] * pc_x[i];

        ta_xy_yyyyz_0[i] = ta_y_yyyyz_0[i] * pa_x[i] - ta_y_yyyyz_1[i] * pc_x[i];

        ta_xy_yyyzz_0[i] = ta_y_yyyzz_0[i] * pa_x[i] - ta_y_yyyzz_1[i] * pc_x[i];

        ta_xy_yyzzz_0[i] = ta_y_yyzzz_0[i] * pa_x[i] - ta_y_yyzzz_1[i] * pc_x[i];

        ta_xy_yzzzz_0[i] = ta_y_yzzzz_0[i] * pa_x[i] - ta_y_yzzzz_1[i] * pc_x[i];

        ta_xy_zzzzz_0[i] = ta_y_zzzzz_0[i] * pa_x[i] - ta_y_zzzzz_1[i] * pc_x[i];
    }

    // Set up 42-63 components of targeted buffer : DH

    auto ta_xz_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 42);

    auto ta_xz_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 43);

    auto ta_xz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 44);

    auto ta_xz_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 45);

    auto ta_xz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 46);

    auto ta_xz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 47);

    auto ta_xz_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 48);

    auto ta_xz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 49);

    auto ta_xz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 50);

    auto ta_xz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 51);

    auto ta_xz_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 52);

    auto ta_xz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 53);

    auto ta_xz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 54);

    auto ta_xz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 55);

    auto ta_xz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 56);

    auto ta_xz_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 57);

    auto ta_xz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 58);

    auto ta_xz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 59);

    auto ta_xz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 60);

    auto ta_xz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 61);

    auto ta_xz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 62);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_x_xxxxx_0,  \
                             ta_x_xxxxx_1,  \
                             ta_x_xxxxy_0,  \
                             ta_x_xxxxy_1,  \
                             ta_x_xxxyy_0,  \
                             ta_x_xxxyy_1,  \
                             ta_x_xxyyy_0,  \
                             ta_x_xxyyy_1,  \
                             ta_x_xyyyy_0,  \
                             ta_x_xyyyy_1,  \
                             ta_xz_xxxxx_0, \
                             ta_xz_xxxxy_0, \
                             ta_xz_xxxxz_0, \
                             ta_xz_xxxyy_0, \
                             ta_xz_xxxyz_0, \
                             ta_xz_xxxzz_0, \
                             ta_xz_xxyyy_0, \
                             ta_xz_xxyyz_0, \
                             ta_xz_xxyzz_0, \
                             ta_xz_xxzzz_0, \
                             ta_xz_xyyyy_0, \
                             ta_xz_xyyyz_0, \
                             ta_xz_xyyzz_0, \
                             ta_xz_xyzzz_0, \
                             ta_xz_xzzzz_0, \
                             ta_xz_yyyyy_0, \
                             ta_xz_yyyyz_0, \
                             ta_xz_yyyzz_0, \
                             ta_xz_yyzzz_0, \
                             ta_xz_yzzzz_0, \
                             ta_xz_zzzzz_0, \
                             ta_z_xxxxz_0,  \
                             ta_z_xxxxz_1,  \
                             ta_z_xxxyz_0,  \
                             ta_z_xxxyz_1,  \
                             ta_z_xxxz_0,   \
                             ta_z_xxxz_1,   \
                             ta_z_xxxzz_0,  \
                             ta_z_xxxzz_1,  \
                             ta_z_xxyyz_0,  \
                             ta_z_xxyyz_1,  \
                             ta_z_xxyz_0,   \
                             ta_z_xxyz_1,   \
                             ta_z_xxyzz_0,  \
                             ta_z_xxyzz_1,  \
                             ta_z_xxzz_0,   \
                             ta_z_xxzz_1,   \
                             ta_z_xxzzz_0,  \
                             ta_z_xxzzz_1,  \
                             ta_z_xyyyz_0,  \
                             ta_z_xyyyz_1,  \
                             ta_z_xyyz_0,   \
                             ta_z_xyyz_1,   \
                             ta_z_xyyzz_0,  \
                             ta_z_xyyzz_1,  \
                             ta_z_xyzz_0,   \
                             ta_z_xyzz_1,   \
                             ta_z_xyzzz_0,  \
                             ta_z_xyzzz_1,  \
                             ta_z_xzzz_0,   \
                             ta_z_xzzz_1,   \
                             ta_z_xzzzz_0,  \
                             ta_z_xzzzz_1,  \
                             ta_z_yyyyy_0,  \
                             ta_z_yyyyy_1,  \
                             ta_z_yyyyz_0,  \
                             ta_z_yyyyz_1,  \
                             ta_z_yyyz_0,   \
                             ta_z_yyyz_1,   \
                             ta_z_yyyzz_0,  \
                             ta_z_yyyzz_1,  \
                             ta_z_yyzz_0,   \
                             ta_z_yyzz_1,   \
                             ta_z_yyzzz_0,  \
                             ta_z_yyzzz_1,  \
                             ta_z_yzzz_0,   \
                             ta_z_yzzz_1,   \
                             ta_z_yzzzz_0,  \
                             ta_z_yzzzz_1,  \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             ta_z_zzzzz_0,  \
                             ta_z_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xz_xxxxx_0[i] = ta_x_xxxxx_0[i] * pa_z[i] - ta_x_xxxxx_1[i] * pc_z[i];

        ta_xz_xxxxy_0[i] = ta_x_xxxxy_0[i] * pa_z[i] - ta_x_xxxxy_1[i] * pc_z[i];

        ta_xz_xxxxz_0[i] = 4.0 * ta_z_xxxz_0[i] * fe_0 - 4.0 * ta_z_xxxz_1[i] * fe_0 + ta_z_xxxxz_0[i] * pa_x[i] - ta_z_xxxxz_1[i] * pc_x[i];

        ta_xz_xxxyy_0[i] = ta_x_xxxyy_0[i] * pa_z[i] - ta_x_xxxyy_1[i] * pc_z[i];

        ta_xz_xxxyz_0[i] = 3.0 * ta_z_xxyz_0[i] * fe_0 - 3.0 * ta_z_xxyz_1[i] * fe_0 + ta_z_xxxyz_0[i] * pa_x[i] - ta_z_xxxyz_1[i] * pc_x[i];

        ta_xz_xxxzz_0[i] = 3.0 * ta_z_xxzz_0[i] * fe_0 - 3.0 * ta_z_xxzz_1[i] * fe_0 + ta_z_xxxzz_0[i] * pa_x[i] - ta_z_xxxzz_1[i] * pc_x[i];

        ta_xz_xxyyy_0[i] = ta_x_xxyyy_0[i] * pa_z[i] - ta_x_xxyyy_1[i] * pc_z[i];

        ta_xz_xxyyz_0[i] = 2.0 * ta_z_xyyz_0[i] * fe_0 - 2.0 * ta_z_xyyz_1[i] * fe_0 + ta_z_xxyyz_0[i] * pa_x[i] - ta_z_xxyyz_1[i] * pc_x[i];

        ta_xz_xxyzz_0[i] = 2.0 * ta_z_xyzz_0[i] * fe_0 - 2.0 * ta_z_xyzz_1[i] * fe_0 + ta_z_xxyzz_0[i] * pa_x[i] - ta_z_xxyzz_1[i] * pc_x[i];

        ta_xz_xxzzz_0[i] = 2.0 * ta_z_xzzz_0[i] * fe_0 - 2.0 * ta_z_xzzz_1[i] * fe_0 + ta_z_xxzzz_0[i] * pa_x[i] - ta_z_xxzzz_1[i] * pc_x[i];

        ta_xz_xyyyy_0[i] = ta_x_xyyyy_0[i] * pa_z[i] - ta_x_xyyyy_1[i] * pc_z[i];

        ta_xz_xyyyz_0[i] = ta_z_yyyz_0[i] * fe_0 - ta_z_yyyz_1[i] * fe_0 + ta_z_xyyyz_0[i] * pa_x[i] - ta_z_xyyyz_1[i] * pc_x[i];

        ta_xz_xyyzz_0[i] = ta_z_yyzz_0[i] * fe_0 - ta_z_yyzz_1[i] * fe_0 + ta_z_xyyzz_0[i] * pa_x[i] - ta_z_xyyzz_1[i] * pc_x[i];

        ta_xz_xyzzz_0[i] = ta_z_yzzz_0[i] * fe_0 - ta_z_yzzz_1[i] * fe_0 + ta_z_xyzzz_0[i] * pa_x[i] - ta_z_xyzzz_1[i] * pc_x[i];

        ta_xz_xzzzz_0[i] = ta_z_zzzz_0[i] * fe_0 - ta_z_zzzz_1[i] * fe_0 + ta_z_xzzzz_0[i] * pa_x[i] - ta_z_xzzzz_1[i] * pc_x[i];

        ta_xz_yyyyy_0[i] = ta_z_yyyyy_0[i] * pa_x[i] - ta_z_yyyyy_1[i] * pc_x[i];

        ta_xz_yyyyz_0[i] = ta_z_yyyyz_0[i] * pa_x[i] - ta_z_yyyyz_1[i] * pc_x[i];

        ta_xz_yyyzz_0[i] = ta_z_yyyzz_0[i] * pa_x[i] - ta_z_yyyzz_1[i] * pc_x[i];

        ta_xz_yyzzz_0[i] = ta_z_yyzzz_0[i] * pa_x[i] - ta_z_yyzzz_1[i] * pc_x[i];

        ta_xz_yzzzz_0[i] = ta_z_yzzzz_0[i] * pa_x[i] - ta_z_yzzzz_1[i] * pc_x[i];

        ta_xz_zzzzz_0[i] = ta_z_zzzzz_0[i] * pa_x[i] - ta_z_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : DH

    auto ta_yy_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 63);

    auto ta_yy_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 64);

    auto ta_yy_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 65);

    auto ta_yy_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 66);

    auto ta_yy_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 67);

    auto ta_yy_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 68);

    auto ta_yy_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 69);

    auto ta_yy_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 70);

    auto ta_yy_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 71);

    auto ta_yy_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 72);

    auto ta_yy_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 73);

    auto ta_yy_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 74);

    auto ta_yy_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 75);

    auto ta_yy_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 76);

    auto ta_yy_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 77);

    auto ta_yy_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 78);

    auto ta_yy_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 79);

    auto ta_yy_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 80);

    auto ta_yy_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 81);

    auto ta_yy_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 82);

    auto ta_yy_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 83);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_0_xxxxx_0,  \
                             ta_0_xxxxx_1,  \
                             ta_0_xxxxy_0,  \
                             ta_0_xxxxy_1,  \
                             ta_0_xxxxz_0,  \
                             ta_0_xxxxz_1,  \
                             ta_0_xxxyy_0,  \
                             ta_0_xxxyy_1,  \
                             ta_0_xxxyz_0,  \
                             ta_0_xxxyz_1,  \
                             ta_0_xxxzz_0,  \
                             ta_0_xxxzz_1,  \
                             ta_0_xxyyy_0,  \
                             ta_0_xxyyy_1,  \
                             ta_0_xxyyz_0,  \
                             ta_0_xxyyz_1,  \
                             ta_0_xxyzz_0,  \
                             ta_0_xxyzz_1,  \
                             ta_0_xxzzz_0,  \
                             ta_0_xxzzz_1,  \
                             ta_0_xyyyy_0,  \
                             ta_0_xyyyy_1,  \
                             ta_0_xyyyz_0,  \
                             ta_0_xyyyz_1,  \
                             ta_0_xyyzz_0,  \
                             ta_0_xyyzz_1,  \
                             ta_0_xyzzz_0,  \
                             ta_0_xyzzz_1,  \
                             ta_0_xzzzz_0,  \
                             ta_0_xzzzz_1,  \
                             ta_0_yyyyy_0,  \
                             ta_0_yyyyy_1,  \
                             ta_0_yyyyz_0,  \
                             ta_0_yyyyz_1,  \
                             ta_0_yyyzz_0,  \
                             ta_0_yyyzz_1,  \
                             ta_0_yyzzz_0,  \
                             ta_0_yyzzz_1,  \
                             ta_0_yzzzz_0,  \
                             ta_0_yzzzz_1,  \
                             ta_0_zzzzz_0,  \
                             ta_0_zzzzz_1,  \
                             ta_y_xxxx_0,   \
                             ta_y_xxxx_1,   \
                             ta_y_xxxxx_0,  \
                             ta_y_xxxxx_1,  \
                             ta_y_xxxxy_0,  \
                             ta_y_xxxxy_1,  \
                             ta_y_xxxxz_0,  \
                             ta_y_xxxxz_1,  \
                             ta_y_xxxy_0,   \
                             ta_y_xxxy_1,   \
                             ta_y_xxxyy_0,  \
                             ta_y_xxxyy_1,  \
                             ta_y_xxxyz_0,  \
                             ta_y_xxxyz_1,  \
                             ta_y_xxxz_0,   \
                             ta_y_xxxz_1,   \
                             ta_y_xxxzz_0,  \
                             ta_y_xxxzz_1,  \
                             ta_y_xxyy_0,   \
                             ta_y_xxyy_1,   \
                             ta_y_xxyyy_0,  \
                             ta_y_xxyyy_1,  \
                             ta_y_xxyyz_0,  \
                             ta_y_xxyyz_1,  \
                             ta_y_xxyz_0,   \
                             ta_y_xxyz_1,   \
                             ta_y_xxyzz_0,  \
                             ta_y_xxyzz_1,  \
                             ta_y_xxzz_0,   \
                             ta_y_xxzz_1,   \
                             ta_y_xxzzz_0,  \
                             ta_y_xxzzz_1,  \
                             ta_y_xyyy_0,   \
                             ta_y_xyyy_1,   \
                             ta_y_xyyyy_0,  \
                             ta_y_xyyyy_1,  \
                             ta_y_xyyyz_0,  \
                             ta_y_xyyyz_1,  \
                             ta_y_xyyz_0,   \
                             ta_y_xyyz_1,   \
                             ta_y_xyyzz_0,  \
                             ta_y_xyyzz_1,  \
                             ta_y_xyzz_0,   \
                             ta_y_xyzz_1,   \
                             ta_y_xyzzz_0,  \
                             ta_y_xyzzz_1,  \
                             ta_y_xzzz_0,   \
                             ta_y_xzzz_1,   \
                             ta_y_xzzzz_0,  \
                             ta_y_xzzzz_1,  \
                             ta_y_yyyy_0,   \
                             ta_y_yyyy_1,   \
                             ta_y_yyyyy_0,  \
                             ta_y_yyyyy_1,  \
                             ta_y_yyyyz_0,  \
                             ta_y_yyyyz_1,  \
                             ta_y_yyyz_0,   \
                             ta_y_yyyz_1,   \
                             ta_y_yyyzz_0,  \
                             ta_y_yyyzz_1,  \
                             ta_y_yyzz_0,   \
                             ta_y_yyzz_1,   \
                             ta_y_yyzzz_0,  \
                             ta_y_yyzzz_1,  \
                             ta_y_yzzz_0,   \
                             ta_y_yzzz_1,   \
                             ta_y_yzzzz_0,  \
                             ta_y_yzzzz_1,  \
                             ta_y_zzzz_0,   \
                             ta_y_zzzz_1,   \
                             ta_y_zzzzz_0,  \
                             ta_y_zzzzz_1,  \
                             ta_yy_xxxxx_0, \
                             ta_yy_xxxxy_0, \
                             ta_yy_xxxxz_0, \
                             ta_yy_xxxyy_0, \
                             ta_yy_xxxyz_0, \
                             ta_yy_xxxzz_0, \
                             ta_yy_xxyyy_0, \
                             ta_yy_xxyyz_0, \
                             ta_yy_xxyzz_0, \
                             ta_yy_xxzzz_0, \
                             ta_yy_xyyyy_0, \
                             ta_yy_xyyyz_0, \
                             ta_yy_xyyzz_0, \
                             ta_yy_xyzzz_0, \
                             ta_yy_xzzzz_0, \
                             ta_yy_yyyyy_0, \
                             ta_yy_yyyyz_0, \
                             ta_yy_yyyzz_0, \
                             ta_yy_yyzzz_0, \
                             ta_yy_yzzzz_0, \
                             ta_yy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yy_xxxxx_0[i] = ta_0_xxxxx_0[i] * fe_0 - ta_0_xxxxx_1[i] * fe_0 + ta_y_xxxxx_0[i] * pa_y[i] - ta_y_xxxxx_1[i] * pc_y[i];

        ta_yy_xxxxy_0[i] = ta_0_xxxxy_0[i] * fe_0 - ta_0_xxxxy_1[i] * fe_0 + ta_y_xxxx_0[i] * fe_0 - ta_y_xxxx_1[i] * fe_0 +
                           ta_y_xxxxy_0[i] * pa_y[i] - ta_y_xxxxy_1[i] * pc_y[i];

        ta_yy_xxxxz_0[i] = ta_0_xxxxz_0[i] * fe_0 - ta_0_xxxxz_1[i] * fe_0 + ta_y_xxxxz_0[i] * pa_y[i] - ta_y_xxxxz_1[i] * pc_y[i];

        ta_yy_xxxyy_0[i] = ta_0_xxxyy_0[i] * fe_0 - ta_0_xxxyy_1[i] * fe_0 + 2.0 * ta_y_xxxy_0[i] * fe_0 - 2.0 * ta_y_xxxy_1[i] * fe_0 +
                           ta_y_xxxyy_0[i] * pa_y[i] - ta_y_xxxyy_1[i] * pc_y[i];

        ta_yy_xxxyz_0[i] = ta_0_xxxyz_0[i] * fe_0 - ta_0_xxxyz_1[i] * fe_0 + ta_y_xxxz_0[i] * fe_0 - ta_y_xxxz_1[i] * fe_0 +
                           ta_y_xxxyz_0[i] * pa_y[i] - ta_y_xxxyz_1[i] * pc_y[i];

        ta_yy_xxxzz_0[i] = ta_0_xxxzz_0[i] * fe_0 - ta_0_xxxzz_1[i] * fe_0 + ta_y_xxxzz_0[i] * pa_y[i] - ta_y_xxxzz_1[i] * pc_y[i];

        ta_yy_xxyyy_0[i] = ta_0_xxyyy_0[i] * fe_0 - ta_0_xxyyy_1[i] * fe_0 + 3.0 * ta_y_xxyy_0[i] * fe_0 - 3.0 * ta_y_xxyy_1[i] * fe_0 +
                           ta_y_xxyyy_0[i] * pa_y[i] - ta_y_xxyyy_1[i] * pc_y[i];

        ta_yy_xxyyz_0[i] = ta_0_xxyyz_0[i] * fe_0 - ta_0_xxyyz_1[i] * fe_0 + 2.0 * ta_y_xxyz_0[i] * fe_0 - 2.0 * ta_y_xxyz_1[i] * fe_0 +
                           ta_y_xxyyz_0[i] * pa_y[i] - ta_y_xxyyz_1[i] * pc_y[i];

        ta_yy_xxyzz_0[i] = ta_0_xxyzz_0[i] * fe_0 - ta_0_xxyzz_1[i] * fe_0 + ta_y_xxzz_0[i] * fe_0 - ta_y_xxzz_1[i] * fe_0 +
                           ta_y_xxyzz_0[i] * pa_y[i] - ta_y_xxyzz_1[i] * pc_y[i];

        ta_yy_xxzzz_0[i] = ta_0_xxzzz_0[i] * fe_0 - ta_0_xxzzz_1[i] * fe_0 + ta_y_xxzzz_0[i] * pa_y[i] - ta_y_xxzzz_1[i] * pc_y[i];

        ta_yy_xyyyy_0[i] = ta_0_xyyyy_0[i] * fe_0 - ta_0_xyyyy_1[i] * fe_0 + 4.0 * ta_y_xyyy_0[i] * fe_0 - 4.0 * ta_y_xyyy_1[i] * fe_0 +
                           ta_y_xyyyy_0[i] * pa_y[i] - ta_y_xyyyy_1[i] * pc_y[i];

        ta_yy_xyyyz_0[i] = ta_0_xyyyz_0[i] * fe_0 - ta_0_xyyyz_1[i] * fe_0 + 3.0 * ta_y_xyyz_0[i] * fe_0 - 3.0 * ta_y_xyyz_1[i] * fe_0 +
                           ta_y_xyyyz_0[i] * pa_y[i] - ta_y_xyyyz_1[i] * pc_y[i];

        ta_yy_xyyzz_0[i] = ta_0_xyyzz_0[i] * fe_0 - ta_0_xyyzz_1[i] * fe_0 + 2.0 * ta_y_xyzz_0[i] * fe_0 - 2.0 * ta_y_xyzz_1[i] * fe_0 +
                           ta_y_xyyzz_0[i] * pa_y[i] - ta_y_xyyzz_1[i] * pc_y[i];

        ta_yy_xyzzz_0[i] = ta_0_xyzzz_0[i] * fe_0 - ta_0_xyzzz_1[i] * fe_0 + ta_y_xzzz_0[i] * fe_0 - ta_y_xzzz_1[i] * fe_0 +
                           ta_y_xyzzz_0[i] * pa_y[i] - ta_y_xyzzz_1[i] * pc_y[i];

        ta_yy_xzzzz_0[i] = ta_0_xzzzz_0[i] * fe_0 - ta_0_xzzzz_1[i] * fe_0 + ta_y_xzzzz_0[i] * pa_y[i] - ta_y_xzzzz_1[i] * pc_y[i];

        ta_yy_yyyyy_0[i] = ta_0_yyyyy_0[i] * fe_0 - ta_0_yyyyy_1[i] * fe_0 + 5.0 * ta_y_yyyy_0[i] * fe_0 - 5.0 * ta_y_yyyy_1[i] * fe_0 +
                           ta_y_yyyyy_0[i] * pa_y[i] - ta_y_yyyyy_1[i] * pc_y[i];

        ta_yy_yyyyz_0[i] = ta_0_yyyyz_0[i] * fe_0 - ta_0_yyyyz_1[i] * fe_0 + 4.0 * ta_y_yyyz_0[i] * fe_0 - 4.0 * ta_y_yyyz_1[i] * fe_0 +
                           ta_y_yyyyz_0[i] * pa_y[i] - ta_y_yyyyz_1[i] * pc_y[i];

        ta_yy_yyyzz_0[i] = ta_0_yyyzz_0[i] * fe_0 - ta_0_yyyzz_1[i] * fe_0 + 3.0 * ta_y_yyzz_0[i] * fe_0 - 3.0 * ta_y_yyzz_1[i] * fe_0 +
                           ta_y_yyyzz_0[i] * pa_y[i] - ta_y_yyyzz_1[i] * pc_y[i];

        ta_yy_yyzzz_0[i] = ta_0_yyzzz_0[i] * fe_0 - ta_0_yyzzz_1[i] * fe_0 + 2.0 * ta_y_yzzz_0[i] * fe_0 - 2.0 * ta_y_yzzz_1[i] * fe_0 +
                           ta_y_yyzzz_0[i] * pa_y[i] - ta_y_yyzzz_1[i] * pc_y[i];

        ta_yy_yzzzz_0[i] = ta_0_yzzzz_0[i] * fe_0 - ta_0_yzzzz_1[i] * fe_0 + ta_y_zzzz_0[i] * fe_0 - ta_y_zzzz_1[i] * fe_0 +
                           ta_y_yzzzz_0[i] * pa_y[i] - ta_y_yzzzz_1[i] * pc_y[i];

        ta_yy_zzzzz_0[i] = ta_0_zzzzz_0[i] * fe_0 - ta_0_zzzzz_1[i] * fe_0 + ta_y_zzzzz_0[i] * pa_y[i] - ta_y_zzzzz_1[i] * pc_y[i];
    }

    // Set up 84-105 components of targeted buffer : DH

    auto ta_yz_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 84);

    auto ta_yz_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 85);

    auto ta_yz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 86);

    auto ta_yz_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 87);

    auto ta_yz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 88);

    auto ta_yz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 89);

    auto ta_yz_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 90);

    auto ta_yz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 91);

    auto ta_yz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 92);

    auto ta_yz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 93);

    auto ta_yz_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 94);

    auto ta_yz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 95);

    auto ta_yz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 96);

    auto ta_yz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 97);

    auto ta_yz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 98);

    auto ta_yz_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 99);

    auto ta_yz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 100);

    auto ta_yz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 101);

    auto ta_yz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 102);

    auto ta_yz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 103);

    auto ta_yz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 104);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_y_xxxxy_0,  \
                             ta_y_xxxxy_1,  \
                             ta_y_xxxyy_0,  \
                             ta_y_xxxyy_1,  \
                             ta_y_xxyyy_0,  \
                             ta_y_xxyyy_1,  \
                             ta_y_xyyyy_0,  \
                             ta_y_xyyyy_1,  \
                             ta_y_yyyyy_0,  \
                             ta_y_yyyyy_1,  \
                             ta_yz_xxxxx_0, \
                             ta_yz_xxxxy_0, \
                             ta_yz_xxxxz_0, \
                             ta_yz_xxxyy_0, \
                             ta_yz_xxxyz_0, \
                             ta_yz_xxxzz_0, \
                             ta_yz_xxyyy_0, \
                             ta_yz_xxyyz_0, \
                             ta_yz_xxyzz_0, \
                             ta_yz_xxzzz_0, \
                             ta_yz_xyyyy_0, \
                             ta_yz_xyyyz_0, \
                             ta_yz_xyyzz_0, \
                             ta_yz_xyzzz_0, \
                             ta_yz_xzzzz_0, \
                             ta_yz_yyyyy_0, \
                             ta_yz_yyyyz_0, \
                             ta_yz_yyyzz_0, \
                             ta_yz_yyzzz_0, \
                             ta_yz_yzzzz_0, \
                             ta_yz_zzzzz_0, \
                             ta_z_xxxxx_0,  \
                             ta_z_xxxxx_1,  \
                             ta_z_xxxxz_0,  \
                             ta_z_xxxxz_1,  \
                             ta_z_xxxyz_0,  \
                             ta_z_xxxyz_1,  \
                             ta_z_xxxz_0,   \
                             ta_z_xxxz_1,   \
                             ta_z_xxxzz_0,  \
                             ta_z_xxxzz_1,  \
                             ta_z_xxyyz_0,  \
                             ta_z_xxyyz_1,  \
                             ta_z_xxyz_0,   \
                             ta_z_xxyz_1,   \
                             ta_z_xxyzz_0,  \
                             ta_z_xxyzz_1,  \
                             ta_z_xxzz_0,   \
                             ta_z_xxzz_1,   \
                             ta_z_xxzzz_0,  \
                             ta_z_xxzzz_1,  \
                             ta_z_xyyyz_0,  \
                             ta_z_xyyyz_1,  \
                             ta_z_xyyz_0,   \
                             ta_z_xyyz_1,   \
                             ta_z_xyyzz_0,  \
                             ta_z_xyyzz_1,  \
                             ta_z_xyzz_0,   \
                             ta_z_xyzz_1,   \
                             ta_z_xyzzz_0,  \
                             ta_z_xyzzz_1,  \
                             ta_z_xzzz_0,   \
                             ta_z_xzzz_1,   \
                             ta_z_xzzzz_0,  \
                             ta_z_xzzzz_1,  \
                             ta_z_yyyyz_0,  \
                             ta_z_yyyyz_1,  \
                             ta_z_yyyz_0,   \
                             ta_z_yyyz_1,   \
                             ta_z_yyyzz_0,  \
                             ta_z_yyyzz_1,  \
                             ta_z_yyzz_0,   \
                             ta_z_yyzz_1,   \
                             ta_z_yyzzz_0,  \
                             ta_z_yyzzz_1,  \
                             ta_z_yzzz_0,   \
                             ta_z_yzzz_1,   \
                             ta_z_yzzzz_0,  \
                             ta_z_yzzzz_1,  \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             ta_z_zzzzz_0,  \
                             ta_z_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yz_xxxxx_0[i] = ta_z_xxxxx_0[i] * pa_y[i] - ta_z_xxxxx_1[i] * pc_y[i];

        ta_yz_xxxxy_0[i] = ta_y_xxxxy_0[i] * pa_z[i] - ta_y_xxxxy_1[i] * pc_z[i];

        ta_yz_xxxxz_0[i] = ta_z_xxxxz_0[i] * pa_y[i] - ta_z_xxxxz_1[i] * pc_y[i];

        ta_yz_xxxyy_0[i] = ta_y_xxxyy_0[i] * pa_z[i] - ta_y_xxxyy_1[i] * pc_z[i];

        ta_yz_xxxyz_0[i] = ta_z_xxxz_0[i] * fe_0 - ta_z_xxxz_1[i] * fe_0 + ta_z_xxxyz_0[i] * pa_y[i] - ta_z_xxxyz_1[i] * pc_y[i];

        ta_yz_xxxzz_0[i] = ta_z_xxxzz_0[i] * pa_y[i] - ta_z_xxxzz_1[i] * pc_y[i];

        ta_yz_xxyyy_0[i] = ta_y_xxyyy_0[i] * pa_z[i] - ta_y_xxyyy_1[i] * pc_z[i];

        ta_yz_xxyyz_0[i] = 2.0 * ta_z_xxyz_0[i] * fe_0 - 2.0 * ta_z_xxyz_1[i] * fe_0 + ta_z_xxyyz_0[i] * pa_y[i] - ta_z_xxyyz_1[i] * pc_y[i];

        ta_yz_xxyzz_0[i] = ta_z_xxzz_0[i] * fe_0 - ta_z_xxzz_1[i] * fe_0 + ta_z_xxyzz_0[i] * pa_y[i] - ta_z_xxyzz_1[i] * pc_y[i];

        ta_yz_xxzzz_0[i] = ta_z_xxzzz_0[i] * pa_y[i] - ta_z_xxzzz_1[i] * pc_y[i];

        ta_yz_xyyyy_0[i] = ta_y_xyyyy_0[i] * pa_z[i] - ta_y_xyyyy_1[i] * pc_z[i];

        ta_yz_xyyyz_0[i] = 3.0 * ta_z_xyyz_0[i] * fe_0 - 3.0 * ta_z_xyyz_1[i] * fe_0 + ta_z_xyyyz_0[i] * pa_y[i] - ta_z_xyyyz_1[i] * pc_y[i];

        ta_yz_xyyzz_0[i] = 2.0 * ta_z_xyzz_0[i] * fe_0 - 2.0 * ta_z_xyzz_1[i] * fe_0 + ta_z_xyyzz_0[i] * pa_y[i] - ta_z_xyyzz_1[i] * pc_y[i];

        ta_yz_xyzzz_0[i] = ta_z_xzzz_0[i] * fe_0 - ta_z_xzzz_1[i] * fe_0 + ta_z_xyzzz_0[i] * pa_y[i] - ta_z_xyzzz_1[i] * pc_y[i];

        ta_yz_xzzzz_0[i] = ta_z_xzzzz_0[i] * pa_y[i] - ta_z_xzzzz_1[i] * pc_y[i];

        ta_yz_yyyyy_0[i] = ta_y_yyyyy_0[i] * pa_z[i] - ta_y_yyyyy_1[i] * pc_z[i];

        ta_yz_yyyyz_0[i] = 4.0 * ta_z_yyyz_0[i] * fe_0 - 4.0 * ta_z_yyyz_1[i] * fe_0 + ta_z_yyyyz_0[i] * pa_y[i] - ta_z_yyyyz_1[i] * pc_y[i];

        ta_yz_yyyzz_0[i] = 3.0 * ta_z_yyzz_0[i] * fe_0 - 3.0 * ta_z_yyzz_1[i] * fe_0 + ta_z_yyyzz_0[i] * pa_y[i] - ta_z_yyyzz_1[i] * pc_y[i];

        ta_yz_yyzzz_0[i] = 2.0 * ta_z_yzzz_0[i] * fe_0 - 2.0 * ta_z_yzzz_1[i] * fe_0 + ta_z_yyzzz_0[i] * pa_y[i] - ta_z_yyzzz_1[i] * pc_y[i];

        ta_yz_yzzzz_0[i] = ta_z_zzzz_0[i] * fe_0 - ta_z_zzzz_1[i] * fe_0 + ta_z_yzzzz_0[i] * pa_y[i] - ta_z_yzzzz_1[i] * pc_y[i];

        ta_yz_zzzzz_0[i] = ta_z_zzzzz_0[i] * pa_y[i] - ta_z_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : DH

    auto ta_zz_xxxxx_0 = pbuffer.data(idx_npot_0_dh + 105);

    auto ta_zz_xxxxy_0 = pbuffer.data(idx_npot_0_dh + 106);

    auto ta_zz_xxxxz_0 = pbuffer.data(idx_npot_0_dh + 107);

    auto ta_zz_xxxyy_0 = pbuffer.data(idx_npot_0_dh + 108);

    auto ta_zz_xxxyz_0 = pbuffer.data(idx_npot_0_dh + 109);

    auto ta_zz_xxxzz_0 = pbuffer.data(idx_npot_0_dh + 110);

    auto ta_zz_xxyyy_0 = pbuffer.data(idx_npot_0_dh + 111);

    auto ta_zz_xxyyz_0 = pbuffer.data(idx_npot_0_dh + 112);

    auto ta_zz_xxyzz_0 = pbuffer.data(idx_npot_0_dh + 113);

    auto ta_zz_xxzzz_0 = pbuffer.data(idx_npot_0_dh + 114);

    auto ta_zz_xyyyy_0 = pbuffer.data(idx_npot_0_dh + 115);

    auto ta_zz_xyyyz_0 = pbuffer.data(idx_npot_0_dh + 116);

    auto ta_zz_xyyzz_0 = pbuffer.data(idx_npot_0_dh + 117);

    auto ta_zz_xyzzz_0 = pbuffer.data(idx_npot_0_dh + 118);

    auto ta_zz_xzzzz_0 = pbuffer.data(idx_npot_0_dh + 119);

    auto ta_zz_yyyyy_0 = pbuffer.data(idx_npot_0_dh + 120);

    auto ta_zz_yyyyz_0 = pbuffer.data(idx_npot_0_dh + 121);

    auto ta_zz_yyyzz_0 = pbuffer.data(idx_npot_0_dh + 122);

    auto ta_zz_yyzzz_0 = pbuffer.data(idx_npot_0_dh + 123);

    auto ta_zz_yzzzz_0 = pbuffer.data(idx_npot_0_dh + 124);

    auto ta_zz_zzzzz_0 = pbuffer.data(idx_npot_0_dh + 125);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta_0_xxxxx_0,  \
                             ta_0_xxxxx_1,  \
                             ta_0_xxxxy_0,  \
                             ta_0_xxxxy_1,  \
                             ta_0_xxxxz_0,  \
                             ta_0_xxxxz_1,  \
                             ta_0_xxxyy_0,  \
                             ta_0_xxxyy_1,  \
                             ta_0_xxxyz_0,  \
                             ta_0_xxxyz_1,  \
                             ta_0_xxxzz_0,  \
                             ta_0_xxxzz_1,  \
                             ta_0_xxyyy_0,  \
                             ta_0_xxyyy_1,  \
                             ta_0_xxyyz_0,  \
                             ta_0_xxyyz_1,  \
                             ta_0_xxyzz_0,  \
                             ta_0_xxyzz_1,  \
                             ta_0_xxzzz_0,  \
                             ta_0_xxzzz_1,  \
                             ta_0_xyyyy_0,  \
                             ta_0_xyyyy_1,  \
                             ta_0_xyyyz_0,  \
                             ta_0_xyyyz_1,  \
                             ta_0_xyyzz_0,  \
                             ta_0_xyyzz_1,  \
                             ta_0_xyzzz_0,  \
                             ta_0_xyzzz_1,  \
                             ta_0_xzzzz_0,  \
                             ta_0_xzzzz_1,  \
                             ta_0_yyyyy_0,  \
                             ta_0_yyyyy_1,  \
                             ta_0_yyyyz_0,  \
                             ta_0_yyyyz_1,  \
                             ta_0_yyyzz_0,  \
                             ta_0_yyyzz_1,  \
                             ta_0_yyzzz_0,  \
                             ta_0_yyzzz_1,  \
                             ta_0_yzzzz_0,  \
                             ta_0_yzzzz_1,  \
                             ta_0_zzzzz_0,  \
                             ta_0_zzzzz_1,  \
                             ta_z_xxxx_0,   \
                             ta_z_xxxx_1,   \
                             ta_z_xxxxx_0,  \
                             ta_z_xxxxx_1,  \
                             ta_z_xxxxy_0,  \
                             ta_z_xxxxy_1,  \
                             ta_z_xxxxz_0,  \
                             ta_z_xxxxz_1,  \
                             ta_z_xxxy_0,   \
                             ta_z_xxxy_1,   \
                             ta_z_xxxyy_0,  \
                             ta_z_xxxyy_1,  \
                             ta_z_xxxyz_0,  \
                             ta_z_xxxyz_1,  \
                             ta_z_xxxz_0,   \
                             ta_z_xxxz_1,   \
                             ta_z_xxxzz_0,  \
                             ta_z_xxxzz_1,  \
                             ta_z_xxyy_0,   \
                             ta_z_xxyy_1,   \
                             ta_z_xxyyy_0,  \
                             ta_z_xxyyy_1,  \
                             ta_z_xxyyz_0,  \
                             ta_z_xxyyz_1,  \
                             ta_z_xxyz_0,   \
                             ta_z_xxyz_1,   \
                             ta_z_xxyzz_0,  \
                             ta_z_xxyzz_1,  \
                             ta_z_xxzz_0,   \
                             ta_z_xxzz_1,   \
                             ta_z_xxzzz_0,  \
                             ta_z_xxzzz_1,  \
                             ta_z_xyyy_0,   \
                             ta_z_xyyy_1,   \
                             ta_z_xyyyy_0,  \
                             ta_z_xyyyy_1,  \
                             ta_z_xyyyz_0,  \
                             ta_z_xyyyz_1,  \
                             ta_z_xyyz_0,   \
                             ta_z_xyyz_1,   \
                             ta_z_xyyzz_0,  \
                             ta_z_xyyzz_1,  \
                             ta_z_xyzz_0,   \
                             ta_z_xyzz_1,   \
                             ta_z_xyzzz_0,  \
                             ta_z_xyzzz_1,  \
                             ta_z_xzzz_0,   \
                             ta_z_xzzz_1,   \
                             ta_z_xzzzz_0,  \
                             ta_z_xzzzz_1,  \
                             ta_z_yyyy_0,   \
                             ta_z_yyyy_1,   \
                             ta_z_yyyyy_0,  \
                             ta_z_yyyyy_1,  \
                             ta_z_yyyyz_0,  \
                             ta_z_yyyyz_1,  \
                             ta_z_yyyz_0,   \
                             ta_z_yyyz_1,   \
                             ta_z_yyyzz_0,  \
                             ta_z_yyyzz_1,  \
                             ta_z_yyzz_0,   \
                             ta_z_yyzz_1,   \
                             ta_z_yyzzz_0,  \
                             ta_z_yyzzz_1,  \
                             ta_z_yzzz_0,   \
                             ta_z_yzzz_1,   \
                             ta_z_yzzzz_0,  \
                             ta_z_yzzzz_1,  \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             ta_z_zzzzz_0,  \
                             ta_z_zzzzz_1,  \
                             ta_zz_xxxxx_0, \
                             ta_zz_xxxxy_0, \
                             ta_zz_xxxxz_0, \
                             ta_zz_xxxyy_0, \
                             ta_zz_xxxyz_0, \
                             ta_zz_xxxzz_0, \
                             ta_zz_xxyyy_0, \
                             ta_zz_xxyyz_0, \
                             ta_zz_xxyzz_0, \
                             ta_zz_xxzzz_0, \
                             ta_zz_xyyyy_0, \
                             ta_zz_xyyyz_0, \
                             ta_zz_xyyzz_0, \
                             ta_zz_xyzzz_0, \
                             ta_zz_xzzzz_0, \
                             ta_zz_yyyyy_0, \
                             ta_zz_yyyyz_0, \
                             ta_zz_yyyzz_0, \
                             ta_zz_yyzzz_0, \
                             ta_zz_yzzzz_0, \
                             ta_zz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zz_xxxxx_0[i] = ta_0_xxxxx_0[i] * fe_0 - ta_0_xxxxx_1[i] * fe_0 + ta_z_xxxxx_0[i] * pa_z[i] - ta_z_xxxxx_1[i] * pc_z[i];

        ta_zz_xxxxy_0[i] = ta_0_xxxxy_0[i] * fe_0 - ta_0_xxxxy_1[i] * fe_0 + ta_z_xxxxy_0[i] * pa_z[i] - ta_z_xxxxy_1[i] * pc_z[i];

        ta_zz_xxxxz_0[i] = ta_0_xxxxz_0[i] * fe_0 - ta_0_xxxxz_1[i] * fe_0 + ta_z_xxxx_0[i] * fe_0 - ta_z_xxxx_1[i] * fe_0 +
                           ta_z_xxxxz_0[i] * pa_z[i] - ta_z_xxxxz_1[i] * pc_z[i];

        ta_zz_xxxyy_0[i] = ta_0_xxxyy_0[i] * fe_0 - ta_0_xxxyy_1[i] * fe_0 + ta_z_xxxyy_0[i] * pa_z[i] - ta_z_xxxyy_1[i] * pc_z[i];

        ta_zz_xxxyz_0[i] = ta_0_xxxyz_0[i] * fe_0 - ta_0_xxxyz_1[i] * fe_0 + ta_z_xxxy_0[i] * fe_0 - ta_z_xxxy_1[i] * fe_0 +
                           ta_z_xxxyz_0[i] * pa_z[i] - ta_z_xxxyz_1[i] * pc_z[i];

        ta_zz_xxxzz_0[i] = ta_0_xxxzz_0[i] * fe_0 - ta_0_xxxzz_1[i] * fe_0 + 2.0 * ta_z_xxxz_0[i] * fe_0 - 2.0 * ta_z_xxxz_1[i] * fe_0 +
                           ta_z_xxxzz_0[i] * pa_z[i] - ta_z_xxxzz_1[i] * pc_z[i];

        ta_zz_xxyyy_0[i] = ta_0_xxyyy_0[i] * fe_0 - ta_0_xxyyy_1[i] * fe_0 + ta_z_xxyyy_0[i] * pa_z[i] - ta_z_xxyyy_1[i] * pc_z[i];

        ta_zz_xxyyz_0[i] = ta_0_xxyyz_0[i] * fe_0 - ta_0_xxyyz_1[i] * fe_0 + ta_z_xxyy_0[i] * fe_0 - ta_z_xxyy_1[i] * fe_0 +
                           ta_z_xxyyz_0[i] * pa_z[i] - ta_z_xxyyz_1[i] * pc_z[i];

        ta_zz_xxyzz_0[i] = ta_0_xxyzz_0[i] * fe_0 - ta_0_xxyzz_1[i] * fe_0 + 2.0 * ta_z_xxyz_0[i] * fe_0 - 2.0 * ta_z_xxyz_1[i] * fe_0 +
                           ta_z_xxyzz_0[i] * pa_z[i] - ta_z_xxyzz_1[i] * pc_z[i];

        ta_zz_xxzzz_0[i] = ta_0_xxzzz_0[i] * fe_0 - ta_0_xxzzz_1[i] * fe_0 + 3.0 * ta_z_xxzz_0[i] * fe_0 - 3.0 * ta_z_xxzz_1[i] * fe_0 +
                           ta_z_xxzzz_0[i] * pa_z[i] - ta_z_xxzzz_1[i] * pc_z[i];

        ta_zz_xyyyy_0[i] = ta_0_xyyyy_0[i] * fe_0 - ta_0_xyyyy_1[i] * fe_0 + ta_z_xyyyy_0[i] * pa_z[i] - ta_z_xyyyy_1[i] * pc_z[i];

        ta_zz_xyyyz_0[i] = ta_0_xyyyz_0[i] * fe_0 - ta_0_xyyyz_1[i] * fe_0 + ta_z_xyyy_0[i] * fe_0 - ta_z_xyyy_1[i] * fe_0 +
                           ta_z_xyyyz_0[i] * pa_z[i] - ta_z_xyyyz_1[i] * pc_z[i];

        ta_zz_xyyzz_0[i] = ta_0_xyyzz_0[i] * fe_0 - ta_0_xyyzz_1[i] * fe_0 + 2.0 * ta_z_xyyz_0[i] * fe_0 - 2.0 * ta_z_xyyz_1[i] * fe_0 +
                           ta_z_xyyzz_0[i] * pa_z[i] - ta_z_xyyzz_1[i] * pc_z[i];

        ta_zz_xyzzz_0[i] = ta_0_xyzzz_0[i] * fe_0 - ta_0_xyzzz_1[i] * fe_0 + 3.0 * ta_z_xyzz_0[i] * fe_0 - 3.0 * ta_z_xyzz_1[i] * fe_0 +
                           ta_z_xyzzz_0[i] * pa_z[i] - ta_z_xyzzz_1[i] * pc_z[i];

        ta_zz_xzzzz_0[i] = ta_0_xzzzz_0[i] * fe_0 - ta_0_xzzzz_1[i] * fe_0 + 4.0 * ta_z_xzzz_0[i] * fe_0 - 4.0 * ta_z_xzzz_1[i] * fe_0 +
                           ta_z_xzzzz_0[i] * pa_z[i] - ta_z_xzzzz_1[i] * pc_z[i];

        ta_zz_yyyyy_0[i] = ta_0_yyyyy_0[i] * fe_0 - ta_0_yyyyy_1[i] * fe_0 + ta_z_yyyyy_0[i] * pa_z[i] - ta_z_yyyyy_1[i] * pc_z[i];

        ta_zz_yyyyz_0[i] = ta_0_yyyyz_0[i] * fe_0 - ta_0_yyyyz_1[i] * fe_0 + ta_z_yyyy_0[i] * fe_0 - ta_z_yyyy_1[i] * fe_0 +
                           ta_z_yyyyz_0[i] * pa_z[i] - ta_z_yyyyz_1[i] * pc_z[i];

        ta_zz_yyyzz_0[i] = ta_0_yyyzz_0[i] * fe_0 - ta_0_yyyzz_1[i] * fe_0 + 2.0 * ta_z_yyyz_0[i] * fe_0 - 2.0 * ta_z_yyyz_1[i] * fe_0 +
                           ta_z_yyyzz_0[i] * pa_z[i] - ta_z_yyyzz_1[i] * pc_z[i];

        ta_zz_yyzzz_0[i] = ta_0_yyzzz_0[i] * fe_0 - ta_0_yyzzz_1[i] * fe_0 + 3.0 * ta_z_yyzz_0[i] * fe_0 - 3.0 * ta_z_yyzz_1[i] * fe_0 +
                           ta_z_yyzzz_0[i] * pa_z[i] - ta_z_yyzzz_1[i] * pc_z[i];

        ta_zz_yzzzz_0[i] = ta_0_yzzzz_0[i] * fe_0 - ta_0_yzzzz_1[i] * fe_0 + 4.0 * ta_z_yzzz_0[i] * fe_0 - 4.0 * ta_z_yzzz_1[i] * fe_0 +
                           ta_z_yzzzz_0[i] * pa_z[i] - ta_z_yzzzz_1[i] * pc_z[i];

        ta_zz_zzzzz_0[i] = ta_0_zzzzz_0[i] * fe_0 - ta_0_zzzzz_1[i] * fe_0 + 5.0 * ta_z_zzzz_0[i] * fe_0 - 5.0 * ta_z_zzzz_1[i] * fe_0 +
                           ta_z_zzzzz_0[i] * pa_z[i] - ta_z_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
