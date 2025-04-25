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

#include "NuclearPotentialGeom010PrimRecDH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_dh(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_dh,
                                        const size_t              idx_npot_geom_010_0_sh,
                                        const size_t              idx_npot_geom_010_1_sh,
                                        const size_t              idx_npot_geom_010_0_pg,
                                        const size_t              idx_npot_geom_010_1_pg,
                                        const size_t              idx_npot_1_ph,
                                        const size_t              idx_npot_geom_010_0_ph,
                                        const size_t              idx_npot_geom_010_1_ph,
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

    auto ta1_x_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh);

    auto ta1_x_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 1);

    auto ta1_x_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 2);

    auto ta1_x_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 3);

    auto ta1_x_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 4);

    auto ta1_x_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 5);

    auto ta1_x_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 6);

    auto ta1_x_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 7);

    auto ta1_x_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 8);

    auto ta1_x_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 9);

    auto ta1_x_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 10);

    auto ta1_x_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 11);

    auto ta1_x_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 12);

    auto ta1_x_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 13);

    auto ta1_x_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 14);

    auto ta1_x_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 15);

    auto ta1_x_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 16);

    auto ta1_x_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 17);

    auto ta1_x_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 18);

    auto ta1_x_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 19);

    auto ta1_x_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 20);

    auto ta1_y_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 21);

    auto ta1_y_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 22);

    auto ta1_y_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 23);

    auto ta1_y_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 24);

    auto ta1_y_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 25);

    auto ta1_y_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 26);

    auto ta1_y_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 27);

    auto ta1_y_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 28);

    auto ta1_y_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 29);

    auto ta1_y_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 30);

    auto ta1_y_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 31);

    auto ta1_y_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 32);

    auto ta1_y_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 33);

    auto ta1_y_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 34);

    auto ta1_y_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 35);

    auto ta1_y_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 36);

    auto ta1_y_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 37);

    auto ta1_y_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 38);

    auto ta1_y_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 39);

    auto ta1_y_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 40);

    auto ta1_y_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 41);

    auto ta1_z_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 42);

    auto ta1_z_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 43);

    auto ta1_z_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 44);

    auto ta1_z_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 45);

    auto ta1_z_0_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 46);

    auto ta1_z_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 47);

    auto ta1_z_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 48);

    auto ta1_z_0_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 49);

    auto ta1_z_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 50);

    auto ta1_z_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 51);

    auto ta1_z_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 52);

    auto ta1_z_0_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 53);

    auto ta1_z_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 54);

    auto ta1_z_0_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 55);

    auto ta1_z_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 56);

    auto ta1_z_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 57);

    auto ta1_z_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 58);

    auto ta1_z_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 59);

    auto ta1_z_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 60);

    auto ta1_z_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 61);

    auto ta1_z_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 62);

    // Set up components of auxiliary buffer : SH

    auto ta1_x_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh);

    auto ta1_x_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 1);

    auto ta1_x_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 2);

    auto ta1_x_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 3);

    auto ta1_x_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 4);

    auto ta1_x_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 5);

    auto ta1_x_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 6);

    auto ta1_x_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 7);

    auto ta1_x_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 8);

    auto ta1_x_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 9);

    auto ta1_x_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 10);

    auto ta1_x_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 11);

    auto ta1_x_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 12);

    auto ta1_x_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 13);

    auto ta1_x_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 14);

    auto ta1_x_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 15);

    auto ta1_x_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 16);

    auto ta1_x_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 17);

    auto ta1_x_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 18);

    auto ta1_x_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 19);

    auto ta1_x_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 20);

    auto ta1_y_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 21);

    auto ta1_y_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 22);

    auto ta1_y_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 23);

    auto ta1_y_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 24);

    auto ta1_y_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 25);

    auto ta1_y_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 26);

    auto ta1_y_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 27);

    auto ta1_y_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 28);

    auto ta1_y_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 29);

    auto ta1_y_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 30);

    auto ta1_y_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 31);

    auto ta1_y_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 32);

    auto ta1_y_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 33);

    auto ta1_y_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 34);

    auto ta1_y_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 35);

    auto ta1_y_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 36);

    auto ta1_y_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 37);

    auto ta1_y_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 38);

    auto ta1_y_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 39);

    auto ta1_y_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 40);

    auto ta1_y_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 41);

    auto ta1_z_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 42);

    auto ta1_z_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 43);

    auto ta1_z_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 44);

    auto ta1_z_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 45);

    auto ta1_z_0_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 46);

    auto ta1_z_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 47);

    auto ta1_z_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 48);

    auto ta1_z_0_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 49);

    auto ta1_z_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 50);

    auto ta1_z_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 51);

    auto ta1_z_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 52);

    auto ta1_z_0_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 53);

    auto ta1_z_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 54);

    auto ta1_z_0_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 55);

    auto ta1_z_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 56);

    auto ta1_z_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 57);

    auto ta1_z_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 58);

    auto ta1_z_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 59);

    auto ta1_z_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 60);

    auto ta1_z_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 61);

    auto ta1_z_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 62);

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg);

    auto ta1_x_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 1);

    auto ta1_x_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 2);

    auto ta1_x_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 3);

    auto ta1_x_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 4);

    auto ta1_x_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 5);

    auto ta1_x_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 6);

    auto ta1_x_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 7);

    auto ta1_x_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 8);

    auto ta1_x_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 9);

    auto ta1_x_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 10);

    auto ta1_x_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 11);

    auto ta1_x_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 12);

    auto ta1_x_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 13);

    auto ta1_x_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 14);

    auto ta1_x_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 15);

    auto ta1_x_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 16);

    auto ta1_x_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 17);

    auto ta1_x_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 18);

    auto ta1_x_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 19);

    auto ta1_x_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 20);

    auto ta1_x_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 21);

    auto ta1_x_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 22);

    auto ta1_x_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 23);

    auto ta1_x_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 24);

    auto ta1_x_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 25);

    auto ta1_x_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 26);

    auto ta1_x_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 27);

    auto ta1_x_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 28);

    auto ta1_x_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 29);

    auto ta1_x_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 30);

    auto ta1_x_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 31);

    auto ta1_x_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 32);

    auto ta1_x_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 33);

    auto ta1_x_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 34);

    auto ta1_x_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 35);

    auto ta1_x_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 36);

    auto ta1_x_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 37);

    auto ta1_x_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 38);

    auto ta1_x_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 39);

    auto ta1_x_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 40);

    auto ta1_x_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 41);

    auto ta1_x_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 42);

    auto ta1_x_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 43);

    auto ta1_x_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 44);

    auto ta1_y_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 45);

    auto ta1_y_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 46);

    auto ta1_y_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 47);

    auto ta1_y_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 48);

    auto ta1_y_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 49);

    auto ta1_y_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 50);

    auto ta1_y_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 51);

    auto ta1_y_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 52);

    auto ta1_y_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 53);

    auto ta1_y_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 54);

    auto ta1_y_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 55);

    auto ta1_y_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 56);

    auto ta1_y_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 57);

    auto ta1_y_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 58);

    auto ta1_y_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 59);

    auto ta1_y_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 60);

    auto ta1_y_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 61);

    auto ta1_y_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 62);

    auto ta1_y_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 63);

    auto ta1_y_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 64);

    auto ta1_y_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 65);

    auto ta1_y_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 66);

    auto ta1_y_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 67);

    auto ta1_y_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 68);

    auto ta1_y_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 69);

    auto ta1_y_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 70);

    auto ta1_y_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 71);

    auto ta1_y_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 72);

    auto ta1_y_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 73);

    auto ta1_y_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 74);

    auto ta1_y_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 75);

    auto ta1_y_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 76);

    auto ta1_y_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 77);

    auto ta1_y_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 78);

    auto ta1_y_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 79);

    auto ta1_y_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 80);

    auto ta1_y_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 81);

    auto ta1_y_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 82);

    auto ta1_y_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 83);

    auto ta1_y_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 84);

    auto ta1_y_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 85);

    auto ta1_y_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 86);

    auto ta1_y_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 87);

    auto ta1_y_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 88);

    auto ta1_y_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 89);

    auto ta1_z_x_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 90);

    auto ta1_z_x_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 91);

    auto ta1_z_x_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 92);

    auto ta1_z_x_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 93);

    auto ta1_z_x_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 94);

    auto ta1_z_x_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 95);

    auto ta1_z_x_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 96);

    auto ta1_z_x_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 97);

    auto ta1_z_x_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 98);

    auto ta1_z_x_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 99);

    auto ta1_z_x_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 100);

    auto ta1_z_x_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 101);

    auto ta1_z_x_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 102);

    auto ta1_z_x_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 103);

    auto ta1_z_x_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 104);

    auto ta1_z_y_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 105);

    auto ta1_z_y_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 106);

    auto ta1_z_y_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 107);

    auto ta1_z_y_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 108);

    auto ta1_z_y_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 109);

    auto ta1_z_y_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 110);

    auto ta1_z_y_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 111);

    auto ta1_z_y_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 112);

    auto ta1_z_y_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 113);

    auto ta1_z_y_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 114);

    auto ta1_z_y_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 115);

    auto ta1_z_y_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 116);

    auto ta1_z_y_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 117);

    auto ta1_z_y_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 118);

    auto ta1_z_y_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 119);

    auto ta1_z_z_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_pg + 120);

    auto ta1_z_z_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 121);

    auto ta1_z_z_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 122);

    auto ta1_z_z_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 123);

    auto ta1_z_z_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 124);

    auto ta1_z_z_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 125);

    auto ta1_z_z_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 126);

    auto ta1_z_z_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 127);

    auto ta1_z_z_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 128);

    auto ta1_z_z_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 129);

    auto ta1_z_z_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_pg + 130);

    auto ta1_z_z_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 131);

    auto ta1_z_z_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 132);

    auto ta1_z_z_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 133);

    auto ta1_z_z_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_pg + 134);

    // Set up components of auxiliary buffer : PG

    auto ta1_x_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg);

    auto ta1_x_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 1);

    auto ta1_x_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 2);

    auto ta1_x_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 3);

    auto ta1_x_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 4);

    auto ta1_x_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 5);

    auto ta1_x_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 6);

    auto ta1_x_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 7);

    auto ta1_x_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 8);

    auto ta1_x_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 9);

    auto ta1_x_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 10);

    auto ta1_x_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 11);

    auto ta1_x_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 12);

    auto ta1_x_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 13);

    auto ta1_x_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 14);

    auto ta1_x_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 15);

    auto ta1_x_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 16);

    auto ta1_x_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 17);

    auto ta1_x_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 18);

    auto ta1_x_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 19);

    auto ta1_x_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 20);

    auto ta1_x_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 21);

    auto ta1_x_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 22);

    auto ta1_x_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 23);

    auto ta1_x_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 24);

    auto ta1_x_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 25);

    auto ta1_x_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 26);

    auto ta1_x_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 27);

    auto ta1_x_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 28);

    auto ta1_x_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 29);

    auto ta1_x_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 30);

    auto ta1_x_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 31);

    auto ta1_x_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 32);

    auto ta1_x_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 33);

    auto ta1_x_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 34);

    auto ta1_x_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 35);

    auto ta1_x_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 36);

    auto ta1_x_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 37);

    auto ta1_x_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 38);

    auto ta1_x_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 39);

    auto ta1_x_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 40);

    auto ta1_x_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 41);

    auto ta1_x_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 42);

    auto ta1_x_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 43);

    auto ta1_x_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 44);

    auto ta1_y_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 45);

    auto ta1_y_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 46);

    auto ta1_y_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 47);

    auto ta1_y_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 48);

    auto ta1_y_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 49);

    auto ta1_y_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 50);

    auto ta1_y_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 51);

    auto ta1_y_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 52);

    auto ta1_y_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 53);

    auto ta1_y_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 54);

    auto ta1_y_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 55);

    auto ta1_y_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 56);

    auto ta1_y_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 57);

    auto ta1_y_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 58);

    auto ta1_y_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 59);

    auto ta1_y_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 60);

    auto ta1_y_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 61);

    auto ta1_y_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 62);

    auto ta1_y_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 63);

    auto ta1_y_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 64);

    auto ta1_y_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 65);

    auto ta1_y_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 66);

    auto ta1_y_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 67);

    auto ta1_y_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 68);

    auto ta1_y_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 69);

    auto ta1_y_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 70);

    auto ta1_y_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 71);

    auto ta1_y_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 72);

    auto ta1_y_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 73);

    auto ta1_y_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 74);

    auto ta1_y_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 75);

    auto ta1_y_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 76);

    auto ta1_y_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 77);

    auto ta1_y_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 78);

    auto ta1_y_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 79);

    auto ta1_y_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 80);

    auto ta1_y_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 81);

    auto ta1_y_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 82);

    auto ta1_y_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 83);

    auto ta1_y_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 84);

    auto ta1_y_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 85);

    auto ta1_y_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 86);

    auto ta1_y_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 87);

    auto ta1_y_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 88);

    auto ta1_y_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 89);

    auto ta1_z_x_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 90);

    auto ta1_z_x_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 91);

    auto ta1_z_x_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 92);

    auto ta1_z_x_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 93);

    auto ta1_z_x_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 94);

    auto ta1_z_x_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 95);

    auto ta1_z_x_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 96);

    auto ta1_z_x_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 97);

    auto ta1_z_x_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 98);

    auto ta1_z_x_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 99);

    auto ta1_z_x_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 100);

    auto ta1_z_x_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 101);

    auto ta1_z_x_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 102);

    auto ta1_z_x_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 103);

    auto ta1_z_x_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 104);

    auto ta1_z_y_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 105);

    auto ta1_z_y_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 106);

    auto ta1_z_y_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 107);

    auto ta1_z_y_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 108);

    auto ta1_z_y_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 109);

    auto ta1_z_y_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 110);

    auto ta1_z_y_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 111);

    auto ta1_z_y_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 112);

    auto ta1_z_y_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 113);

    auto ta1_z_y_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 114);

    auto ta1_z_y_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 115);

    auto ta1_z_y_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 116);

    auto ta1_z_y_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 117);

    auto ta1_z_y_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 118);

    auto ta1_z_y_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 119);

    auto ta1_z_z_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_pg + 120);

    auto ta1_z_z_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 121);

    auto ta1_z_z_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 122);

    auto ta1_z_z_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 123);

    auto ta1_z_z_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 124);

    auto ta1_z_z_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 125);

    auto ta1_z_z_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 126);

    auto ta1_z_z_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 127);

    auto ta1_z_z_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 128);

    auto ta1_z_z_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 129);

    auto ta1_z_z_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_pg + 130);

    auto ta1_z_z_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 131);

    auto ta1_z_z_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 132);

    auto ta1_z_z_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 133);

    auto ta1_z_z_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_pg + 134);

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

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph);

    auto ta1_x_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 1);

    auto ta1_x_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 2);

    auto ta1_x_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 3);

    auto ta1_x_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 4);

    auto ta1_x_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 5);

    auto ta1_x_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 6);

    auto ta1_x_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 7);

    auto ta1_x_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 8);

    auto ta1_x_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 9);

    auto ta1_x_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 10);

    auto ta1_x_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 11);

    auto ta1_x_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 12);

    auto ta1_x_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 13);

    auto ta1_x_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 14);

    auto ta1_x_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 15);

    auto ta1_x_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 16);

    auto ta1_x_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 17);

    auto ta1_x_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 18);

    auto ta1_x_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 19);

    auto ta1_x_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 20);

    auto ta1_x_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 21);

    auto ta1_x_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 22);

    auto ta1_x_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 23);

    auto ta1_x_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 24);

    auto ta1_x_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 25);

    auto ta1_x_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 26);

    auto ta1_x_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 27);

    auto ta1_x_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 28);

    auto ta1_x_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 29);

    auto ta1_x_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 30);

    auto ta1_x_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 31);

    auto ta1_x_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 32);

    auto ta1_x_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 33);

    auto ta1_x_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 34);

    auto ta1_x_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 35);

    auto ta1_x_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 36);

    auto ta1_x_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 37);

    auto ta1_x_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 38);

    auto ta1_x_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 39);

    auto ta1_x_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 40);

    auto ta1_x_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 41);

    auto ta1_x_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 42);

    auto ta1_x_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 43);

    auto ta1_x_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 44);

    auto ta1_x_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 45);

    auto ta1_x_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 46);

    auto ta1_x_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 47);

    auto ta1_x_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 48);

    auto ta1_x_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 49);

    auto ta1_x_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 50);

    auto ta1_x_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 51);

    auto ta1_x_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 52);

    auto ta1_x_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 53);

    auto ta1_x_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 54);

    auto ta1_x_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 55);

    auto ta1_x_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 56);

    auto ta1_x_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 57);

    auto ta1_x_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 58);

    auto ta1_x_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 59);

    auto ta1_x_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 60);

    auto ta1_x_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 61);

    auto ta1_x_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 62);

    auto ta1_y_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 63);

    auto ta1_y_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 64);

    auto ta1_y_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 65);

    auto ta1_y_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 66);

    auto ta1_y_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 67);

    auto ta1_y_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 68);

    auto ta1_y_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 69);

    auto ta1_y_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 70);

    auto ta1_y_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 71);

    auto ta1_y_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 72);

    auto ta1_y_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 73);

    auto ta1_y_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 74);

    auto ta1_y_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 75);

    auto ta1_y_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 76);

    auto ta1_y_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 77);

    auto ta1_y_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 78);

    auto ta1_y_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 79);

    auto ta1_y_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 80);

    auto ta1_y_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 81);

    auto ta1_y_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 82);

    auto ta1_y_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 83);

    auto ta1_y_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 84);

    auto ta1_y_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 85);

    auto ta1_y_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 86);

    auto ta1_y_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 87);

    auto ta1_y_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 88);

    auto ta1_y_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 89);

    auto ta1_y_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 90);

    auto ta1_y_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 91);

    auto ta1_y_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 92);

    auto ta1_y_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 93);

    auto ta1_y_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 94);

    auto ta1_y_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 95);

    auto ta1_y_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 96);

    auto ta1_y_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 97);

    auto ta1_y_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 98);

    auto ta1_y_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 99);

    auto ta1_y_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 100);

    auto ta1_y_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 101);

    auto ta1_y_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 102);

    auto ta1_y_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 103);

    auto ta1_y_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 104);

    auto ta1_y_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 105);

    auto ta1_y_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 106);

    auto ta1_y_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 107);

    auto ta1_y_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 108);

    auto ta1_y_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 109);

    auto ta1_y_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 110);

    auto ta1_y_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 111);

    auto ta1_y_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 112);

    auto ta1_y_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 113);

    auto ta1_y_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 114);

    auto ta1_y_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 115);

    auto ta1_y_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 116);

    auto ta1_y_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 117);

    auto ta1_y_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 118);

    auto ta1_y_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 119);

    auto ta1_y_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 120);

    auto ta1_y_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 121);

    auto ta1_y_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 122);

    auto ta1_y_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 123);

    auto ta1_y_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 124);

    auto ta1_y_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 125);

    auto ta1_z_x_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 126);

    auto ta1_z_x_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 127);

    auto ta1_z_x_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 128);

    auto ta1_z_x_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 129);

    auto ta1_z_x_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 130);

    auto ta1_z_x_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 131);

    auto ta1_z_x_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 132);

    auto ta1_z_x_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 133);

    auto ta1_z_x_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 134);

    auto ta1_z_x_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 135);

    auto ta1_z_x_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 136);

    auto ta1_z_x_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 137);

    auto ta1_z_x_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 138);

    auto ta1_z_x_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 139);

    auto ta1_z_x_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 140);

    auto ta1_z_x_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 141);

    auto ta1_z_x_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 142);

    auto ta1_z_x_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 143);

    auto ta1_z_x_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 144);

    auto ta1_z_x_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 145);

    auto ta1_z_x_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 146);

    auto ta1_z_y_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 147);

    auto ta1_z_y_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 148);

    auto ta1_z_y_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 149);

    auto ta1_z_y_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 150);

    auto ta1_z_y_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 151);

    auto ta1_z_y_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 152);

    auto ta1_z_y_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 153);

    auto ta1_z_y_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 154);

    auto ta1_z_y_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 155);

    auto ta1_z_y_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 156);

    auto ta1_z_y_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 157);

    auto ta1_z_y_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 158);

    auto ta1_z_y_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 159);

    auto ta1_z_y_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 160);

    auto ta1_z_y_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 161);

    auto ta1_z_y_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 162);

    auto ta1_z_y_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 163);

    auto ta1_z_y_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 164);

    auto ta1_z_y_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 165);

    auto ta1_z_y_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 166);

    auto ta1_z_y_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 167);

    auto ta1_z_z_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_ph + 168);

    auto ta1_z_z_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 169);

    auto ta1_z_z_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 170);

    auto ta1_z_z_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 171);

    auto ta1_z_z_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 172);

    auto ta1_z_z_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 173);

    auto ta1_z_z_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 174);

    auto ta1_z_z_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 175);

    auto ta1_z_z_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 176);

    auto ta1_z_z_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 177);

    auto ta1_z_z_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 178);

    auto ta1_z_z_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 179);

    auto ta1_z_z_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 180);

    auto ta1_z_z_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 181);

    auto ta1_z_z_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 182);

    auto ta1_z_z_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_ph + 183);

    auto ta1_z_z_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 184);

    auto ta1_z_z_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 185);

    auto ta1_z_z_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 186);

    auto ta1_z_z_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 187);

    auto ta1_z_z_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_ph + 188);

    // Set up components of auxiliary buffer : PH

    auto ta1_x_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph);

    auto ta1_x_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 1);

    auto ta1_x_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 2);

    auto ta1_x_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 3);

    auto ta1_x_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 4);

    auto ta1_x_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 5);

    auto ta1_x_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 6);

    auto ta1_x_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 7);

    auto ta1_x_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 8);

    auto ta1_x_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 9);

    auto ta1_x_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 10);

    auto ta1_x_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 11);

    auto ta1_x_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 12);

    auto ta1_x_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 13);

    auto ta1_x_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 14);

    auto ta1_x_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 15);

    auto ta1_x_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 16);

    auto ta1_x_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 17);

    auto ta1_x_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 18);

    auto ta1_x_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 19);

    auto ta1_x_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 20);

    auto ta1_x_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 21);

    auto ta1_x_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 22);

    auto ta1_x_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 23);

    auto ta1_x_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 24);

    auto ta1_x_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 25);

    auto ta1_x_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 26);

    auto ta1_x_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 27);

    auto ta1_x_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 28);

    auto ta1_x_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 29);

    auto ta1_x_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 30);

    auto ta1_x_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 31);

    auto ta1_x_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 32);

    auto ta1_x_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 33);

    auto ta1_x_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 34);

    auto ta1_x_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 35);

    auto ta1_x_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 36);

    auto ta1_x_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 37);

    auto ta1_x_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 38);

    auto ta1_x_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 39);

    auto ta1_x_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 40);

    auto ta1_x_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 41);

    auto ta1_x_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 42);

    auto ta1_x_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 43);

    auto ta1_x_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 44);

    auto ta1_x_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 45);

    auto ta1_x_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 46);

    auto ta1_x_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 47);

    auto ta1_x_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 48);

    auto ta1_x_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 49);

    auto ta1_x_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 50);

    auto ta1_x_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 51);

    auto ta1_x_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 52);

    auto ta1_x_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 53);

    auto ta1_x_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 54);

    auto ta1_x_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 55);

    auto ta1_x_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 56);

    auto ta1_x_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 57);

    auto ta1_x_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 58);

    auto ta1_x_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 59);

    auto ta1_x_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 60);

    auto ta1_x_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 61);

    auto ta1_x_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 62);

    auto ta1_y_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 63);

    auto ta1_y_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 64);

    auto ta1_y_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 65);

    auto ta1_y_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 66);

    auto ta1_y_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 67);

    auto ta1_y_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 68);

    auto ta1_y_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 69);

    auto ta1_y_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 70);

    auto ta1_y_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 71);

    auto ta1_y_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 72);

    auto ta1_y_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 73);

    auto ta1_y_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 74);

    auto ta1_y_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 75);

    auto ta1_y_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 76);

    auto ta1_y_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 77);

    auto ta1_y_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 78);

    auto ta1_y_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 79);

    auto ta1_y_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 80);

    auto ta1_y_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 81);

    auto ta1_y_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 82);

    auto ta1_y_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 83);

    auto ta1_y_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 84);

    auto ta1_y_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 85);

    auto ta1_y_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 86);

    auto ta1_y_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 87);

    auto ta1_y_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 88);

    auto ta1_y_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 89);

    auto ta1_y_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 90);

    auto ta1_y_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 91);

    auto ta1_y_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 92);

    auto ta1_y_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 93);

    auto ta1_y_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 94);

    auto ta1_y_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 95);

    auto ta1_y_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 96);

    auto ta1_y_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 97);

    auto ta1_y_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 98);

    auto ta1_y_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 99);

    auto ta1_y_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 100);

    auto ta1_y_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 101);

    auto ta1_y_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 102);

    auto ta1_y_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 103);

    auto ta1_y_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 104);

    auto ta1_y_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 105);

    auto ta1_y_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 106);

    auto ta1_y_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 107);

    auto ta1_y_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 108);

    auto ta1_y_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 109);

    auto ta1_y_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 110);

    auto ta1_y_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 111);

    auto ta1_y_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 112);

    auto ta1_y_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 113);

    auto ta1_y_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 114);

    auto ta1_y_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 115);

    auto ta1_y_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 116);

    auto ta1_y_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 117);

    auto ta1_y_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 118);

    auto ta1_y_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 119);

    auto ta1_y_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 120);

    auto ta1_y_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 121);

    auto ta1_y_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 122);

    auto ta1_y_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 123);

    auto ta1_y_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 124);

    auto ta1_y_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 125);

    auto ta1_z_x_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 126);

    auto ta1_z_x_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 127);

    auto ta1_z_x_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 128);

    auto ta1_z_x_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 129);

    auto ta1_z_x_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 130);

    auto ta1_z_x_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 131);

    auto ta1_z_x_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 132);

    auto ta1_z_x_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 133);

    auto ta1_z_x_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 134);

    auto ta1_z_x_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 135);

    auto ta1_z_x_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 136);

    auto ta1_z_x_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 137);

    auto ta1_z_x_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 138);

    auto ta1_z_x_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 139);

    auto ta1_z_x_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 140);

    auto ta1_z_x_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 141);

    auto ta1_z_x_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 142);

    auto ta1_z_x_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 143);

    auto ta1_z_x_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 144);

    auto ta1_z_x_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 145);

    auto ta1_z_x_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 146);

    auto ta1_z_y_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 147);

    auto ta1_z_y_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 148);

    auto ta1_z_y_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 149);

    auto ta1_z_y_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 150);

    auto ta1_z_y_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 151);

    auto ta1_z_y_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 152);

    auto ta1_z_y_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 153);

    auto ta1_z_y_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 154);

    auto ta1_z_y_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 155);

    auto ta1_z_y_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 156);

    auto ta1_z_y_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 157);

    auto ta1_z_y_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 158);

    auto ta1_z_y_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 159);

    auto ta1_z_y_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 160);

    auto ta1_z_y_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 161);

    auto ta1_z_y_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 162);

    auto ta1_z_y_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 163);

    auto ta1_z_y_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 164);

    auto ta1_z_y_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 165);

    auto ta1_z_y_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 166);

    auto ta1_z_y_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 167);

    auto ta1_z_z_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_ph + 168);

    auto ta1_z_z_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 169);

    auto ta1_z_z_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 170);

    auto ta1_z_z_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 171);

    auto ta1_z_z_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 172);

    auto ta1_z_z_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 173);

    auto ta1_z_z_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 174);

    auto ta1_z_z_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 175);

    auto ta1_z_z_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 176);

    auto ta1_z_z_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 177);

    auto ta1_z_z_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 178);

    auto ta1_z_z_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 179);

    auto ta1_z_z_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 180);

    auto ta1_z_z_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 181);

    auto ta1_z_z_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 182);

    auto ta1_z_z_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_ph + 183);

    auto ta1_z_z_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 184);

    auto ta1_z_z_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 185);

    auto ta1_z_z_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 186);

    auto ta1_z_z_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 187);

    auto ta1_z_z_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_ph + 188);

    // Set up 0-21 components of targeted buffer : DH

    auto ta1_x_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh);

    auto ta1_x_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 1);

    auto ta1_x_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 2);

    auto ta1_x_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 3);

    auto ta1_x_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 4);

    auto ta1_x_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 5);

    auto ta1_x_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 6);

    auto ta1_x_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 7);

    auto ta1_x_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 8);

    auto ta1_x_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 9);

    auto ta1_x_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 10);

    auto ta1_x_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 11);

    auto ta1_x_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 12);

    auto ta1_x_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 13);

    auto ta1_x_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 14);

    auto ta1_x_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 15);

    auto ta1_x_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 16);

    auto ta1_x_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 17);

    auto ta1_x_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 18);

    auto ta1_x_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 19);

    auto ta1_x_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 20);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxxx_0,  \
                             ta1_x_x_xxxxx_1,  \
                             ta1_x_x_xxxxy_0,  \
                             ta1_x_x_xxxxy_1,  \
                             ta1_x_x_xxxxz_0,  \
                             ta1_x_x_xxxxz_1,  \
                             ta1_x_x_xxxy_0,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxxyy_0,  \
                             ta1_x_x_xxxyy_1,  \
                             ta1_x_x_xxxyz_0,  \
                             ta1_x_x_xxxyz_1,  \
                             ta1_x_x_xxxz_0,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxxzz_0,  \
                             ta1_x_x_xxxzz_1,  \
                             ta1_x_x_xxyy_0,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xxyyy_0,  \
                             ta1_x_x_xxyyy_1,  \
                             ta1_x_x_xxyyz_0,  \
                             ta1_x_x_xxyyz_1,  \
                             ta1_x_x_xxyz_0,   \
                             ta1_x_x_xxyz_1,   \
                             ta1_x_x_xxyzz_0,  \
                             ta1_x_x_xxyzz_1,  \
                             ta1_x_x_xxzz_0,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xxzzz_0,  \
                             ta1_x_x_xxzzz_1,  \
                             ta1_x_x_xyyy_0,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_x_xyyyy_0,  \
                             ta1_x_x_xyyyy_1,  \
                             ta1_x_x_xyyyz_0,  \
                             ta1_x_x_xyyyz_1,  \
                             ta1_x_x_xyyz_0,   \
                             ta1_x_x_xyyz_1,   \
                             ta1_x_x_xyyzz_0,  \
                             ta1_x_x_xyyzz_1,  \
                             ta1_x_x_xyzz_0,   \
                             ta1_x_x_xyzz_1,   \
                             ta1_x_x_xyzzz_0,  \
                             ta1_x_x_xyzzz_1,  \
                             ta1_x_x_xzzz_0,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_x_xzzzz_0,  \
                             ta1_x_x_xzzzz_1,  \
                             ta1_x_x_yyyy_0,   \
                             ta1_x_x_yyyy_1,   \
                             ta1_x_x_yyyyy_0,  \
                             ta1_x_x_yyyyy_1,  \
                             ta1_x_x_yyyyz_0,  \
                             ta1_x_x_yyyyz_1,  \
                             ta1_x_x_yyyz_0,   \
                             ta1_x_x_yyyz_1,   \
                             ta1_x_x_yyyzz_0,  \
                             ta1_x_x_yyyzz_1,  \
                             ta1_x_x_yyzz_0,   \
                             ta1_x_x_yyzz_1,   \
                             ta1_x_x_yyzzz_0,  \
                             ta1_x_x_yyzzz_1,  \
                             ta1_x_x_yzzz_0,   \
                             ta1_x_x_yzzz_1,   \
                             ta1_x_x_yzzzz_0,  \
                             ta1_x_x_yzzzz_1,  \
                             ta1_x_x_zzzz_0,   \
                             ta1_x_x_zzzz_1,   \
                             ta1_x_x_zzzzz_0,  \
                             ta1_x_x_zzzzz_1,  \
                             ta1_x_xx_xxxxx_0, \
                             ta1_x_xx_xxxxy_0, \
                             ta1_x_xx_xxxxz_0, \
                             ta1_x_xx_xxxyy_0, \
                             ta1_x_xx_xxxyz_0, \
                             ta1_x_xx_xxxzz_0, \
                             ta1_x_xx_xxyyy_0, \
                             ta1_x_xx_xxyyz_0, \
                             ta1_x_xx_xxyzz_0, \
                             ta1_x_xx_xxzzz_0, \
                             ta1_x_xx_xyyyy_0, \
                             ta1_x_xx_xyyyz_0, \
                             ta1_x_xx_xyyzz_0, \
                             ta1_x_xx_xyzzz_0, \
                             ta1_x_xx_xzzzz_0, \
                             ta1_x_xx_yyyyy_0, \
                             ta1_x_xx_yyyyz_0, \
                             ta1_x_xx_yyyzz_0, \
                             ta1_x_xx_yyzzz_0, \
                             ta1_x_xx_yzzzz_0, \
                             ta1_x_xx_zzzzz_0, \
                             ta_x_xxxxx_1,     \
                             ta_x_xxxxy_1,     \
                             ta_x_xxxxz_1,     \
                             ta_x_xxxyy_1,     \
                             ta_x_xxxyz_1,     \
                             ta_x_xxxzz_1,     \
                             ta_x_xxyyy_1,     \
                             ta_x_xxyyz_1,     \
                             ta_x_xxyzz_1,     \
                             ta_x_xxzzz_1,     \
                             ta_x_xyyyy_1,     \
                             ta_x_xyyyz_1,     \
                             ta_x_xyyzz_1,     \
                             ta_x_xyzzz_1,     \
                             ta_x_xzzzz_1,     \
                             ta_x_yyyyy_1,     \
                             ta_x_yyyyz_1,     \
                             ta_x_yyyzz_1,     \
                             ta_x_yyzzz_1,     \
                             ta_x_yzzzz_1,     \
                             ta_x_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xx_xxxxx_0[i] = ta1_x_0_xxxxx_0[i] * fe_0 - ta1_x_0_xxxxx_1[i] * fe_0 + 5.0 * ta1_x_x_xxxx_0[i] * fe_0 -
                              5.0 * ta1_x_x_xxxx_1[i] * fe_0 + ta_x_xxxxx_1[i] + ta1_x_x_xxxxx_0[i] * pa_x[i] - ta1_x_x_xxxxx_1[i] * pc_x[i];

        ta1_x_xx_xxxxy_0[i] = ta1_x_0_xxxxy_0[i] * fe_0 - ta1_x_0_xxxxy_1[i] * fe_0 + 4.0 * ta1_x_x_xxxy_0[i] * fe_0 -
                              4.0 * ta1_x_x_xxxy_1[i] * fe_0 + ta_x_xxxxy_1[i] + ta1_x_x_xxxxy_0[i] * pa_x[i] - ta1_x_x_xxxxy_1[i] * pc_x[i];

        ta1_x_xx_xxxxz_0[i] = ta1_x_0_xxxxz_0[i] * fe_0 - ta1_x_0_xxxxz_1[i] * fe_0 + 4.0 * ta1_x_x_xxxz_0[i] * fe_0 -
                              4.0 * ta1_x_x_xxxz_1[i] * fe_0 + ta_x_xxxxz_1[i] + ta1_x_x_xxxxz_0[i] * pa_x[i] - ta1_x_x_xxxxz_1[i] * pc_x[i];

        ta1_x_xx_xxxyy_0[i] = ta1_x_0_xxxyy_0[i] * fe_0 - ta1_x_0_xxxyy_1[i] * fe_0 + 3.0 * ta1_x_x_xxyy_0[i] * fe_0 -
                              3.0 * ta1_x_x_xxyy_1[i] * fe_0 + ta_x_xxxyy_1[i] + ta1_x_x_xxxyy_0[i] * pa_x[i] - ta1_x_x_xxxyy_1[i] * pc_x[i];

        ta1_x_xx_xxxyz_0[i] = ta1_x_0_xxxyz_0[i] * fe_0 - ta1_x_0_xxxyz_1[i] * fe_0 + 3.0 * ta1_x_x_xxyz_0[i] * fe_0 -
                              3.0 * ta1_x_x_xxyz_1[i] * fe_0 + ta_x_xxxyz_1[i] + ta1_x_x_xxxyz_0[i] * pa_x[i] - ta1_x_x_xxxyz_1[i] * pc_x[i];

        ta1_x_xx_xxxzz_0[i] = ta1_x_0_xxxzz_0[i] * fe_0 - ta1_x_0_xxxzz_1[i] * fe_0 + 3.0 * ta1_x_x_xxzz_0[i] * fe_0 -
                              3.0 * ta1_x_x_xxzz_1[i] * fe_0 + ta_x_xxxzz_1[i] + ta1_x_x_xxxzz_0[i] * pa_x[i] - ta1_x_x_xxxzz_1[i] * pc_x[i];

        ta1_x_xx_xxyyy_0[i] = ta1_x_0_xxyyy_0[i] * fe_0 - ta1_x_0_xxyyy_1[i] * fe_0 + 2.0 * ta1_x_x_xyyy_0[i] * fe_0 -
                              2.0 * ta1_x_x_xyyy_1[i] * fe_0 + ta_x_xxyyy_1[i] + ta1_x_x_xxyyy_0[i] * pa_x[i] - ta1_x_x_xxyyy_1[i] * pc_x[i];

        ta1_x_xx_xxyyz_0[i] = ta1_x_0_xxyyz_0[i] * fe_0 - ta1_x_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_x_xyyz_0[i] * fe_0 -
                              2.0 * ta1_x_x_xyyz_1[i] * fe_0 + ta_x_xxyyz_1[i] + ta1_x_x_xxyyz_0[i] * pa_x[i] - ta1_x_x_xxyyz_1[i] * pc_x[i];

        ta1_x_xx_xxyzz_0[i] = ta1_x_0_xxyzz_0[i] * fe_0 - ta1_x_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_x_xyzz_0[i] * fe_0 -
                              2.0 * ta1_x_x_xyzz_1[i] * fe_0 + ta_x_xxyzz_1[i] + ta1_x_x_xxyzz_0[i] * pa_x[i] - ta1_x_x_xxyzz_1[i] * pc_x[i];

        ta1_x_xx_xxzzz_0[i] = ta1_x_0_xxzzz_0[i] * fe_0 - ta1_x_0_xxzzz_1[i] * fe_0 + 2.0 * ta1_x_x_xzzz_0[i] * fe_0 -
                              2.0 * ta1_x_x_xzzz_1[i] * fe_0 + ta_x_xxzzz_1[i] + ta1_x_x_xxzzz_0[i] * pa_x[i] - ta1_x_x_xxzzz_1[i] * pc_x[i];

        ta1_x_xx_xyyyy_0[i] = ta1_x_0_xyyyy_0[i] * fe_0 - ta1_x_0_xyyyy_1[i] * fe_0 + ta1_x_x_yyyy_0[i] * fe_0 - ta1_x_x_yyyy_1[i] * fe_0 +
                              ta_x_xyyyy_1[i] + ta1_x_x_xyyyy_0[i] * pa_x[i] - ta1_x_x_xyyyy_1[i] * pc_x[i];

        ta1_x_xx_xyyyz_0[i] = ta1_x_0_xyyyz_0[i] * fe_0 - ta1_x_0_xyyyz_1[i] * fe_0 + ta1_x_x_yyyz_0[i] * fe_0 - ta1_x_x_yyyz_1[i] * fe_0 +
                              ta_x_xyyyz_1[i] + ta1_x_x_xyyyz_0[i] * pa_x[i] - ta1_x_x_xyyyz_1[i] * pc_x[i];

        ta1_x_xx_xyyzz_0[i] = ta1_x_0_xyyzz_0[i] * fe_0 - ta1_x_0_xyyzz_1[i] * fe_0 + ta1_x_x_yyzz_0[i] * fe_0 - ta1_x_x_yyzz_1[i] * fe_0 +
                              ta_x_xyyzz_1[i] + ta1_x_x_xyyzz_0[i] * pa_x[i] - ta1_x_x_xyyzz_1[i] * pc_x[i];

        ta1_x_xx_xyzzz_0[i] = ta1_x_0_xyzzz_0[i] * fe_0 - ta1_x_0_xyzzz_1[i] * fe_0 + ta1_x_x_yzzz_0[i] * fe_0 - ta1_x_x_yzzz_1[i] * fe_0 +
                              ta_x_xyzzz_1[i] + ta1_x_x_xyzzz_0[i] * pa_x[i] - ta1_x_x_xyzzz_1[i] * pc_x[i];

        ta1_x_xx_xzzzz_0[i] = ta1_x_0_xzzzz_0[i] * fe_0 - ta1_x_0_xzzzz_1[i] * fe_0 + ta1_x_x_zzzz_0[i] * fe_0 - ta1_x_x_zzzz_1[i] * fe_0 +
                              ta_x_xzzzz_1[i] + ta1_x_x_xzzzz_0[i] * pa_x[i] - ta1_x_x_xzzzz_1[i] * pc_x[i];

        ta1_x_xx_yyyyy_0[i] =
            ta1_x_0_yyyyy_0[i] * fe_0 - ta1_x_0_yyyyy_1[i] * fe_0 + ta_x_yyyyy_1[i] + ta1_x_x_yyyyy_0[i] * pa_x[i] - ta1_x_x_yyyyy_1[i] * pc_x[i];

        ta1_x_xx_yyyyz_0[i] =
            ta1_x_0_yyyyz_0[i] * fe_0 - ta1_x_0_yyyyz_1[i] * fe_0 + ta_x_yyyyz_1[i] + ta1_x_x_yyyyz_0[i] * pa_x[i] - ta1_x_x_yyyyz_1[i] * pc_x[i];

        ta1_x_xx_yyyzz_0[i] =
            ta1_x_0_yyyzz_0[i] * fe_0 - ta1_x_0_yyyzz_1[i] * fe_0 + ta_x_yyyzz_1[i] + ta1_x_x_yyyzz_0[i] * pa_x[i] - ta1_x_x_yyyzz_1[i] * pc_x[i];

        ta1_x_xx_yyzzz_0[i] =
            ta1_x_0_yyzzz_0[i] * fe_0 - ta1_x_0_yyzzz_1[i] * fe_0 + ta_x_yyzzz_1[i] + ta1_x_x_yyzzz_0[i] * pa_x[i] - ta1_x_x_yyzzz_1[i] * pc_x[i];

        ta1_x_xx_yzzzz_0[i] =
            ta1_x_0_yzzzz_0[i] * fe_0 - ta1_x_0_yzzzz_1[i] * fe_0 + ta_x_yzzzz_1[i] + ta1_x_x_yzzzz_0[i] * pa_x[i] - ta1_x_x_yzzzz_1[i] * pc_x[i];

        ta1_x_xx_zzzzz_0[i] =
            ta1_x_0_zzzzz_0[i] * fe_0 - ta1_x_0_zzzzz_1[i] * fe_0 + ta_x_zzzzz_1[i] + ta1_x_x_zzzzz_0[i] * pa_x[i] - ta1_x_x_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : DH

    auto ta1_x_xy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 21);

    auto ta1_x_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 22);

    auto ta1_x_xy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 23);

    auto ta1_x_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 24);

    auto ta1_x_xy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 25);

    auto ta1_x_xy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 26);

    auto ta1_x_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 27);

    auto ta1_x_xy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 28);

    auto ta1_x_xy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 29);

    auto ta1_x_xy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 30);

    auto ta1_x_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 31);

    auto ta1_x_xy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 32);

    auto ta1_x_xy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 33);

    auto ta1_x_xy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 34);

    auto ta1_x_xy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 35);

    auto ta1_x_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 36);

    auto ta1_x_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 37);

    auto ta1_x_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 38);

    auto ta1_x_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 39);

    auto ta1_x_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 40);

    auto ta1_x_xy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 41);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxxx_0,  \
                             ta1_x_x_xxxxx_1,  \
                             ta1_x_x_xxxxy_0,  \
                             ta1_x_x_xxxxy_1,  \
                             ta1_x_x_xxxxz_0,  \
                             ta1_x_x_xxxxz_1,  \
                             ta1_x_x_xxxy_0,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxxyy_0,  \
                             ta1_x_x_xxxyy_1,  \
                             ta1_x_x_xxxyz_0,  \
                             ta1_x_x_xxxyz_1,  \
                             ta1_x_x_xxxz_0,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxxzz_0,  \
                             ta1_x_x_xxxzz_1,  \
                             ta1_x_x_xxyy_0,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xxyyy_0,  \
                             ta1_x_x_xxyyy_1,  \
                             ta1_x_x_xxyyz_0,  \
                             ta1_x_x_xxyyz_1,  \
                             ta1_x_x_xxyz_0,   \
                             ta1_x_x_xxyz_1,   \
                             ta1_x_x_xxyzz_0,  \
                             ta1_x_x_xxyzz_1,  \
                             ta1_x_x_xxzz_0,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xxzzz_0,  \
                             ta1_x_x_xxzzz_1,  \
                             ta1_x_x_xyyy_0,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_x_xyyyy_0,  \
                             ta1_x_x_xyyyy_1,  \
                             ta1_x_x_xyyyz_0,  \
                             ta1_x_x_xyyyz_1,  \
                             ta1_x_x_xyyz_0,   \
                             ta1_x_x_xyyz_1,   \
                             ta1_x_x_xyyzz_0,  \
                             ta1_x_x_xyyzz_1,  \
                             ta1_x_x_xyzz_0,   \
                             ta1_x_x_xyzz_1,   \
                             ta1_x_x_xyzzz_0,  \
                             ta1_x_x_xyzzz_1,  \
                             ta1_x_x_xzzz_0,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_x_xzzzz_0,  \
                             ta1_x_x_xzzzz_1,  \
                             ta1_x_x_zzzzz_0,  \
                             ta1_x_x_zzzzz_1,  \
                             ta1_x_xy_xxxxx_0, \
                             ta1_x_xy_xxxxy_0, \
                             ta1_x_xy_xxxxz_0, \
                             ta1_x_xy_xxxyy_0, \
                             ta1_x_xy_xxxyz_0, \
                             ta1_x_xy_xxxzz_0, \
                             ta1_x_xy_xxyyy_0, \
                             ta1_x_xy_xxyyz_0, \
                             ta1_x_xy_xxyzz_0, \
                             ta1_x_xy_xxzzz_0, \
                             ta1_x_xy_xyyyy_0, \
                             ta1_x_xy_xyyyz_0, \
                             ta1_x_xy_xyyzz_0, \
                             ta1_x_xy_xyzzz_0, \
                             ta1_x_xy_xzzzz_0, \
                             ta1_x_xy_yyyyy_0, \
                             ta1_x_xy_yyyyz_0, \
                             ta1_x_xy_yyyzz_0, \
                             ta1_x_xy_yyzzz_0, \
                             ta1_x_xy_yzzzz_0, \
                             ta1_x_xy_zzzzz_0, \
                             ta1_x_y_yyyyy_0,  \
                             ta1_x_y_yyyyy_1,  \
                             ta1_x_y_yyyyz_0,  \
                             ta1_x_y_yyyyz_1,  \
                             ta1_x_y_yyyzz_0,  \
                             ta1_x_y_yyyzz_1,  \
                             ta1_x_y_yyzzz_0,  \
                             ta1_x_y_yyzzz_1,  \
                             ta1_x_y_yzzzz_0,  \
                             ta1_x_y_yzzzz_1,  \
                             ta_y_yyyyy_1,     \
                             ta_y_yyyyz_1,     \
                             ta_y_yyyzz_1,     \
                             ta_y_yyzzz_1,     \
                             ta_y_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xy_xxxxx_0[i] = ta1_x_x_xxxxx_0[i] * pa_y[i] - ta1_x_x_xxxxx_1[i] * pc_y[i];

        ta1_x_xy_xxxxy_0[i] = ta1_x_x_xxxx_0[i] * fe_0 - ta1_x_x_xxxx_1[i] * fe_0 + ta1_x_x_xxxxy_0[i] * pa_y[i] - ta1_x_x_xxxxy_1[i] * pc_y[i];

        ta1_x_xy_xxxxz_0[i] = ta1_x_x_xxxxz_0[i] * pa_y[i] - ta1_x_x_xxxxz_1[i] * pc_y[i];

        ta1_x_xy_xxxyy_0[i] =
            2.0 * ta1_x_x_xxxy_0[i] * fe_0 - 2.0 * ta1_x_x_xxxy_1[i] * fe_0 + ta1_x_x_xxxyy_0[i] * pa_y[i] - ta1_x_x_xxxyy_1[i] * pc_y[i];

        ta1_x_xy_xxxyz_0[i] = ta1_x_x_xxxz_0[i] * fe_0 - ta1_x_x_xxxz_1[i] * fe_0 + ta1_x_x_xxxyz_0[i] * pa_y[i] - ta1_x_x_xxxyz_1[i] * pc_y[i];

        ta1_x_xy_xxxzz_0[i] = ta1_x_x_xxxzz_0[i] * pa_y[i] - ta1_x_x_xxxzz_1[i] * pc_y[i];

        ta1_x_xy_xxyyy_0[i] =
            3.0 * ta1_x_x_xxyy_0[i] * fe_0 - 3.0 * ta1_x_x_xxyy_1[i] * fe_0 + ta1_x_x_xxyyy_0[i] * pa_y[i] - ta1_x_x_xxyyy_1[i] * pc_y[i];

        ta1_x_xy_xxyyz_0[i] =
            2.0 * ta1_x_x_xxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyz_1[i] * fe_0 + ta1_x_x_xxyyz_0[i] * pa_y[i] - ta1_x_x_xxyyz_1[i] * pc_y[i];

        ta1_x_xy_xxyzz_0[i] = ta1_x_x_xxzz_0[i] * fe_0 - ta1_x_x_xxzz_1[i] * fe_0 + ta1_x_x_xxyzz_0[i] * pa_y[i] - ta1_x_x_xxyzz_1[i] * pc_y[i];

        ta1_x_xy_xxzzz_0[i] = ta1_x_x_xxzzz_0[i] * pa_y[i] - ta1_x_x_xxzzz_1[i] * pc_y[i];

        ta1_x_xy_xyyyy_0[i] =
            4.0 * ta1_x_x_xyyy_0[i] * fe_0 - 4.0 * ta1_x_x_xyyy_1[i] * fe_0 + ta1_x_x_xyyyy_0[i] * pa_y[i] - ta1_x_x_xyyyy_1[i] * pc_y[i];

        ta1_x_xy_xyyyz_0[i] =
            3.0 * ta1_x_x_xyyz_0[i] * fe_0 - 3.0 * ta1_x_x_xyyz_1[i] * fe_0 + ta1_x_x_xyyyz_0[i] * pa_y[i] - ta1_x_x_xyyyz_1[i] * pc_y[i];

        ta1_x_xy_xyyzz_0[i] =
            2.0 * ta1_x_x_xyzz_0[i] * fe_0 - 2.0 * ta1_x_x_xyzz_1[i] * fe_0 + ta1_x_x_xyyzz_0[i] * pa_y[i] - ta1_x_x_xyyzz_1[i] * pc_y[i];

        ta1_x_xy_xyzzz_0[i] = ta1_x_x_xzzz_0[i] * fe_0 - ta1_x_x_xzzz_1[i] * fe_0 + ta1_x_x_xyzzz_0[i] * pa_y[i] - ta1_x_x_xyzzz_1[i] * pc_y[i];

        ta1_x_xy_xzzzz_0[i] = ta1_x_x_xzzzz_0[i] * pa_y[i] - ta1_x_x_xzzzz_1[i] * pc_y[i];

        ta1_x_xy_yyyyy_0[i] = ta_y_yyyyy_1[i] + ta1_x_y_yyyyy_0[i] * pa_x[i] - ta1_x_y_yyyyy_1[i] * pc_x[i];

        ta1_x_xy_yyyyz_0[i] = ta_y_yyyyz_1[i] + ta1_x_y_yyyyz_0[i] * pa_x[i] - ta1_x_y_yyyyz_1[i] * pc_x[i];

        ta1_x_xy_yyyzz_0[i] = ta_y_yyyzz_1[i] + ta1_x_y_yyyzz_0[i] * pa_x[i] - ta1_x_y_yyyzz_1[i] * pc_x[i];

        ta1_x_xy_yyzzz_0[i] = ta_y_yyzzz_1[i] + ta1_x_y_yyzzz_0[i] * pa_x[i] - ta1_x_y_yyzzz_1[i] * pc_x[i];

        ta1_x_xy_yzzzz_0[i] = ta_y_yzzzz_1[i] + ta1_x_y_yzzzz_0[i] * pa_x[i] - ta1_x_y_yzzzz_1[i] * pc_x[i];

        ta1_x_xy_zzzzz_0[i] = ta1_x_x_zzzzz_0[i] * pa_y[i] - ta1_x_x_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : DH

    auto ta1_x_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 42);

    auto ta1_x_xz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 43);

    auto ta1_x_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 44);

    auto ta1_x_xz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 45);

    auto ta1_x_xz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 46);

    auto ta1_x_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 47);

    auto ta1_x_xz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 48);

    auto ta1_x_xz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 49);

    auto ta1_x_xz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 50);

    auto ta1_x_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 51);

    auto ta1_x_xz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 52);

    auto ta1_x_xz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 53);

    auto ta1_x_xz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 54);

    auto ta1_x_xz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 55);

    auto ta1_x_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 56);

    auto ta1_x_xz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 57);

    auto ta1_x_xz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 58);

    auto ta1_x_xz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 59);

    auto ta1_x_xz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 60);

    auto ta1_x_xz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 61);

    auto ta1_x_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 62);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_x_xxxx_0,   \
                             ta1_x_x_xxxx_1,   \
                             ta1_x_x_xxxxx_0,  \
                             ta1_x_x_xxxxx_1,  \
                             ta1_x_x_xxxxy_0,  \
                             ta1_x_x_xxxxy_1,  \
                             ta1_x_x_xxxxz_0,  \
                             ta1_x_x_xxxxz_1,  \
                             ta1_x_x_xxxy_0,   \
                             ta1_x_x_xxxy_1,   \
                             ta1_x_x_xxxyy_0,  \
                             ta1_x_x_xxxyy_1,  \
                             ta1_x_x_xxxyz_0,  \
                             ta1_x_x_xxxyz_1,  \
                             ta1_x_x_xxxz_0,   \
                             ta1_x_x_xxxz_1,   \
                             ta1_x_x_xxxzz_0,  \
                             ta1_x_x_xxxzz_1,  \
                             ta1_x_x_xxyy_0,   \
                             ta1_x_x_xxyy_1,   \
                             ta1_x_x_xxyyy_0,  \
                             ta1_x_x_xxyyy_1,  \
                             ta1_x_x_xxyyz_0,  \
                             ta1_x_x_xxyyz_1,  \
                             ta1_x_x_xxyz_0,   \
                             ta1_x_x_xxyz_1,   \
                             ta1_x_x_xxyzz_0,  \
                             ta1_x_x_xxyzz_1,  \
                             ta1_x_x_xxzz_0,   \
                             ta1_x_x_xxzz_1,   \
                             ta1_x_x_xxzzz_0,  \
                             ta1_x_x_xxzzz_1,  \
                             ta1_x_x_xyyy_0,   \
                             ta1_x_x_xyyy_1,   \
                             ta1_x_x_xyyyy_0,  \
                             ta1_x_x_xyyyy_1,  \
                             ta1_x_x_xyyyz_0,  \
                             ta1_x_x_xyyyz_1,  \
                             ta1_x_x_xyyz_0,   \
                             ta1_x_x_xyyz_1,   \
                             ta1_x_x_xyyzz_0,  \
                             ta1_x_x_xyyzz_1,  \
                             ta1_x_x_xyzz_0,   \
                             ta1_x_x_xyzz_1,   \
                             ta1_x_x_xyzzz_0,  \
                             ta1_x_x_xyzzz_1,  \
                             ta1_x_x_xzzz_0,   \
                             ta1_x_x_xzzz_1,   \
                             ta1_x_x_xzzzz_0,  \
                             ta1_x_x_xzzzz_1,  \
                             ta1_x_x_yyyyy_0,  \
                             ta1_x_x_yyyyy_1,  \
                             ta1_x_xz_xxxxx_0, \
                             ta1_x_xz_xxxxy_0, \
                             ta1_x_xz_xxxxz_0, \
                             ta1_x_xz_xxxyy_0, \
                             ta1_x_xz_xxxyz_0, \
                             ta1_x_xz_xxxzz_0, \
                             ta1_x_xz_xxyyy_0, \
                             ta1_x_xz_xxyyz_0, \
                             ta1_x_xz_xxyzz_0, \
                             ta1_x_xz_xxzzz_0, \
                             ta1_x_xz_xyyyy_0, \
                             ta1_x_xz_xyyyz_0, \
                             ta1_x_xz_xyyzz_0, \
                             ta1_x_xz_xyzzz_0, \
                             ta1_x_xz_xzzzz_0, \
                             ta1_x_xz_yyyyy_0, \
                             ta1_x_xz_yyyyz_0, \
                             ta1_x_xz_yyyzz_0, \
                             ta1_x_xz_yyzzz_0, \
                             ta1_x_xz_yzzzz_0, \
                             ta1_x_xz_zzzzz_0, \
                             ta1_x_z_yyyyz_0,  \
                             ta1_x_z_yyyyz_1,  \
                             ta1_x_z_yyyzz_0,  \
                             ta1_x_z_yyyzz_1,  \
                             ta1_x_z_yyzzz_0,  \
                             ta1_x_z_yyzzz_1,  \
                             ta1_x_z_yzzzz_0,  \
                             ta1_x_z_yzzzz_1,  \
                             ta1_x_z_zzzzz_0,  \
                             ta1_x_z_zzzzz_1,  \
                             ta_z_yyyyz_1,     \
                             ta_z_yyyzz_1,     \
                             ta_z_yyzzz_1,     \
                             ta_z_yzzzz_1,     \
                             ta_z_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xz_xxxxx_0[i] = ta1_x_x_xxxxx_0[i] * pa_z[i] - ta1_x_x_xxxxx_1[i] * pc_z[i];

        ta1_x_xz_xxxxy_0[i] = ta1_x_x_xxxxy_0[i] * pa_z[i] - ta1_x_x_xxxxy_1[i] * pc_z[i];

        ta1_x_xz_xxxxz_0[i] = ta1_x_x_xxxx_0[i] * fe_0 - ta1_x_x_xxxx_1[i] * fe_0 + ta1_x_x_xxxxz_0[i] * pa_z[i] - ta1_x_x_xxxxz_1[i] * pc_z[i];

        ta1_x_xz_xxxyy_0[i] = ta1_x_x_xxxyy_0[i] * pa_z[i] - ta1_x_x_xxxyy_1[i] * pc_z[i];

        ta1_x_xz_xxxyz_0[i] = ta1_x_x_xxxy_0[i] * fe_0 - ta1_x_x_xxxy_1[i] * fe_0 + ta1_x_x_xxxyz_0[i] * pa_z[i] - ta1_x_x_xxxyz_1[i] * pc_z[i];

        ta1_x_xz_xxxzz_0[i] =
            2.0 * ta1_x_x_xxxz_0[i] * fe_0 - 2.0 * ta1_x_x_xxxz_1[i] * fe_0 + ta1_x_x_xxxzz_0[i] * pa_z[i] - ta1_x_x_xxxzz_1[i] * pc_z[i];

        ta1_x_xz_xxyyy_0[i] = ta1_x_x_xxyyy_0[i] * pa_z[i] - ta1_x_x_xxyyy_1[i] * pc_z[i];

        ta1_x_xz_xxyyz_0[i] = ta1_x_x_xxyy_0[i] * fe_0 - ta1_x_x_xxyy_1[i] * fe_0 + ta1_x_x_xxyyz_0[i] * pa_z[i] - ta1_x_x_xxyyz_1[i] * pc_z[i];

        ta1_x_xz_xxyzz_0[i] =
            2.0 * ta1_x_x_xxyz_0[i] * fe_0 - 2.0 * ta1_x_x_xxyz_1[i] * fe_0 + ta1_x_x_xxyzz_0[i] * pa_z[i] - ta1_x_x_xxyzz_1[i] * pc_z[i];

        ta1_x_xz_xxzzz_0[i] =
            3.0 * ta1_x_x_xxzz_0[i] * fe_0 - 3.0 * ta1_x_x_xxzz_1[i] * fe_0 + ta1_x_x_xxzzz_0[i] * pa_z[i] - ta1_x_x_xxzzz_1[i] * pc_z[i];

        ta1_x_xz_xyyyy_0[i] = ta1_x_x_xyyyy_0[i] * pa_z[i] - ta1_x_x_xyyyy_1[i] * pc_z[i];

        ta1_x_xz_xyyyz_0[i] = ta1_x_x_xyyy_0[i] * fe_0 - ta1_x_x_xyyy_1[i] * fe_0 + ta1_x_x_xyyyz_0[i] * pa_z[i] - ta1_x_x_xyyyz_1[i] * pc_z[i];

        ta1_x_xz_xyyzz_0[i] =
            2.0 * ta1_x_x_xyyz_0[i] * fe_0 - 2.0 * ta1_x_x_xyyz_1[i] * fe_0 + ta1_x_x_xyyzz_0[i] * pa_z[i] - ta1_x_x_xyyzz_1[i] * pc_z[i];

        ta1_x_xz_xyzzz_0[i] =
            3.0 * ta1_x_x_xyzz_0[i] * fe_0 - 3.0 * ta1_x_x_xyzz_1[i] * fe_0 + ta1_x_x_xyzzz_0[i] * pa_z[i] - ta1_x_x_xyzzz_1[i] * pc_z[i];

        ta1_x_xz_xzzzz_0[i] =
            4.0 * ta1_x_x_xzzz_0[i] * fe_0 - 4.0 * ta1_x_x_xzzz_1[i] * fe_0 + ta1_x_x_xzzzz_0[i] * pa_z[i] - ta1_x_x_xzzzz_1[i] * pc_z[i];

        ta1_x_xz_yyyyy_0[i] = ta1_x_x_yyyyy_0[i] * pa_z[i] - ta1_x_x_yyyyy_1[i] * pc_z[i];

        ta1_x_xz_yyyyz_0[i] = ta_z_yyyyz_1[i] + ta1_x_z_yyyyz_0[i] * pa_x[i] - ta1_x_z_yyyyz_1[i] * pc_x[i];

        ta1_x_xz_yyyzz_0[i] = ta_z_yyyzz_1[i] + ta1_x_z_yyyzz_0[i] * pa_x[i] - ta1_x_z_yyyzz_1[i] * pc_x[i];

        ta1_x_xz_yyzzz_0[i] = ta_z_yyzzz_1[i] + ta1_x_z_yyzzz_0[i] * pa_x[i] - ta1_x_z_yyzzz_1[i] * pc_x[i];

        ta1_x_xz_yzzzz_0[i] = ta_z_yzzzz_1[i] + ta1_x_z_yzzzz_0[i] * pa_x[i] - ta1_x_z_yzzzz_1[i] * pc_x[i];

        ta1_x_xz_zzzzz_0[i] = ta_z_zzzzz_1[i] + ta1_x_z_zzzzz_0[i] * pa_x[i] - ta1_x_z_zzzzz_1[i] * pc_x[i];
    }

    // Set up 63-84 components of targeted buffer : DH

    auto ta1_x_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 63);

    auto ta1_x_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 64);

    auto ta1_x_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 65);

    auto ta1_x_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 66);

    auto ta1_x_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 67);

    auto ta1_x_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 68);

    auto ta1_x_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 69);

    auto ta1_x_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 70);

    auto ta1_x_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 71);

    auto ta1_x_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 72);

    auto ta1_x_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 73);

    auto ta1_x_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 74);

    auto ta1_x_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 75);

    auto ta1_x_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 76);

    auto ta1_x_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 77);

    auto ta1_x_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 78);

    auto ta1_x_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 79);

    auto ta1_x_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 80);

    auto ta1_x_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 81);

    auto ta1_x_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 82);

    auto ta1_x_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 83);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_y_xxxx_0,   \
                             ta1_x_y_xxxx_1,   \
                             ta1_x_y_xxxxx_0,  \
                             ta1_x_y_xxxxx_1,  \
                             ta1_x_y_xxxxy_0,  \
                             ta1_x_y_xxxxy_1,  \
                             ta1_x_y_xxxxz_0,  \
                             ta1_x_y_xxxxz_1,  \
                             ta1_x_y_xxxy_0,   \
                             ta1_x_y_xxxy_1,   \
                             ta1_x_y_xxxyy_0,  \
                             ta1_x_y_xxxyy_1,  \
                             ta1_x_y_xxxyz_0,  \
                             ta1_x_y_xxxyz_1,  \
                             ta1_x_y_xxxz_0,   \
                             ta1_x_y_xxxz_1,   \
                             ta1_x_y_xxxzz_0,  \
                             ta1_x_y_xxxzz_1,  \
                             ta1_x_y_xxyy_0,   \
                             ta1_x_y_xxyy_1,   \
                             ta1_x_y_xxyyy_0,  \
                             ta1_x_y_xxyyy_1,  \
                             ta1_x_y_xxyyz_0,  \
                             ta1_x_y_xxyyz_1,  \
                             ta1_x_y_xxyz_0,   \
                             ta1_x_y_xxyz_1,   \
                             ta1_x_y_xxyzz_0,  \
                             ta1_x_y_xxyzz_1,  \
                             ta1_x_y_xxzz_0,   \
                             ta1_x_y_xxzz_1,   \
                             ta1_x_y_xxzzz_0,  \
                             ta1_x_y_xxzzz_1,  \
                             ta1_x_y_xyyy_0,   \
                             ta1_x_y_xyyy_1,   \
                             ta1_x_y_xyyyy_0,  \
                             ta1_x_y_xyyyy_1,  \
                             ta1_x_y_xyyyz_0,  \
                             ta1_x_y_xyyyz_1,  \
                             ta1_x_y_xyyz_0,   \
                             ta1_x_y_xyyz_1,   \
                             ta1_x_y_xyyzz_0,  \
                             ta1_x_y_xyyzz_1,  \
                             ta1_x_y_xyzz_0,   \
                             ta1_x_y_xyzz_1,   \
                             ta1_x_y_xyzzz_0,  \
                             ta1_x_y_xyzzz_1,  \
                             ta1_x_y_xzzz_0,   \
                             ta1_x_y_xzzz_1,   \
                             ta1_x_y_xzzzz_0,  \
                             ta1_x_y_xzzzz_1,  \
                             ta1_x_y_yyyy_0,   \
                             ta1_x_y_yyyy_1,   \
                             ta1_x_y_yyyyy_0,  \
                             ta1_x_y_yyyyy_1,  \
                             ta1_x_y_yyyyz_0,  \
                             ta1_x_y_yyyyz_1,  \
                             ta1_x_y_yyyz_0,   \
                             ta1_x_y_yyyz_1,   \
                             ta1_x_y_yyyzz_0,  \
                             ta1_x_y_yyyzz_1,  \
                             ta1_x_y_yyzz_0,   \
                             ta1_x_y_yyzz_1,   \
                             ta1_x_y_yyzzz_0,  \
                             ta1_x_y_yyzzz_1,  \
                             ta1_x_y_yzzz_0,   \
                             ta1_x_y_yzzz_1,   \
                             ta1_x_y_yzzzz_0,  \
                             ta1_x_y_yzzzz_1,  \
                             ta1_x_y_zzzz_0,   \
                             ta1_x_y_zzzz_1,   \
                             ta1_x_y_zzzzz_0,  \
                             ta1_x_y_zzzzz_1,  \
                             ta1_x_yy_xxxxx_0, \
                             ta1_x_yy_xxxxy_0, \
                             ta1_x_yy_xxxxz_0, \
                             ta1_x_yy_xxxyy_0, \
                             ta1_x_yy_xxxyz_0, \
                             ta1_x_yy_xxxzz_0, \
                             ta1_x_yy_xxyyy_0, \
                             ta1_x_yy_xxyyz_0, \
                             ta1_x_yy_xxyzz_0, \
                             ta1_x_yy_xxzzz_0, \
                             ta1_x_yy_xyyyy_0, \
                             ta1_x_yy_xyyyz_0, \
                             ta1_x_yy_xyyzz_0, \
                             ta1_x_yy_xyzzz_0, \
                             ta1_x_yy_xzzzz_0, \
                             ta1_x_yy_yyyyy_0, \
                             ta1_x_yy_yyyyz_0, \
                             ta1_x_yy_yyyzz_0, \
                             ta1_x_yy_yyzzz_0, \
                             ta1_x_yy_yzzzz_0, \
                             ta1_x_yy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yy_xxxxx_0[i] = ta1_x_0_xxxxx_0[i] * fe_0 - ta1_x_0_xxxxx_1[i] * fe_0 + ta1_x_y_xxxxx_0[i] * pa_y[i] - ta1_x_y_xxxxx_1[i] * pc_y[i];

        ta1_x_yy_xxxxy_0[i] = ta1_x_0_xxxxy_0[i] * fe_0 - ta1_x_0_xxxxy_1[i] * fe_0 + ta1_x_y_xxxx_0[i] * fe_0 - ta1_x_y_xxxx_1[i] * fe_0 +
                              ta1_x_y_xxxxy_0[i] * pa_y[i] - ta1_x_y_xxxxy_1[i] * pc_y[i];

        ta1_x_yy_xxxxz_0[i] = ta1_x_0_xxxxz_0[i] * fe_0 - ta1_x_0_xxxxz_1[i] * fe_0 + ta1_x_y_xxxxz_0[i] * pa_y[i] - ta1_x_y_xxxxz_1[i] * pc_y[i];

        ta1_x_yy_xxxyy_0[i] = ta1_x_0_xxxyy_0[i] * fe_0 - ta1_x_0_xxxyy_1[i] * fe_0 + 2.0 * ta1_x_y_xxxy_0[i] * fe_0 -
                              2.0 * ta1_x_y_xxxy_1[i] * fe_0 + ta1_x_y_xxxyy_0[i] * pa_y[i] - ta1_x_y_xxxyy_1[i] * pc_y[i];

        ta1_x_yy_xxxyz_0[i] = ta1_x_0_xxxyz_0[i] * fe_0 - ta1_x_0_xxxyz_1[i] * fe_0 + ta1_x_y_xxxz_0[i] * fe_0 - ta1_x_y_xxxz_1[i] * fe_0 +
                              ta1_x_y_xxxyz_0[i] * pa_y[i] - ta1_x_y_xxxyz_1[i] * pc_y[i];

        ta1_x_yy_xxxzz_0[i] = ta1_x_0_xxxzz_0[i] * fe_0 - ta1_x_0_xxxzz_1[i] * fe_0 + ta1_x_y_xxxzz_0[i] * pa_y[i] - ta1_x_y_xxxzz_1[i] * pc_y[i];

        ta1_x_yy_xxyyy_0[i] = ta1_x_0_xxyyy_0[i] * fe_0 - ta1_x_0_xxyyy_1[i] * fe_0 + 3.0 * ta1_x_y_xxyy_0[i] * fe_0 -
                              3.0 * ta1_x_y_xxyy_1[i] * fe_0 + ta1_x_y_xxyyy_0[i] * pa_y[i] - ta1_x_y_xxyyy_1[i] * pc_y[i];

        ta1_x_yy_xxyyz_0[i] = ta1_x_0_xxyyz_0[i] * fe_0 - ta1_x_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_y_xxyz_0[i] * fe_0 -
                              2.0 * ta1_x_y_xxyz_1[i] * fe_0 + ta1_x_y_xxyyz_0[i] * pa_y[i] - ta1_x_y_xxyyz_1[i] * pc_y[i];

        ta1_x_yy_xxyzz_0[i] = ta1_x_0_xxyzz_0[i] * fe_0 - ta1_x_0_xxyzz_1[i] * fe_0 + ta1_x_y_xxzz_0[i] * fe_0 - ta1_x_y_xxzz_1[i] * fe_0 +
                              ta1_x_y_xxyzz_0[i] * pa_y[i] - ta1_x_y_xxyzz_1[i] * pc_y[i];

        ta1_x_yy_xxzzz_0[i] = ta1_x_0_xxzzz_0[i] * fe_0 - ta1_x_0_xxzzz_1[i] * fe_0 + ta1_x_y_xxzzz_0[i] * pa_y[i] - ta1_x_y_xxzzz_1[i] * pc_y[i];

        ta1_x_yy_xyyyy_0[i] = ta1_x_0_xyyyy_0[i] * fe_0 - ta1_x_0_xyyyy_1[i] * fe_0 + 4.0 * ta1_x_y_xyyy_0[i] * fe_0 -
                              4.0 * ta1_x_y_xyyy_1[i] * fe_0 + ta1_x_y_xyyyy_0[i] * pa_y[i] - ta1_x_y_xyyyy_1[i] * pc_y[i];

        ta1_x_yy_xyyyz_0[i] = ta1_x_0_xyyyz_0[i] * fe_0 - ta1_x_0_xyyyz_1[i] * fe_0 + 3.0 * ta1_x_y_xyyz_0[i] * fe_0 -
                              3.0 * ta1_x_y_xyyz_1[i] * fe_0 + ta1_x_y_xyyyz_0[i] * pa_y[i] - ta1_x_y_xyyyz_1[i] * pc_y[i];

        ta1_x_yy_xyyzz_0[i] = ta1_x_0_xyyzz_0[i] * fe_0 - ta1_x_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_y_xyzz_0[i] * fe_0 -
                              2.0 * ta1_x_y_xyzz_1[i] * fe_0 + ta1_x_y_xyyzz_0[i] * pa_y[i] - ta1_x_y_xyyzz_1[i] * pc_y[i];

        ta1_x_yy_xyzzz_0[i] = ta1_x_0_xyzzz_0[i] * fe_0 - ta1_x_0_xyzzz_1[i] * fe_0 + ta1_x_y_xzzz_0[i] * fe_0 - ta1_x_y_xzzz_1[i] * fe_0 +
                              ta1_x_y_xyzzz_0[i] * pa_y[i] - ta1_x_y_xyzzz_1[i] * pc_y[i];

        ta1_x_yy_xzzzz_0[i] = ta1_x_0_xzzzz_0[i] * fe_0 - ta1_x_0_xzzzz_1[i] * fe_0 + ta1_x_y_xzzzz_0[i] * pa_y[i] - ta1_x_y_xzzzz_1[i] * pc_y[i];

        ta1_x_yy_yyyyy_0[i] = ta1_x_0_yyyyy_0[i] * fe_0 - ta1_x_0_yyyyy_1[i] * fe_0 + 5.0 * ta1_x_y_yyyy_0[i] * fe_0 -
                              5.0 * ta1_x_y_yyyy_1[i] * fe_0 + ta1_x_y_yyyyy_0[i] * pa_y[i] - ta1_x_y_yyyyy_1[i] * pc_y[i];

        ta1_x_yy_yyyyz_0[i] = ta1_x_0_yyyyz_0[i] * fe_0 - ta1_x_0_yyyyz_1[i] * fe_0 + 4.0 * ta1_x_y_yyyz_0[i] * fe_0 -
                              4.0 * ta1_x_y_yyyz_1[i] * fe_0 + ta1_x_y_yyyyz_0[i] * pa_y[i] - ta1_x_y_yyyyz_1[i] * pc_y[i];

        ta1_x_yy_yyyzz_0[i] = ta1_x_0_yyyzz_0[i] * fe_0 - ta1_x_0_yyyzz_1[i] * fe_0 + 3.0 * ta1_x_y_yyzz_0[i] * fe_0 -
                              3.0 * ta1_x_y_yyzz_1[i] * fe_0 + ta1_x_y_yyyzz_0[i] * pa_y[i] - ta1_x_y_yyyzz_1[i] * pc_y[i];

        ta1_x_yy_yyzzz_0[i] = ta1_x_0_yyzzz_0[i] * fe_0 - ta1_x_0_yyzzz_1[i] * fe_0 + 2.0 * ta1_x_y_yzzz_0[i] * fe_0 -
                              2.0 * ta1_x_y_yzzz_1[i] * fe_0 + ta1_x_y_yyzzz_0[i] * pa_y[i] - ta1_x_y_yyzzz_1[i] * pc_y[i];

        ta1_x_yy_yzzzz_0[i] = ta1_x_0_yzzzz_0[i] * fe_0 - ta1_x_0_yzzzz_1[i] * fe_0 + ta1_x_y_zzzz_0[i] * fe_0 - ta1_x_y_zzzz_1[i] * fe_0 +
                              ta1_x_y_yzzzz_0[i] * pa_y[i] - ta1_x_y_yzzzz_1[i] * pc_y[i];

        ta1_x_yy_zzzzz_0[i] = ta1_x_0_zzzzz_0[i] * fe_0 - ta1_x_0_zzzzz_1[i] * fe_0 + ta1_x_y_zzzzz_0[i] * pa_y[i] - ta1_x_y_zzzzz_1[i] * pc_y[i];
    }

    // Set up 84-105 components of targeted buffer : DH

    auto ta1_x_yz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 84);

    auto ta1_x_yz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 85);

    auto ta1_x_yz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 86);

    auto ta1_x_yz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 87);

    auto ta1_x_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 88);

    auto ta1_x_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 89);

    auto ta1_x_yz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 90);

    auto ta1_x_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 91);

    auto ta1_x_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 92);

    auto ta1_x_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 93);

    auto ta1_x_yz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 94);

    auto ta1_x_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 95);

    auto ta1_x_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 96);

    auto ta1_x_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 97);

    auto ta1_x_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 98);

    auto ta1_x_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 99);

    auto ta1_x_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 100);

    auto ta1_x_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 101);

    auto ta1_x_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 102);

    auto ta1_x_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 103);

    auto ta1_x_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 104);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_y_xxxxy_0,  \
                             ta1_x_y_xxxxy_1,  \
                             ta1_x_y_xxxyy_0,  \
                             ta1_x_y_xxxyy_1,  \
                             ta1_x_y_xxyyy_0,  \
                             ta1_x_y_xxyyy_1,  \
                             ta1_x_y_xyyyy_0,  \
                             ta1_x_y_xyyyy_1,  \
                             ta1_x_y_yyyyy_0,  \
                             ta1_x_y_yyyyy_1,  \
                             ta1_x_yz_xxxxx_0, \
                             ta1_x_yz_xxxxy_0, \
                             ta1_x_yz_xxxxz_0, \
                             ta1_x_yz_xxxyy_0, \
                             ta1_x_yz_xxxyz_0, \
                             ta1_x_yz_xxxzz_0, \
                             ta1_x_yz_xxyyy_0, \
                             ta1_x_yz_xxyyz_0, \
                             ta1_x_yz_xxyzz_0, \
                             ta1_x_yz_xxzzz_0, \
                             ta1_x_yz_xyyyy_0, \
                             ta1_x_yz_xyyyz_0, \
                             ta1_x_yz_xyyzz_0, \
                             ta1_x_yz_xyzzz_0, \
                             ta1_x_yz_xzzzz_0, \
                             ta1_x_yz_yyyyy_0, \
                             ta1_x_yz_yyyyz_0, \
                             ta1_x_yz_yyyzz_0, \
                             ta1_x_yz_yyzzz_0, \
                             ta1_x_yz_yzzzz_0, \
                             ta1_x_yz_zzzzz_0, \
                             ta1_x_z_xxxxx_0,  \
                             ta1_x_z_xxxxx_1,  \
                             ta1_x_z_xxxxz_0,  \
                             ta1_x_z_xxxxz_1,  \
                             ta1_x_z_xxxyz_0,  \
                             ta1_x_z_xxxyz_1,  \
                             ta1_x_z_xxxz_0,   \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxxzz_0,  \
                             ta1_x_z_xxxzz_1,  \
                             ta1_x_z_xxyyz_0,  \
                             ta1_x_z_xxyyz_1,  \
                             ta1_x_z_xxyz_0,   \
                             ta1_x_z_xxyz_1,   \
                             ta1_x_z_xxyzz_0,  \
                             ta1_x_z_xxyzz_1,  \
                             ta1_x_z_xxzz_0,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xxzzz_0,  \
                             ta1_x_z_xxzzz_1,  \
                             ta1_x_z_xyyyz_0,  \
                             ta1_x_z_xyyyz_1,  \
                             ta1_x_z_xyyz_0,   \
                             ta1_x_z_xyyz_1,   \
                             ta1_x_z_xyyzz_0,  \
                             ta1_x_z_xyyzz_1,  \
                             ta1_x_z_xyzz_0,   \
                             ta1_x_z_xyzz_1,   \
                             ta1_x_z_xyzzz_0,  \
                             ta1_x_z_xyzzz_1,  \
                             ta1_x_z_xzzz_0,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_xzzzz_0,  \
                             ta1_x_z_xzzzz_1,  \
                             ta1_x_z_yyyyz_0,  \
                             ta1_x_z_yyyyz_1,  \
                             ta1_x_z_yyyz_0,   \
                             ta1_x_z_yyyz_1,   \
                             ta1_x_z_yyyzz_0,  \
                             ta1_x_z_yyyzz_1,  \
                             ta1_x_z_yyzz_0,   \
                             ta1_x_z_yyzz_1,   \
                             ta1_x_z_yyzzz_0,  \
                             ta1_x_z_yyzzz_1,  \
                             ta1_x_z_yzzz_0,   \
                             ta1_x_z_yzzz_1,   \
                             ta1_x_z_yzzzz_0,  \
                             ta1_x_z_yzzzz_1,  \
                             ta1_x_z_zzzz_0,   \
                             ta1_x_z_zzzz_1,   \
                             ta1_x_z_zzzzz_0,  \
                             ta1_x_z_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yz_xxxxx_0[i] = ta1_x_z_xxxxx_0[i] * pa_y[i] - ta1_x_z_xxxxx_1[i] * pc_y[i];

        ta1_x_yz_xxxxy_0[i] = ta1_x_y_xxxxy_0[i] * pa_z[i] - ta1_x_y_xxxxy_1[i] * pc_z[i];

        ta1_x_yz_xxxxz_0[i] = ta1_x_z_xxxxz_0[i] * pa_y[i] - ta1_x_z_xxxxz_1[i] * pc_y[i];

        ta1_x_yz_xxxyy_0[i] = ta1_x_y_xxxyy_0[i] * pa_z[i] - ta1_x_y_xxxyy_1[i] * pc_z[i];

        ta1_x_yz_xxxyz_0[i] = ta1_x_z_xxxz_0[i] * fe_0 - ta1_x_z_xxxz_1[i] * fe_0 + ta1_x_z_xxxyz_0[i] * pa_y[i] - ta1_x_z_xxxyz_1[i] * pc_y[i];

        ta1_x_yz_xxxzz_0[i] = ta1_x_z_xxxzz_0[i] * pa_y[i] - ta1_x_z_xxxzz_1[i] * pc_y[i];

        ta1_x_yz_xxyyy_0[i] = ta1_x_y_xxyyy_0[i] * pa_z[i] - ta1_x_y_xxyyy_1[i] * pc_z[i];

        ta1_x_yz_xxyyz_0[i] =
            2.0 * ta1_x_z_xxyz_0[i] * fe_0 - 2.0 * ta1_x_z_xxyz_1[i] * fe_0 + ta1_x_z_xxyyz_0[i] * pa_y[i] - ta1_x_z_xxyyz_1[i] * pc_y[i];

        ta1_x_yz_xxyzz_0[i] = ta1_x_z_xxzz_0[i] * fe_0 - ta1_x_z_xxzz_1[i] * fe_0 + ta1_x_z_xxyzz_0[i] * pa_y[i] - ta1_x_z_xxyzz_1[i] * pc_y[i];

        ta1_x_yz_xxzzz_0[i] = ta1_x_z_xxzzz_0[i] * pa_y[i] - ta1_x_z_xxzzz_1[i] * pc_y[i];

        ta1_x_yz_xyyyy_0[i] = ta1_x_y_xyyyy_0[i] * pa_z[i] - ta1_x_y_xyyyy_1[i] * pc_z[i];

        ta1_x_yz_xyyyz_0[i] =
            3.0 * ta1_x_z_xyyz_0[i] * fe_0 - 3.0 * ta1_x_z_xyyz_1[i] * fe_0 + ta1_x_z_xyyyz_0[i] * pa_y[i] - ta1_x_z_xyyyz_1[i] * pc_y[i];

        ta1_x_yz_xyyzz_0[i] =
            2.0 * ta1_x_z_xyzz_0[i] * fe_0 - 2.0 * ta1_x_z_xyzz_1[i] * fe_0 + ta1_x_z_xyyzz_0[i] * pa_y[i] - ta1_x_z_xyyzz_1[i] * pc_y[i];

        ta1_x_yz_xyzzz_0[i] = ta1_x_z_xzzz_0[i] * fe_0 - ta1_x_z_xzzz_1[i] * fe_0 + ta1_x_z_xyzzz_0[i] * pa_y[i] - ta1_x_z_xyzzz_1[i] * pc_y[i];

        ta1_x_yz_xzzzz_0[i] = ta1_x_z_xzzzz_0[i] * pa_y[i] - ta1_x_z_xzzzz_1[i] * pc_y[i];

        ta1_x_yz_yyyyy_0[i] = ta1_x_y_yyyyy_0[i] * pa_z[i] - ta1_x_y_yyyyy_1[i] * pc_z[i];

        ta1_x_yz_yyyyz_0[i] =
            4.0 * ta1_x_z_yyyz_0[i] * fe_0 - 4.0 * ta1_x_z_yyyz_1[i] * fe_0 + ta1_x_z_yyyyz_0[i] * pa_y[i] - ta1_x_z_yyyyz_1[i] * pc_y[i];

        ta1_x_yz_yyyzz_0[i] =
            3.0 * ta1_x_z_yyzz_0[i] * fe_0 - 3.0 * ta1_x_z_yyzz_1[i] * fe_0 + ta1_x_z_yyyzz_0[i] * pa_y[i] - ta1_x_z_yyyzz_1[i] * pc_y[i];

        ta1_x_yz_yyzzz_0[i] =
            2.0 * ta1_x_z_yzzz_0[i] * fe_0 - 2.0 * ta1_x_z_yzzz_1[i] * fe_0 + ta1_x_z_yyzzz_0[i] * pa_y[i] - ta1_x_z_yyzzz_1[i] * pc_y[i];

        ta1_x_yz_yzzzz_0[i] = ta1_x_z_zzzz_0[i] * fe_0 - ta1_x_z_zzzz_1[i] * fe_0 + ta1_x_z_yzzzz_0[i] * pa_y[i] - ta1_x_z_yzzzz_1[i] * pc_y[i];

        ta1_x_yz_zzzzz_0[i] = ta1_x_z_zzzzz_0[i] * pa_y[i] - ta1_x_z_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : DH

    auto ta1_x_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 105);

    auto ta1_x_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 106);

    auto ta1_x_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 107);

    auto ta1_x_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 108);

    auto ta1_x_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 109);

    auto ta1_x_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 110);

    auto ta1_x_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 111);

    auto ta1_x_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 112);

    auto ta1_x_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 113);

    auto ta1_x_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 114);

    auto ta1_x_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 115);

    auto ta1_x_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 116);

    auto ta1_x_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 117);

    auto ta1_x_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 118);

    auto ta1_x_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 119);

    auto ta1_x_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 120);

    auto ta1_x_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 121);

    auto ta1_x_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 122);

    auto ta1_x_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 123);

    auto ta1_x_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 124);

    auto ta1_x_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 125);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyz_0,  \
                             ta1_x_0_xxxyz_1,  \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyz_0,  \
                             ta1_x_0_xxyyz_1,  \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyz_0,  \
                             ta1_x_0_xyyyz_1,  \
                             ta1_x_0_xyyzz_0,  \
                             ta1_x_0_xyyzz_1,  \
                             ta1_x_0_xyzzz_0,  \
                             ta1_x_0_xyzzz_1,  \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyz_0,  \
                             ta1_x_0_yyyyz_1,  \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_z_xxxx_0,   \
                             ta1_x_z_xxxx_1,   \
                             ta1_x_z_xxxxx_0,  \
                             ta1_x_z_xxxxx_1,  \
                             ta1_x_z_xxxxy_0,  \
                             ta1_x_z_xxxxy_1,  \
                             ta1_x_z_xxxxz_0,  \
                             ta1_x_z_xxxxz_1,  \
                             ta1_x_z_xxxy_0,   \
                             ta1_x_z_xxxy_1,   \
                             ta1_x_z_xxxyy_0,  \
                             ta1_x_z_xxxyy_1,  \
                             ta1_x_z_xxxyz_0,  \
                             ta1_x_z_xxxyz_1,  \
                             ta1_x_z_xxxz_0,   \
                             ta1_x_z_xxxz_1,   \
                             ta1_x_z_xxxzz_0,  \
                             ta1_x_z_xxxzz_1,  \
                             ta1_x_z_xxyy_0,   \
                             ta1_x_z_xxyy_1,   \
                             ta1_x_z_xxyyy_0,  \
                             ta1_x_z_xxyyy_1,  \
                             ta1_x_z_xxyyz_0,  \
                             ta1_x_z_xxyyz_1,  \
                             ta1_x_z_xxyz_0,   \
                             ta1_x_z_xxyz_1,   \
                             ta1_x_z_xxyzz_0,  \
                             ta1_x_z_xxyzz_1,  \
                             ta1_x_z_xxzz_0,   \
                             ta1_x_z_xxzz_1,   \
                             ta1_x_z_xxzzz_0,  \
                             ta1_x_z_xxzzz_1,  \
                             ta1_x_z_xyyy_0,   \
                             ta1_x_z_xyyy_1,   \
                             ta1_x_z_xyyyy_0,  \
                             ta1_x_z_xyyyy_1,  \
                             ta1_x_z_xyyyz_0,  \
                             ta1_x_z_xyyyz_1,  \
                             ta1_x_z_xyyz_0,   \
                             ta1_x_z_xyyz_1,   \
                             ta1_x_z_xyyzz_0,  \
                             ta1_x_z_xyyzz_1,  \
                             ta1_x_z_xyzz_0,   \
                             ta1_x_z_xyzz_1,   \
                             ta1_x_z_xyzzz_0,  \
                             ta1_x_z_xyzzz_1,  \
                             ta1_x_z_xzzz_0,   \
                             ta1_x_z_xzzz_1,   \
                             ta1_x_z_xzzzz_0,  \
                             ta1_x_z_xzzzz_1,  \
                             ta1_x_z_yyyy_0,   \
                             ta1_x_z_yyyy_1,   \
                             ta1_x_z_yyyyy_0,  \
                             ta1_x_z_yyyyy_1,  \
                             ta1_x_z_yyyyz_0,  \
                             ta1_x_z_yyyyz_1,  \
                             ta1_x_z_yyyz_0,   \
                             ta1_x_z_yyyz_1,   \
                             ta1_x_z_yyyzz_0,  \
                             ta1_x_z_yyyzz_1,  \
                             ta1_x_z_yyzz_0,   \
                             ta1_x_z_yyzz_1,   \
                             ta1_x_z_yyzzz_0,  \
                             ta1_x_z_yyzzz_1,  \
                             ta1_x_z_yzzz_0,   \
                             ta1_x_z_yzzz_1,   \
                             ta1_x_z_yzzzz_0,  \
                             ta1_x_z_yzzzz_1,  \
                             ta1_x_z_zzzz_0,   \
                             ta1_x_z_zzzz_1,   \
                             ta1_x_z_zzzzz_0,  \
                             ta1_x_z_zzzzz_1,  \
                             ta1_x_zz_xxxxx_0, \
                             ta1_x_zz_xxxxy_0, \
                             ta1_x_zz_xxxxz_0, \
                             ta1_x_zz_xxxyy_0, \
                             ta1_x_zz_xxxyz_0, \
                             ta1_x_zz_xxxzz_0, \
                             ta1_x_zz_xxyyy_0, \
                             ta1_x_zz_xxyyz_0, \
                             ta1_x_zz_xxyzz_0, \
                             ta1_x_zz_xxzzz_0, \
                             ta1_x_zz_xyyyy_0, \
                             ta1_x_zz_xyyyz_0, \
                             ta1_x_zz_xyyzz_0, \
                             ta1_x_zz_xyzzz_0, \
                             ta1_x_zz_xzzzz_0, \
                             ta1_x_zz_yyyyy_0, \
                             ta1_x_zz_yyyyz_0, \
                             ta1_x_zz_yyyzz_0, \
                             ta1_x_zz_yyzzz_0, \
                             ta1_x_zz_yzzzz_0, \
                             ta1_x_zz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zz_xxxxx_0[i] = ta1_x_0_xxxxx_0[i] * fe_0 - ta1_x_0_xxxxx_1[i] * fe_0 + ta1_x_z_xxxxx_0[i] * pa_z[i] - ta1_x_z_xxxxx_1[i] * pc_z[i];

        ta1_x_zz_xxxxy_0[i] = ta1_x_0_xxxxy_0[i] * fe_0 - ta1_x_0_xxxxy_1[i] * fe_0 + ta1_x_z_xxxxy_0[i] * pa_z[i] - ta1_x_z_xxxxy_1[i] * pc_z[i];

        ta1_x_zz_xxxxz_0[i] = ta1_x_0_xxxxz_0[i] * fe_0 - ta1_x_0_xxxxz_1[i] * fe_0 + ta1_x_z_xxxx_0[i] * fe_0 - ta1_x_z_xxxx_1[i] * fe_0 +
                              ta1_x_z_xxxxz_0[i] * pa_z[i] - ta1_x_z_xxxxz_1[i] * pc_z[i];

        ta1_x_zz_xxxyy_0[i] = ta1_x_0_xxxyy_0[i] * fe_0 - ta1_x_0_xxxyy_1[i] * fe_0 + ta1_x_z_xxxyy_0[i] * pa_z[i] - ta1_x_z_xxxyy_1[i] * pc_z[i];

        ta1_x_zz_xxxyz_0[i] = ta1_x_0_xxxyz_0[i] * fe_0 - ta1_x_0_xxxyz_1[i] * fe_0 + ta1_x_z_xxxy_0[i] * fe_0 - ta1_x_z_xxxy_1[i] * fe_0 +
                              ta1_x_z_xxxyz_0[i] * pa_z[i] - ta1_x_z_xxxyz_1[i] * pc_z[i];

        ta1_x_zz_xxxzz_0[i] = ta1_x_0_xxxzz_0[i] * fe_0 - ta1_x_0_xxxzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxxz_0[i] * fe_0 -
                              2.0 * ta1_x_z_xxxz_1[i] * fe_0 + ta1_x_z_xxxzz_0[i] * pa_z[i] - ta1_x_z_xxxzz_1[i] * pc_z[i];

        ta1_x_zz_xxyyy_0[i] = ta1_x_0_xxyyy_0[i] * fe_0 - ta1_x_0_xxyyy_1[i] * fe_0 + ta1_x_z_xxyyy_0[i] * pa_z[i] - ta1_x_z_xxyyy_1[i] * pc_z[i];

        ta1_x_zz_xxyyz_0[i] = ta1_x_0_xxyyz_0[i] * fe_0 - ta1_x_0_xxyyz_1[i] * fe_0 + ta1_x_z_xxyy_0[i] * fe_0 - ta1_x_z_xxyy_1[i] * fe_0 +
                              ta1_x_z_xxyyz_0[i] * pa_z[i] - ta1_x_z_xxyyz_1[i] * pc_z[i];

        ta1_x_zz_xxyzz_0[i] = ta1_x_0_xxyzz_0[i] * fe_0 - ta1_x_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xxyz_0[i] * fe_0 -
                              2.0 * ta1_x_z_xxyz_1[i] * fe_0 + ta1_x_z_xxyzz_0[i] * pa_z[i] - ta1_x_z_xxyzz_1[i] * pc_z[i];

        ta1_x_zz_xxzzz_0[i] = ta1_x_0_xxzzz_0[i] * fe_0 - ta1_x_0_xxzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xxzz_0[i] * fe_0 -
                              3.0 * ta1_x_z_xxzz_1[i] * fe_0 + ta1_x_z_xxzzz_0[i] * pa_z[i] - ta1_x_z_xxzzz_1[i] * pc_z[i];

        ta1_x_zz_xyyyy_0[i] = ta1_x_0_xyyyy_0[i] * fe_0 - ta1_x_0_xyyyy_1[i] * fe_0 + ta1_x_z_xyyyy_0[i] * pa_z[i] - ta1_x_z_xyyyy_1[i] * pc_z[i];

        ta1_x_zz_xyyyz_0[i] = ta1_x_0_xyyyz_0[i] * fe_0 - ta1_x_0_xyyyz_1[i] * fe_0 + ta1_x_z_xyyy_0[i] * fe_0 - ta1_x_z_xyyy_1[i] * fe_0 +
                              ta1_x_z_xyyyz_0[i] * pa_z[i] - ta1_x_z_xyyyz_1[i] * pc_z[i];

        ta1_x_zz_xyyzz_0[i] = ta1_x_0_xyyzz_0[i] * fe_0 - ta1_x_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_z_xyyz_0[i] * fe_0 -
                              2.0 * ta1_x_z_xyyz_1[i] * fe_0 + ta1_x_z_xyyzz_0[i] * pa_z[i] - ta1_x_z_xyyzz_1[i] * pc_z[i];

        ta1_x_zz_xyzzz_0[i] = ta1_x_0_xyzzz_0[i] * fe_0 - ta1_x_0_xyzzz_1[i] * fe_0 + 3.0 * ta1_x_z_xyzz_0[i] * fe_0 -
                              3.0 * ta1_x_z_xyzz_1[i] * fe_0 + ta1_x_z_xyzzz_0[i] * pa_z[i] - ta1_x_z_xyzzz_1[i] * pc_z[i];

        ta1_x_zz_xzzzz_0[i] = ta1_x_0_xzzzz_0[i] * fe_0 - ta1_x_0_xzzzz_1[i] * fe_0 + 4.0 * ta1_x_z_xzzz_0[i] * fe_0 -
                              4.0 * ta1_x_z_xzzz_1[i] * fe_0 + ta1_x_z_xzzzz_0[i] * pa_z[i] - ta1_x_z_xzzzz_1[i] * pc_z[i];

        ta1_x_zz_yyyyy_0[i] = ta1_x_0_yyyyy_0[i] * fe_0 - ta1_x_0_yyyyy_1[i] * fe_0 + ta1_x_z_yyyyy_0[i] * pa_z[i] - ta1_x_z_yyyyy_1[i] * pc_z[i];

        ta1_x_zz_yyyyz_0[i] = ta1_x_0_yyyyz_0[i] * fe_0 - ta1_x_0_yyyyz_1[i] * fe_0 + ta1_x_z_yyyy_0[i] * fe_0 - ta1_x_z_yyyy_1[i] * fe_0 +
                              ta1_x_z_yyyyz_0[i] * pa_z[i] - ta1_x_z_yyyyz_1[i] * pc_z[i];

        ta1_x_zz_yyyzz_0[i] = ta1_x_0_yyyzz_0[i] * fe_0 - ta1_x_0_yyyzz_1[i] * fe_0 + 2.0 * ta1_x_z_yyyz_0[i] * fe_0 -
                              2.0 * ta1_x_z_yyyz_1[i] * fe_0 + ta1_x_z_yyyzz_0[i] * pa_z[i] - ta1_x_z_yyyzz_1[i] * pc_z[i];

        ta1_x_zz_yyzzz_0[i] = ta1_x_0_yyzzz_0[i] * fe_0 - ta1_x_0_yyzzz_1[i] * fe_0 + 3.0 * ta1_x_z_yyzz_0[i] * fe_0 -
                              3.0 * ta1_x_z_yyzz_1[i] * fe_0 + ta1_x_z_yyzzz_0[i] * pa_z[i] - ta1_x_z_yyzzz_1[i] * pc_z[i];

        ta1_x_zz_yzzzz_0[i] = ta1_x_0_yzzzz_0[i] * fe_0 - ta1_x_0_yzzzz_1[i] * fe_0 + 4.0 * ta1_x_z_yzzz_0[i] * fe_0 -
                              4.0 * ta1_x_z_yzzz_1[i] * fe_0 + ta1_x_z_yzzzz_0[i] * pa_z[i] - ta1_x_z_yzzzz_1[i] * pc_z[i];

        ta1_x_zz_zzzzz_0[i] = ta1_x_0_zzzzz_0[i] * fe_0 - ta1_x_0_zzzzz_1[i] * fe_0 + 5.0 * ta1_x_z_zzzz_0[i] * fe_0 -
                              5.0 * ta1_x_z_zzzz_1[i] * fe_0 + ta1_x_z_zzzzz_0[i] * pa_z[i] - ta1_x_z_zzzzz_1[i] * pc_z[i];
    }

    // Set up 126-147 components of targeted buffer : DH

    auto ta1_y_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 126);

    auto ta1_y_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 127);

    auto ta1_y_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 128);

    auto ta1_y_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 129);

    auto ta1_y_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 130);

    auto ta1_y_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 131);

    auto ta1_y_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 132);

    auto ta1_y_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 133);

    auto ta1_y_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 134);

    auto ta1_y_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 135);

    auto ta1_y_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 136);

    auto ta1_y_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 137);

    auto ta1_y_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 138);

    auto ta1_y_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 139);

    auto ta1_y_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 140);

    auto ta1_y_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 141);

    auto ta1_y_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 142);

    auto ta1_y_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 143);

    auto ta1_y_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 144);

    auto ta1_y_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 145);

    auto ta1_y_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 146);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_x_xxxx_0,   \
                             ta1_y_x_xxxx_1,   \
                             ta1_y_x_xxxxx_0,  \
                             ta1_y_x_xxxxx_1,  \
                             ta1_y_x_xxxxy_0,  \
                             ta1_y_x_xxxxy_1,  \
                             ta1_y_x_xxxxz_0,  \
                             ta1_y_x_xxxxz_1,  \
                             ta1_y_x_xxxy_0,   \
                             ta1_y_x_xxxy_1,   \
                             ta1_y_x_xxxyy_0,  \
                             ta1_y_x_xxxyy_1,  \
                             ta1_y_x_xxxyz_0,  \
                             ta1_y_x_xxxyz_1,  \
                             ta1_y_x_xxxz_0,   \
                             ta1_y_x_xxxz_1,   \
                             ta1_y_x_xxxzz_0,  \
                             ta1_y_x_xxxzz_1,  \
                             ta1_y_x_xxyy_0,   \
                             ta1_y_x_xxyy_1,   \
                             ta1_y_x_xxyyy_0,  \
                             ta1_y_x_xxyyy_1,  \
                             ta1_y_x_xxyyz_0,  \
                             ta1_y_x_xxyyz_1,  \
                             ta1_y_x_xxyz_0,   \
                             ta1_y_x_xxyz_1,   \
                             ta1_y_x_xxyzz_0,  \
                             ta1_y_x_xxyzz_1,  \
                             ta1_y_x_xxzz_0,   \
                             ta1_y_x_xxzz_1,   \
                             ta1_y_x_xxzzz_0,  \
                             ta1_y_x_xxzzz_1,  \
                             ta1_y_x_xyyy_0,   \
                             ta1_y_x_xyyy_1,   \
                             ta1_y_x_xyyyy_0,  \
                             ta1_y_x_xyyyy_1,  \
                             ta1_y_x_xyyyz_0,  \
                             ta1_y_x_xyyyz_1,  \
                             ta1_y_x_xyyz_0,   \
                             ta1_y_x_xyyz_1,   \
                             ta1_y_x_xyyzz_0,  \
                             ta1_y_x_xyyzz_1,  \
                             ta1_y_x_xyzz_0,   \
                             ta1_y_x_xyzz_1,   \
                             ta1_y_x_xyzzz_0,  \
                             ta1_y_x_xyzzz_1,  \
                             ta1_y_x_xzzz_0,   \
                             ta1_y_x_xzzz_1,   \
                             ta1_y_x_xzzzz_0,  \
                             ta1_y_x_xzzzz_1,  \
                             ta1_y_x_yyyy_0,   \
                             ta1_y_x_yyyy_1,   \
                             ta1_y_x_yyyyy_0,  \
                             ta1_y_x_yyyyy_1,  \
                             ta1_y_x_yyyyz_0,  \
                             ta1_y_x_yyyyz_1,  \
                             ta1_y_x_yyyz_0,   \
                             ta1_y_x_yyyz_1,   \
                             ta1_y_x_yyyzz_0,  \
                             ta1_y_x_yyyzz_1,  \
                             ta1_y_x_yyzz_0,   \
                             ta1_y_x_yyzz_1,   \
                             ta1_y_x_yyzzz_0,  \
                             ta1_y_x_yyzzz_1,  \
                             ta1_y_x_yzzz_0,   \
                             ta1_y_x_yzzz_1,   \
                             ta1_y_x_yzzzz_0,  \
                             ta1_y_x_yzzzz_1,  \
                             ta1_y_x_zzzz_0,   \
                             ta1_y_x_zzzz_1,   \
                             ta1_y_x_zzzzz_0,  \
                             ta1_y_x_zzzzz_1,  \
                             ta1_y_xx_xxxxx_0, \
                             ta1_y_xx_xxxxy_0, \
                             ta1_y_xx_xxxxz_0, \
                             ta1_y_xx_xxxyy_0, \
                             ta1_y_xx_xxxyz_0, \
                             ta1_y_xx_xxxzz_0, \
                             ta1_y_xx_xxyyy_0, \
                             ta1_y_xx_xxyyz_0, \
                             ta1_y_xx_xxyzz_0, \
                             ta1_y_xx_xxzzz_0, \
                             ta1_y_xx_xyyyy_0, \
                             ta1_y_xx_xyyyz_0, \
                             ta1_y_xx_xyyzz_0, \
                             ta1_y_xx_xyzzz_0, \
                             ta1_y_xx_xzzzz_0, \
                             ta1_y_xx_yyyyy_0, \
                             ta1_y_xx_yyyyz_0, \
                             ta1_y_xx_yyyzz_0, \
                             ta1_y_xx_yyzzz_0, \
                             ta1_y_xx_yzzzz_0, \
                             ta1_y_xx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xx_xxxxx_0[i] = ta1_y_0_xxxxx_0[i] * fe_0 - ta1_y_0_xxxxx_1[i] * fe_0 + 5.0 * ta1_y_x_xxxx_0[i] * fe_0 -
                              5.0 * ta1_y_x_xxxx_1[i] * fe_0 + ta1_y_x_xxxxx_0[i] * pa_x[i] - ta1_y_x_xxxxx_1[i] * pc_x[i];

        ta1_y_xx_xxxxy_0[i] = ta1_y_0_xxxxy_0[i] * fe_0 - ta1_y_0_xxxxy_1[i] * fe_0 + 4.0 * ta1_y_x_xxxy_0[i] * fe_0 -
                              4.0 * ta1_y_x_xxxy_1[i] * fe_0 + ta1_y_x_xxxxy_0[i] * pa_x[i] - ta1_y_x_xxxxy_1[i] * pc_x[i];

        ta1_y_xx_xxxxz_0[i] = ta1_y_0_xxxxz_0[i] * fe_0 - ta1_y_0_xxxxz_1[i] * fe_0 + 4.0 * ta1_y_x_xxxz_0[i] * fe_0 -
                              4.0 * ta1_y_x_xxxz_1[i] * fe_0 + ta1_y_x_xxxxz_0[i] * pa_x[i] - ta1_y_x_xxxxz_1[i] * pc_x[i];

        ta1_y_xx_xxxyy_0[i] = ta1_y_0_xxxyy_0[i] * fe_0 - ta1_y_0_xxxyy_1[i] * fe_0 + 3.0 * ta1_y_x_xxyy_0[i] * fe_0 -
                              3.0 * ta1_y_x_xxyy_1[i] * fe_0 + ta1_y_x_xxxyy_0[i] * pa_x[i] - ta1_y_x_xxxyy_1[i] * pc_x[i];

        ta1_y_xx_xxxyz_0[i] = ta1_y_0_xxxyz_0[i] * fe_0 - ta1_y_0_xxxyz_1[i] * fe_0 + 3.0 * ta1_y_x_xxyz_0[i] * fe_0 -
                              3.0 * ta1_y_x_xxyz_1[i] * fe_0 + ta1_y_x_xxxyz_0[i] * pa_x[i] - ta1_y_x_xxxyz_1[i] * pc_x[i];

        ta1_y_xx_xxxzz_0[i] = ta1_y_0_xxxzz_0[i] * fe_0 - ta1_y_0_xxxzz_1[i] * fe_0 + 3.0 * ta1_y_x_xxzz_0[i] * fe_0 -
                              3.0 * ta1_y_x_xxzz_1[i] * fe_0 + ta1_y_x_xxxzz_0[i] * pa_x[i] - ta1_y_x_xxxzz_1[i] * pc_x[i];

        ta1_y_xx_xxyyy_0[i] = ta1_y_0_xxyyy_0[i] * fe_0 - ta1_y_0_xxyyy_1[i] * fe_0 + 2.0 * ta1_y_x_xyyy_0[i] * fe_0 -
                              2.0 * ta1_y_x_xyyy_1[i] * fe_0 + ta1_y_x_xxyyy_0[i] * pa_x[i] - ta1_y_x_xxyyy_1[i] * pc_x[i];

        ta1_y_xx_xxyyz_0[i] = ta1_y_0_xxyyz_0[i] * fe_0 - ta1_y_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_x_xyyz_0[i] * fe_0 -
                              2.0 * ta1_y_x_xyyz_1[i] * fe_0 + ta1_y_x_xxyyz_0[i] * pa_x[i] - ta1_y_x_xxyyz_1[i] * pc_x[i];

        ta1_y_xx_xxyzz_0[i] = ta1_y_0_xxyzz_0[i] * fe_0 - ta1_y_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_x_xyzz_0[i] * fe_0 -
                              2.0 * ta1_y_x_xyzz_1[i] * fe_0 + ta1_y_x_xxyzz_0[i] * pa_x[i] - ta1_y_x_xxyzz_1[i] * pc_x[i];

        ta1_y_xx_xxzzz_0[i] = ta1_y_0_xxzzz_0[i] * fe_0 - ta1_y_0_xxzzz_1[i] * fe_0 + 2.0 * ta1_y_x_xzzz_0[i] * fe_0 -
                              2.0 * ta1_y_x_xzzz_1[i] * fe_0 + ta1_y_x_xxzzz_0[i] * pa_x[i] - ta1_y_x_xxzzz_1[i] * pc_x[i];

        ta1_y_xx_xyyyy_0[i] = ta1_y_0_xyyyy_0[i] * fe_0 - ta1_y_0_xyyyy_1[i] * fe_0 + ta1_y_x_yyyy_0[i] * fe_0 - ta1_y_x_yyyy_1[i] * fe_0 +
                              ta1_y_x_xyyyy_0[i] * pa_x[i] - ta1_y_x_xyyyy_1[i] * pc_x[i];

        ta1_y_xx_xyyyz_0[i] = ta1_y_0_xyyyz_0[i] * fe_0 - ta1_y_0_xyyyz_1[i] * fe_0 + ta1_y_x_yyyz_0[i] * fe_0 - ta1_y_x_yyyz_1[i] * fe_0 +
                              ta1_y_x_xyyyz_0[i] * pa_x[i] - ta1_y_x_xyyyz_1[i] * pc_x[i];

        ta1_y_xx_xyyzz_0[i] = ta1_y_0_xyyzz_0[i] * fe_0 - ta1_y_0_xyyzz_1[i] * fe_0 + ta1_y_x_yyzz_0[i] * fe_0 - ta1_y_x_yyzz_1[i] * fe_0 +
                              ta1_y_x_xyyzz_0[i] * pa_x[i] - ta1_y_x_xyyzz_1[i] * pc_x[i];

        ta1_y_xx_xyzzz_0[i] = ta1_y_0_xyzzz_0[i] * fe_0 - ta1_y_0_xyzzz_1[i] * fe_0 + ta1_y_x_yzzz_0[i] * fe_0 - ta1_y_x_yzzz_1[i] * fe_0 +
                              ta1_y_x_xyzzz_0[i] * pa_x[i] - ta1_y_x_xyzzz_1[i] * pc_x[i];

        ta1_y_xx_xzzzz_0[i] = ta1_y_0_xzzzz_0[i] * fe_0 - ta1_y_0_xzzzz_1[i] * fe_0 + ta1_y_x_zzzz_0[i] * fe_0 - ta1_y_x_zzzz_1[i] * fe_0 +
                              ta1_y_x_xzzzz_0[i] * pa_x[i] - ta1_y_x_xzzzz_1[i] * pc_x[i];

        ta1_y_xx_yyyyy_0[i] = ta1_y_0_yyyyy_0[i] * fe_0 - ta1_y_0_yyyyy_1[i] * fe_0 + ta1_y_x_yyyyy_0[i] * pa_x[i] - ta1_y_x_yyyyy_1[i] * pc_x[i];

        ta1_y_xx_yyyyz_0[i] = ta1_y_0_yyyyz_0[i] * fe_0 - ta1_y_0_yyyyz_1[i] * fe_0 + ta1_y_x_yyyyz_0[i] * pa_x[i] - ta1_y_x_yyyyz_1[i] * pc_x[i];

        ta1_y_xx_yyyzz_0[i] = ta1_y_0_yyyzz_0[i] * fe_0 - ta1_y_0_yyyzz_1[i] * fe_0 + ta1_y_x_yyyzz_0[i] * pa_x[i] - ta1_y_x_yyyzz_1[i] * pc_x[i];

        ta1_y_xx_yyzzz_0[i] = ta1_y_0_yyzzz_0[i] * fe_0 - ta1_y_0_yyzzz_1[i] * fe_0 + ta1_y_x_yyzzz_0[i] * pa_x[i] - ta1_y_x_yyzzz_1[i] * pc_x[i];

        ta1_y_xx_yzzzz_0[i] = ta1_y_0_yzzzz_0[i] * fe_0 - ta1_y_0_yzzzz_1[i] * fe_0 + ta1_y_x_yzzzz_0[i] * pa_x[i] - ta1_y_x_yzzzz_1[i] * pc_x[i];

        ta1_y_xx_zzzzz_0[i] = ta1_y_0_zzzzz_0[i] * fe_0 - ta1_y_0_zzzzz_1[i] * fe_0 + ta1_y_x_zzzzz_0[i] * pa_x[i] - ta1_y_x_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : DH

    auto ta1_y_xy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 147);

    auto ta1_y_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 148);

    auto ta1_y_xy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 149);

    auto ta1_y_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 150);

    auto ta1_y_xy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 151);

    auto ta1_y_xy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 152);

    auto ta1_y_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 153);

    auto ta1_y_xy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 154);

    auto ta1_y_xy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 155);

    auto ta1_y_xy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 156);

    auto ta1_y_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 157);

    auto ta1_y_xy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 158);

    auto ta1_y_xy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 159);

    auto ta1_y_xy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 160);

    auto ta1_y_xy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 161);

    auto ta1_y_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 162);

    auto ta1_y_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 163);

    auto ta1_y_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 164);

    auto ta1_y_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 165);

    auto ta1_y_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 166);

    auto ta1_y_xy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 167);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_x_xxxxx_0,  \
                             ta1_y_x_xxxxx_1,  \
                             ta1_y_x_xxxxz_0,  \
                             ta1_y_x_xxxxz_1,  \
                             ta1_y_x_xxxzz_0,  \
                             ta1_y_x_xxxzz_1,  \
                             ta1_y_x_xxzzz_0,  \
                             ta1_y_x_xxzzz_1,  \
                             ta1_y_x_xzzzz_0,  \
                             ta1_y_x_xzzzz_1,  \
                             ta1_y_xy_xxxxx_0, \
                             ta1_y_xy_xxxxy_0, \
                             ta1_y_xy_xxxxz_0, \
                             ta1_y_xy_xxxyy_0, \
                             ta1_y_xy_xxxyz_0, \
                             ta1_y_xy_xxxzz_0, \
                             ta1_y_xy_xxyyy_0, \
                             ta1_y_xy_xxyyz_0, \
                             ta1_y_xy_xxyzz_0, \
                             ta1_y_xy_xxzzz_0, \
                             ta1_y_xy_xyyyy_0, \
                             ta1_y_xy_xyyyz_0, \
                             ta1_y_xy_xyyzz_0, \
                             ta1_y_xy_xyzzz_0, \
                             ta1_y_xy_xzzzz_0, \
                             ta1_y_xy_yyyyy_0, \
                             ta1_y_xy_yyyyz_0, \
                             ta1_y_xy_yyyzz_0, \
                             ta1_y_xy_yyzzz_0, \
                             ta1_y_xy_yzzzz_0, \
                             ta1_y_xy_zzzzz_0, \
                             ta1_y_y_xxxxy_0,  \
                             ta1_y_y_xxxxy_1,  \
                             ta1_y_y_xxxy_0,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxxyy_0,  \
                             ta1_y_y_xxxyy_1,  \
                             ta1_y_y_xxxyz_0,  \
                             ta1_y_y_xxxyz_1,  \
                             ta1_y_y_xxyy_0,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyyy_0,  \
                             ta1_y_y_xxyyy_1,  \
                             ta1_y_y_xxyyz_0,  \
                             ta1_y_y_xxyyz_1,  \
                             ta1_y_y_xxyz_0,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xxyzz_0,  \
                             ta1_y_y_xxyzz_1,  \
                             ta1_y_y_xyyy_0,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyyy_0,  \
                             ta1_y_y_xyyyy_1,  \
                             ta1_y_y_xyyyz_0,  \
                             ta1_y_y_xyyyz_1,  \
                             ta1_y_y_xyyz_0,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyyzz_0,  \
                             ta1_y_y_xyyzz_1,  \
                             ta1_y_y_xyzz_0,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_xyzzz_0,  \
                             ta1_y_y_xyzzz_1,  \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyyy_0,  \
                             ta1_y_y_yyyyy_1,  \
                             ta1_y_y_yyyyz_0,  \
                             ta1_y_y_yyyyz_1,  \
                             ta1_y_y_yyyz_0,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyyzz_0,  \
                             ta1_y_y_yyyzz_1,  \
                             ta1_y_y_yyzz_0,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yyzzz_0,  \
                             ta1_y_y_yyzzz_1,  \
                             ta1_y_y_yzzz_0,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_yzzzz_0,  \
                             ta1_y_y_yzzzz_1,  \
                             ta1_y_y_zzzzz_0,  \
                             ta1_y_y_zzzzz_1,  \
                             ta_x_xxxxx_1,     \
                             ta_x_xxxxz_1,     \
                             ta_x_xxxzz_1,     \
                             ta_x_xxzzz_1,     \
                             ta_x_xzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xy_xxxxx_0[i] = ta_x_xxxxx_1[i] + ta1_y_x_xxxxx_0[i] * pa_y[i] - ta1_y_x_xxxxx_1[i] * pc_y[i];

        ta1_y_xy_xxxxy_0[i] =
            4.0 * ta1_y_y_xxxy_0[i] * fe_0 - 4.0 * ta1_y_y_xxxy_1[i] * fe_0 + ta1_y_y_xxxxy_0[i] * pa_x[i] - ta1_y_y_xxxxy_1[i] * pc_x[i];

        ta1_y_xy_xxxxz_0[i] = ta_x_xxxxz_1[i] + ta1_y_x_xxxxz_0[i] * pa_y[i] - ta1_y_x_xxxxz_1[i] * pc_y[i];

        ta1_y_xy_xxxyy_0[i] =
            3.0 * ta1_y_y_xxyy_0[i] * fe_0 - 3.0 * ta1_y_y_xxyy_1[i] * fe_0 + ta1_y_y_xxxyy_0[i] * pa_x[i] - ta1_y_y_xxxyy_1[i] * pc_x[i];

        ta1_y_xy_xxxyz_0[i] =
            3.0 * ta1_y_y_xxyz_0[i] * fe_0 - 3.0 * ta1_y_y_xxyz_1[i] * fe_0 + ta1_y_y_xxxyz_0[i] * pa_x[i] - ta1_y_y_xxxyz_1[i] * pc_x[i];

        ta1_y_xy_xxxzz_0[i] = ta_x_xxxzz_1[i] + ta1_y_x_xxxzz_0[i] * pa_y[i] - ta1_y_x_xxxzz_1[i] * pc_y[i];

        ta1_y_xy_xxyyy_0[i] =
            2.0 * ta1_y_y_xyyy_0[i] * fe_0 - 2.0 * ta1_y_y_xyyy_1[i] * fe_0 + ta1_y_y_xxyyy_0[i] * pa_x[i] - ta1_y_y_xxyyy_1[i] * pc_x[i];

        ta1_y_xy_xxyyz_0[i] =
            2.0 * ta1_y_y_xyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyz_1[i] * fe_0 + ta1_y_y_xxyyz_0[i] * pa_x[i] - ta1_y_y_xxyyz_1[i] * pc_x[i];

        ta1_y_xy_xxyzz_0[i] =
            2.0 * ta1_y_y_xyzz_0[i] * fe_0 - 2.0 * ta1_y_y_xyzz_1[i] * fe_0 + ta1_y_y_xxyzz_0[i] * pa_x[i] - ta1_y_y_xxyzz_1[i] * pc_x[i];

        ta1_y_xy_xxzzz_0[i] = ta_x_xxzzz_1[i] + ta1_y_x_xxzzz_0[i] * pa_y[i] - ta1_y_x_xxzzz_1[i] * pc_y[i];

        ta1_y_xy_xyyyy_0[i] = ta1_y_y_yyyy_0[i] * fe_0 - ta1_y_y_yyyy_1[i] * fe_0 + ta1_y_y_xyyyy_0[i] * pa_x[i] - ta1_y_y_xyyyy_1[i] * pc_x[i];

        ta1_y_xy_xyyyz_0[i] = ta1_y_y_yyyz_0[i] * fe_0 - ta1_y_y_yyyz_1[i] * fe_0 + ta1_y_y_xyyyz_0[i] * pa_x[i] - ta1_y_y_xyyyz_1[i] * pc_x[i];

        ta1_y_xy_xyyzz_0[i] = ta1_y_y_yyzz_0[i] * fe_0 - ta1_y_y_yyzz_1[i] * fe_0 + ta1_y_y_xyyzz_0[i] * pa_x[i] - ta1_y_y_xyyzz_1[i] * pc_x[i];

        ta1_y_xy_xyzzz_0[i] = ta1_y_y_yzzz_0[i] * fe_0 - ta1_y_y_yzzz_1[i] * fe_0 + ta1_y_y_xyzzz_0[i] * pa_x[i] - ta1_y_y_xyzzz_1[i] * pc_x[i];

        ta1_y_xy_xzzzz_0[i] = ta_x_xzzzz_1[i] + ta1_y_x_xzzzz_0[i] * pa_y[i] - ta1_y_x_xzzzz_1[i] * pc_y[i];

        ta1_y_xy_yyyyy_0[i] = ta1_y_y_yyyyy_0[i] * pa_x[i] - ta1_y_y_yyyyy_1[i] * pc_x[i];

        ta1_y_xy_yyyyz_0[i] = ta1_y_y_yyyyz_0[i] * pa_x[i] - ta1_y_y_yyyyz_1[i] * pc_x[i];

        ta1_y_xy_yyyzz_0[i] = ta1_y_y_yyyzz_0[i] * pa_x[i] - ta1_y_y_yyyzz_1[i] * pc_x[i];

        ta1_y_xy_yyzzz_0[i] = ta1_y_y_yyzzz_0[i] * pa_x[i] - ta1_y_y_yyzzz_1[i] * pc_x[i];

        ta1_y_xy_yzzzz_0[i] = ta1_y_y_yzzzz_0[i] * pa_x[i] - ta1_y_y_yzzzz_1[i] * pc_x[i];

        ta1_y_xy_zzzzz_0[i] = ta1_y_y_zzzzz_0[i] * pa_x[i] - ta1_y_y_zzzzz_1[i] * pc_x[i];
    }

    // Set up 168-189 components of targeted buffer : DH

    auto ta1_y_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 168);

    auto ta1_y_xz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 169);

    auto ta1_y_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 170);

    auto ta1_y_xz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 171);

    auto ta1_y_xz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 172);

    auto ta1_y_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 173);

    auto ta1_y_xz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 174);

    auto ta1_y_xz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 175);

    auto ta1_y_xz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 176);

    auto ta1_y_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 177);

    auto ta1_y_xz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 178);

    auto ta1_y_xz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 179);

    auto ta1_y_xz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 180);

    auto ta1_y_xz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 181);

    auto ta1_y_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 182);

    auto ta1_y_xz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 183);

    auto ta1_y_xz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 184);

    auto ta1_y_xz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 185);

    auto ta1_y_xz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 186);

    auto ta1_y_xz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 187);

    auto ta1_y_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 188);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_x_xxxxx_0,  \
                             ta1_y_x_xxxxx_1,  \
                             ta1_y_x_xxxxy_0,  \
                             ta1_y_x_xxxxy_1,  \
                             ta1_y_x_xxxyy_0,  \
                             ta1_y_x_xxxyy_1,  \
                             ta1_y_x_xxyyy_0,  \
                             ta1_y_x_xxyyy_1,  \
                             ta1_y_x_xyyyy_0,  \
                             ta1_y_x_xyyyy_1,  \
                             ta1_y_xz_xxxxx_0, \
                             ta1_y_xz_xxxxy_0, \
                             ta1_y_xz_xxxxz_0, \
                             ta1_y_xz_xxxyy_0, \
                             ta1_y_xz_xxxyz_0, \
                             ta1_y_xz_xxxzz_0, \
                             ta1_y_xz_xxyyy_0, \
                             ta1_y_xz_xxyyz_0, \
                             ta1_y_xz_xxyzz_0, \
                             ta1_y_xz_xxzzz_0, \
                             ta1_y_xz_xyyyy_0, \
                             ta1_y_xz_xyyyz_0, \
                             ta1_y_xz_xyyzz_0, \
                             ta1_y_xz_xyzzz_0, \
                             ta1_y_xz_xzzzz_0, \
                             ta1_y_xz_yyyyy_0, \
                             ta1_y_xz_yyyyz_0, \
                             ta1_y_xz_yyyzz_0, \
                             ta1_y_xz_yyzzz_0, \
                             ta1_y_xz_yzzzz_0, \
                             ta1_y_xz_zzzzz_0, \
                             ta1_y_z_xxxxz_0,  \
                             ta1_y_z_xxxxz_1,  \
                             ta1_y_z_xxxyz_0,  \
                             ta1_y_z_xxxyz_1,  \
                             ta1_y_z_xxxz_0,   \
                             ta1_y_z_xxxz_1,   \
                             ta1_y_z_xxxzz_0,  \
                             ta1_y_z_xxxzz_1,  \
                             ta1_y_z_xxyyz_0,  \
                             ta1_y_z_xxyyz_1,  \
                             ta1_y_z_xxyz_0,   \
                             ta1_y_z_xxyz_1,   \
                             ta1_y_z_xxyzz_0,  \
                             ta1_y_z_xxyzz_1,  \
                             ta1_y_z_xxzz_0,   \
                             ta1_y_z_xxzz_1,   \
                             ta1_y_z_xxzzz_0,  \
                             ta1_y_z_xxzzz_1,  \
                             ta1_y_z_xyyyz_0,  \
                             ta1_y_z_xyyyz_1,  \
                             ta1_y_z_xyyz_0,   \
                             ta1_y_z_xyyz_1,   \
                             ta1_y_z_xyyzz_0,  \
                             ta1_y_z_xyyzz_1,  \
                             ta1_y_z_xyzz_0,   \
                             ta1_y_z_xyzz_1,   \
                             ta1_y_z_xyzzz_0,  \
                             ta1_y_z_xyzzz_1,  \
                             ta1_y_z_xzzz_0,   \
                             ta1_y_z_xzzz_1,   \
                             ta1_y_z_xzzzz_0,  \
                             ta1_y_z_xzzzz_1,  \
                             ta1_y_z_yyyyy_0,  \
                             ta1_y_z_yyyyy_1,  \
                             ta1_y_z_yyyyz_0,  \
                             ta1_y_z_yyyyz_1,  \
                             ta1_y_z_yyyz_0,   \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyyzz_0,  \
                             ta1_y_z_yyyzz_1,  \
                             ta1_y_z_yyzz_0,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yyzzz_0,  \
                             ta1_y_z_yyzzz_1,  \
                             ta1_y_z_yzzz_0,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_yzzzz_0,  \
                             ta1_y_z_yzzzz_1,  \
                             ta1_y_z_zzzz_0,   \
                             ta1_y_z_zzzz_1,   \
                             ta1_y_z_zzzzz_0,  \
                             ta1_y_z_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xz_xxxxx_0[i] = ta1_y_x_xxxxx_0[i] * pa_z[i] - ta1_y_x_xxxxx_1[i] * pc_z[i];

        ta1_y_xz_xxxxy_0[i] = ta1_y_x_xxxxy_0[i] * pa_z[i] - ta1_y_x_xxxxy_1[i] * pc_z[i];

        ta1_y_xz_xxxxz_0[i] =
            4.0 * ta1_y_z_xxxz_0[i] * fe_0 - 4.0 * ta1_y_z_xxxz_1[i] * fe_0 + ta1_y_z_xxxxz_0[i] * pa_x[i] - ta1_y_z_xxxxz_1[i] * pc_x[i];

        ta1_y_xz_xxxyy_0[i] = ta1_y_x_xxxyy_0[i] * pa_z[i] - ta1_y_x_xxxyy_1[i] * pc_z[i];

        ta1_y_xz_xxxyz_0[i] =
            3.0 * ta1_y_z_xxyz_0[i] * fe_0 - 3.0 * ta1_y_z_xxyz_1[i] * fe_0 + ta1_y_z_xxxyz_0[i] * pa_x[i] - ta1_y_z_xxxyz_1[i] * pc_x[i];

        ta1_y_xz_xxxzz_0[i] =
            3.0 * ta1_y_z_xxzz_0[i] * fe_0 - 3.0 * ta1_y_z_xxzz_1[i] * fe_0 + ta1_y_z_xxxzz_0[i] * pa_x[i] - ta1_y_z_xxxzz_1[i] * pc_x[i];

        ta1_y_xz_xxyyy_0[i] = ta1_y_x_xxyyy_0[i] * pa_z[i] - ta1_y_x_xxyyy_1[i] * pc_z[i];

        ta1_y_xz_xxyyz_0[i] =
            2.0 * ta1_y_z_xyyz_0[i] * fe_0 - 2.0 * ta1_y_z_xyyz_1[i] * fe_0 + ta1_y_z_xxyyz_0[i] * pa_x[i] - ta1_y_z_xxyyz_1[i] * pc_x[i];

        ta1_y_xz_xxyzz_0[i] =
            2.0 * ta1_y_z_xyzz_0[i] * fe_0 - 2.0 * ta1_y_z_xyzz_1[i] * fe_0 + ta1_y_z_xxyzz_0[i] * pa_x[i] - ta1_y_z_xxyzz_1[i] * pc_x[i];

        ta1_y_xz_xxzzz_0[i] =
            2.0 * ta1_y_z_xzzz_0[i] * fe_0 - 2.0 * ta1_y_z_xzzz_1[i] * fe_0 + ta1_y_z_xxzzz_0[i] * pa_x[i] - ta1_y_z_xxzzz_1[i] * pc_x[i];

        ta1_y_xz_xyyyy_0[i] = ta1_y_x_xyyyy_0[i] * pa_z[i] - ta1_y_x_xyyyy_1[i] * pc_z[i];

        ta1_y_xz_xyyyz_0[i] = ta1_y_z_yyyz_0[i] * fe_0 - ta1_y_z_yyyz_1[i] * fe_0 + ta1_y_z_xyyyz_0[i] * pa_x[i] - ta1_y_z_xyyyz_1[i] * pc_x[i];

        ta1_y_xz_xyyzz_0[i] = ta1_y_z_yyzz_0[i] * fe_0 - ta1_y_z_yyzz_1[i] * fe_0 + ta1_y_z_xyyzz_0[i] * pa_x[i] - ta1_y_z_xyyzz_1[i] * pc_x[i];

        ta1_y_xz_xyzzz_0[i] = ta1_y_z_yzzz_0[i] * fe_0 - ta1_y_z_yzzz_1[i] * fe_0 + ta1_y_z_xyzzz_0[i] * pa_x[i] - ta1_y_z_xyzzz_1[i] * pc_x[i];

        ta1_y_xz_xzzzz_0[i] = ta1_y_z_zzzz_0[i] * fe_0 - ta1_y_z_zzzz_1[i] * fe_0 + ta1_y_z_xzzzz_0[i] * pa_x[i] - ta1_y_z_xzzzz_1[i] * pc_x[i];

        ta1_y_xz_yyyyy_0[i] = ta1_y_z_yyyyy_0[i] * pa_x[i] - ta1_y_z_yyyyy_1[i] * pc_x[i];

        ta1_y_xz_yyyyz_0[i] = ta1_y_z_yyyyz_0[i] * pa_x[i] - ta1_y_z_yyyyz_1[i] * pc_x[i];

        ta1_y_xz_yyyzz_0[i] = ta1_y_z_yyyzz_0[i] * pa_x[i] - ta1_y_z_yyyzz_1[i] * pc_x[i];

        ta1_y_xz_yyzzz_0[i] = ta1_y_z_yyzzz_0[i] * pa_x[i] - ta1_y_z_yyzzz_1[i] * pc_x[i];

        ta1_y_xz_yzzzz_0[i] = ta1_y_z_yzzzz_0[i] * pa_x[i] - ta1_y_z_yzzzz_1[i] * pc_x[i];

        ta1_y_xz_zzzzz_0[i] = ta1_y_z_zzzzz_0[i] * pa_x[i] - ta1_y_z_zzzzz_1[i] * pc_x[i];
    }

    // Set up 189-210 components of targeted buffer : DH

    auto ta1_y_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 189);

    auto ta1_y_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 190);

    auto ta1_y_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 191);

    auto ta1_y_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 192);

    auto ta1_y_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 193);

    auto ta1_y_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 194);

    auto ta1_y_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 195);

    auto ta1_y_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 196);

    auto ta1_y_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 197);

    auto ta1_y_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 198);

    auto ta1_y_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 199);

    auto ta1_y_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 200);

    auto ta1_y_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 201);

    auto ta1_y_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 202);

    auto ta1_y_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 203);

    auto ta1_y_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 204);

    auto ta1_y_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 205);

    auto ta1_y_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 206);

    auto ta1_y_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 207);

    auto ta1_y_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 208);

    auto ta1_y_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 209);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_y_xxxx_0,   \
                             ta1_y_y_xxxx_1,   \
                             ta1_y_y_xxxxx_0,  \
                             ta1_y_y_xxxxx_1,  \
                             ta1_y_y_xxxxy_0,  \
                             ta1_y_y_xxxxy_1,  \
                             ta1_y_y_xxxxz_0,  \
                             ta1_y_y_xxxxz_1,  \
                             ta1_y_y_xxxy_0,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxxyy_0,  \
                             ta1_y_y_xxxyy_1,  \
                             ta1_y_y_xxxyz_0,  \
                             ta1_y_y_xxxyz_1,  \
                             ta1_y_y_xxxz_0,   \
                             ta1_y_y_xxxz_1,   \
                             ta1_y_y_xxxzz_0,  \
                             ta1_y_y_xxxzz_1,  \
                             ta1_y_y_xxyy_0,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyyy_0,  \
                             ta1_y_y_xxyyy_1,  \
                             ta1_y_y_xxyyz_0,  \
                             ta1_y_y_xxyyz_1,  \
                             ta1_y_y_xxyz_0,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xxyzz_0,  \
                             ta1_y_y_xxyzz_1,  \
                             ta1_y_y_xxzz_0,   \
                             ta1_y_y_xxzz_1,   \
                             ta1_y_y_xxzzz_0,  \
                             ta1_y_y_xxzzz_1,  \
                             ta1_y_y_xyyy_0,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyyy_0,  \
                             ta1_y_y_xyyyy_1,  \
                             ta1_y_y_xyyyz_0,  \
                             ta1_y_y_xyyyz_1,  \
                             ta1_y_y_xyyz_0,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyyzz_0,  \
                             ta1_y_y_xyyzz_1,  \
                             ta1_y_y_xyzz_0,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_xyzzz_0,  \
                             ta1_y_y_xyzzz_1,  \
                             ta1_y_y_xzzz_0,   \
                             ta1_y_y_xzzz_1,   \
                             ta1_y_y_xzzzz_0,  \
                             ta1_y_y_xzzzz_1,  \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyyy_0,  \
                             ta1_y_y_yyyyy_1,  \
                             ta1_y_y_yyyyz_0,  \
                             ta1_y_y_yyyyz_1,  \
                             ta1_y_y_yyyz_0,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyyzz_0,  \
                             ta1_y_y_yyyzz_1,  \
                             ta1_y_y_yyzz_0,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yyzzz_0,  \
                             ta1_y_y_yyzzz_1,  \
                             ta1_y_y_yzzz_0,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_yzzzz_0,  \
                             ta1_y_y_yzzzz_1,  \
                             ta1_y_y_zzzz_0,   \
                             ta1_y_y_zzzz_1,   \
                             ta1_y_y_zzzzz_0,  \
                             ta1_y_y_zzzzz_1,  \
                             ta1_y_yy_xxxxx_0, \
                             ta1_y_yy_xxxxy_0, \
                             ta1_y_yy_xxxxz_0, \
                             ta1_y_yy_xxxyy_0, \
                             ta1_y_yy_xxxyz_0, \
                             ta1_y_yy_xxxzz_0, \
                             ta1_y_yy_xxyyy_0, \
                             ta1_y_yy_xxyyz_0, \
                             ta1_y_yy_xxyzz_0, \
                             ta1_y_yy_xxzzz_0, \
                             ta1_y_yy_xyyyy_0, \
                             ta1_y_yy_xyyyz_0, \
                             ta1_y_yy_xyyzz_0, \
                             ta1_y_yy_xyzzz_0, \
                             ta1_y_yy_xzzzz_0, \
                             ta1_y_yy_yyyyy_0, \
                             ta1_y_yy_yyyyz_0, \
                             ta1_y_yy_yyyzz_0, \
                             ta1_y_yy_yyzzz_0, \
                             ta1_y_yy_yzzzz_0, \
                             ta1_y_yy_zzzzz_0, \
                             ta_y_xxxxx_1,     \
                             ta_y_xxxxy_1,     \
                             ta_y_xxxxz_1,     \
                             ta_y_xxxyy_1,     \
                             ta_y_xxxyz_1,     \
                             ta_y_xxxzz_1,     \
                             ta_y_xxyyy_1,     \
                             ta_y_xxyyz_1,     \
                             ta_y_xxyzz_1,     \
                             ta_y_xxzzz_1,     \
                             ta_y_xyyyy_1,     \
                             ta_y_xyyyz_1,     \
                             ta_y_xyyzz_1,     \
                             ta_y_xyzzz_1,     \
                             ta_y_xzzzz_1,     \
                             ta_y_yyyyy_1,     \
                             ta_y_yyyyz_1,     \
                             ta_y_yyyzz_1,     \
                             ta_y_yyzzz_1,     \
                             ta_y_yzzzz_1,     \
                             ta_y_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yy_xxxxx_0[i] =
            ta1_y_0_xxxxx_0[i] * fe_0 - ta1_y_0_xxxxx_1[i] * fe_0 + ta_y_xxxxx_1[i] + ta1_y_y_xxxxx_0[i] * pa_y[i] - ta1_y_y_xxxxx_1[i] * pc_y[i];

        ta1_y_yy_xxxxy_0[i] = ta1_y_0_xxxxy_0[i] * fe_0 - ta1_y_0_xxxxy_1[i] * fe_0 + ta1_y_y_xxxx_0[i] * fe_0 - ta1_y_y_xxxx_1[i] * fe_0 +
                              ta_y_xxxxy_1[i] + ta1_y_y_xxxxy_0[i] * pa_y[i] - ta1_y_y_xxxxy_1[i] * pc_y[i];

        ta1_y_yy_xxxxz_0[i] =
            ta1_y_0_xxxxz_0[i] * fe_0 - ta1_y_0_xxxxz_1[i] * fe_0 + ta_y_xxxxz_1[i] + ta1_y_y_xxxxz_0[i] * pa_y[i] - ta1_y_y_xxxxz_1[i] * pc_y[i];

        ta1_y_yy_xxxyy_0[i] = ta1_y_0_xxxyy_0[i] * fe_0 - ta1_y_0_xxxyy_1[i] * fe_0 + 2.0 * ta1_y_y_xxxy_0[i] * fe_0 -
                              2.0 * ta1_y_y_xxxy_1[i] * fe_0 + ta_y_xxxyy_1[i] + ta1_y_y_xxxyy_0[i] * pa_y[i] - ta1_y_y_xxxyy_1[i] * pc_y[i];

        ta1_y_yy_xxxyz_0[i] = ta1_y_0_xxxyz_0[i] * fe_0 - ta1_y_0_xxxyz_1[i] * fe_0 + ta1_y_y_xxxz_0[i] * fe_0 - ta1_y_y_xxxz_1[i] * fe_0 +
                              ta_y_xxxyz_1[i] + ta1_y_y_xxxyz_0[i] * pa_y[i] - ta1_y_y_xxxyz_1[i] * pc_y[i];

        ta1_y_yy_xxxzz_0[i] =
            ta1_y_0_xxxzz_0[i] * fe_0 - ta1_y_0_xxxzz_1[i] * fe_0 + ta_y_xxxzz_1[i] + ta1_y_y_xxxzz_0[i] * pa_y[i] - ta1_y_y_xxxzz_1[i] * pc_y[i];

        ta1_y_yy_xxyyy_0[i] = ta1_y_0_xxyyy_0[i] * fe_0 - ta1_y_0_xxyyy_1[i] * fe_0 + 3.0 * ta1_y_y_xxyy_0[i] * fe_0 -
                              3.0 * ta1_y_y_xxyy_1[i] * fe_0 + ta_y_xxyyy_1[i] + ta1_y_y_xxyyy_0[i] * pa_y[i] - ta1_y_y_xxyyy_1[i] * pc_y[i];

        ta1_y_yy_xxyyz_0[i] = ta1_y_0_xxyyz_0[i] * fe_0 - ta1_y_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_y_xxyz_0[i] * fe_0 -
                              2.0 * ta1_y_y_xxyz_1[i] * fe_0 + ta_y_xxyyz_1[i] + ta1_y_y_xxyyz_0[i] * pa_y[i] - ta1_y_y_xxyyz_1[i] * pc_y[i];

        ta1_y_yy_xxyzz_0[i] = ta1_y_0_xxyzz_0[i] * fe_0 - ta1_y_0_xxyzz_1[i] * fe_0 + ta1_y_y_xxzz_0[i] * fe_0 - ta1_y_y_xxzz_1[i] * fe_0 +
                              ta_y_xxyzz_1[i] + ta1_y_y_xxyzz_0[i] * pa_y[i] - ta1_y_y_xxyzz_1[i] * pc_y[i];

        ta1_y_yy_xxzzz_0[i] =
            ta1_y_0_xxzzz_0[i] * fe_0 - ta1_y_0_xxzzz_1[i] * fe_0 + ta_y_xxzzz_1[i] + ta1_y_y_xxzzz_0[i] * pa_y[i] - ta1_y_y_xxzzz_1[i] * pc_y[i];

        ta1_y_yy_xyyyy_0[i] = ta1_y_0_xyyyy_0[i] * fe_0 - ta1_y_0_xyyyy_1[i] * fe_0 + 4.0 * ta1_y_y_xyyy_0[i] * fe_0 -
                              4.0 * ta1_y_y_xyyy_1[i] * fe_0 + ta_y_xyyyy_1[i] + ta1_y_y_xyyyy_0[i] * pa_y[i] - ta1_y_y_xyyyy_1[i] * pc_y[i];

        ta1_y_yy_xyyyz_0[i] = ta1_y_0_xyyyz_0[i] * fe_0 - ta1_y_0_xyyyz_1[i] * fe_0 + 3.0 * ta1_y_y_xyyz_0[i] * fe_0 -
                              3.0 * ta1_y_y_xyyz_1[i] * fe_0 + ta_y_xyyyz_1[i] + ta1_y_y_xyyyz_0[i] * pa_y[i] - ta1_y_y_xyyyz_1[i] * pc_y[i];

        ta1_y_yy_xyyzz_0[i] = ta1_y_0_xyyzz_0[i] * fe_0 - ta1_y_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_y_xyzz_0[i] * fe_0 -
                              2.0 * ta1_y_y_xyzz_1[i] * fe_0 + ta_y_xyyzz_1[i] + ta1_y_y_xyyzz_0[i] * pa_y[i] - ta1_y_y_xyyzz_1[i] * pc_y[i];

        ta1_y_yy_xyzzz_0[i] = ta1_y_0_xyzzz_0[i] * fe_0 - ta1_y_0_xyzzz_1[i] * fe_0 + ta1_y_y_xzzz_0[i] * fe_0 - ta1_y_y_xzzz_1[i] * fe_0 +
                              ta_y_xyzzz_1[i] + ta1_y_y_xyzzz_0[i] * pa_y[i] - ta1_y_y_xyzzz_1[i] * pc_y[i];

        ta1_y_yy_xzzzz_0[i] =
            ta1_y_0_xzzzz_0[i] * fe_0 - ta1_y_0_xzzzz_1[i] * fe_0 + ta_y_xzzzz_1[i] + ta1_y_y_xzzzz_0[i] * pa_y[i] - ta1_y_y_xzzzz_1[i] * pc_y[i];

        ta1_y_yy_yyyyy_0[i] = ta1_y_0_yyyyy_0[i] * fe_0 - ta1_y_0_yyyyy_1[i] * fe_0 + 5.0 * ta1_y_y_yyyy_0[i] * fe_0 -
                              5.0 * ta1_y_y_yyyy_1[i] * fe_0 + ta_y_yyyyy_1[i] + ta1_y_y_yyyyy_0[i] * pa_y[i] - ta1_y_y_yyyyy_1[i] * pc_y[i];

        ta1_y_yy_yyyyz_0[i] = ta1_y_0_yyyyz_0[i] * fe_0 - ta1_y_0_yyyyz_1[i] * fe_0 + 4.0 * ta1_y_y_yyyz_0[i] * fe_0 -
                              4.0 * ta1_y_y_yyyz_1[i] * fe_0 + ta_y_yyyyz_1[i] + ta1_y_y_yyyyz_0[i] * pa_y[i] - ta1_y_y_yyyyz_1[i] * pc_y[i];

        ta1_y_yy_yyyzz_0[i] = ta1_y_0_yyyzz_0[i] * fe_0 - ta1_y_0_yyyzz_1[i] * fe_0 + 3.0 * ta1_y_y_yyzz_0[i] * fe_0 -
                              3.0 * ta1_y_y_yyzz_1[i] * fe_0 + ta_y_yyyzz_1[i] + ta1_y_y_yyyzz_0[i] * pa_y[i] - ta1_y_y_yyyzz_1[i] * pc_y[i];

        ta1_y_yy_yyzzz_0[i] = ta1_y_0_yyzzz_0[i] * fe_0 - ta1_y_0_yyzzz_1[i] * fe_0 + 2.0 * ta1_y_y_yzzz_0[i] * fe_0 -
                              2.0 * ta1_y_y_yzzz_1[i] * fe_0 + ta_y_yyzzz_1[i] + ta1_y_y_yyzzz_0[i] * pa_y[i] - ta1_y_y_yyzzz_1[i] * pc_y[i];

        ta1_y_yy_yzzzz_0[i] = ta1_y_0_yzzzz_0[i] * fe_0 - ta1_y_0_yzzzz_1[i] * fe_0 + ta1_y_y_zzzz_0[i] * fe_0 - ta1_y_y_zzzz_1[i] * fe_0 +
                              ta_y_yzzzz_1[i] + ta1_y_y_yzzzz_0[i] * pa_y[i] - ta1_y_y_yzzzz_1[i] * pc_y[i];

        ta1_y_yy_zzzzz_0[i] =
            ta1_y_0_zzzzz_0[i] * fe_0 - ta1_y_0_zzzzz_1[i] * fe_0 + ta_y_zzzzz_1[i] + ta1_y_y_zzzzz_0[i] * pa_y[i] - ta1_y_y_zzzzz_1[i] * pc_y[i];
    }

    // Set up 210-231 components of targeted buffer : DH

    auto ta1_y_yz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 210);

    auto ta1_y_yz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 211);

    auto ta1_y_yz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 212);

    auto ta1_y_yz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 213);

    auto ta1_y_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 214);

    auto ta1_y_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 215);

    auto ta1_y_yz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 216);

    auto ta1_y_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 217);

    auto ta1_y_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 218);

    auto ta1_y_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 219);

    auto ta1_y_yz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 220);

    auto ta1_y_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 221);

    auto ta1_y_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 222);

    auto ta1_y_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 223);

    auto ta1_y_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 224);

    auto ta1_y_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 225);

    auto ta1_y_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 226);

    auto ta1_y_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 227);

    auto ta1_y_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 228);

    auto ta1_y_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 229);

    auto ta1_y_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 230);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_y_xxxxx_0,  \
                             ta1_y_y_xxxxx_1,  \
                             ta1_y_y_xxxxy_0,  \
                             ta1_y_y_xxxxy_1,  \
                             ta1_y_y_xxxy_0,   \
                             ta1_y_y_xxxy_1,   \
                             ta1_y_y_xxxyy_0,  \
                             ta1_y_y_xxxyy_1,  \
                             ta1_y_y_xxxyz_0,  \
                             ta1_y_y_xxxyz_1,  \
                             ta1_y_y_xxyy_0,   \
                             ta1_y_y_xxyy_1,   \
                             ta1_y_y_xxyyy_0,  \
                             ta1_y_y_xxyyy_1,  \
                             ta1_y_y_xxyyz_0,  \
                             ta1_y_y_xxyyz_1,  \
                             ta1_y_y_xxyz_0,   \
                             ta1_y_y_xxyz_1,   \
                             ta1_y_y_xxyzz_0,  \
                             ta1_y_y_xxyzz_1,  \
                             ta1_y_y_xyyy_0,   \
                             ta1_y_y_xyyy_1,   \
                             ta1_y_y_xyyyy_0,  \
                             ta1_y_y_xyyyy_1,  \
                             ta1_y_y_xyyyz_0,  \
                             ta1_y_y_xyyyz_1,  \
                             ta1_y_y_xyyz_0,   \
                             ta1_y_y_xyyz_1,   \
                             ta1_y_y_xyyzz_0,  \
                             ta1_y_y_xyyzz_1,  \
                             ta1_y_y_xyzz_0,   \
                             ta1_y_y_xyzz_1,   \
                             ta1_y_y_xyzzz_0,  \
                             ta1_y_y_xyzzz_1,  \
                             ta1_y_y_yyyy_0,   \
                             ta1_y_y_yyyy_1,   \
                             ta1_y_y_yyyyy_0,  \
                             ta1_y_y_yyyyy_1,  \
                             ta1_y_y_yyyyz_0,  \
                             ta1_y_y_yyyyz_1,  \
                             ta1_y_y_yyyz_0,   \
                             ta1_y_y_yyyz_1,   \
                             ta1_y_y_yyyzz_0,  \
                             ta1_y_y_yyyzz_1,  \
                             ta1_y_y_yyzz_0,   \
                             ta1_y_y_yyzz_1,   \
                             ta1_y_y_yyzzz_0,  \
                             ta1_y_y_yyzzz_1,  \
                             ta1_y_y_yzzz_0,   \
                             ta1_y_y_yzzz_1,   \
                             ta1_y_y_yzzzz_0,  \
                             ta1_y_y_yzzzz_1,  \
                             ta1_y_yz_xxxxx_0, \
                             ta1_y_yz_xxxxy_0, \
                             ta1_y_yz_xxxxz_0, \
                             ta1_y_yz_xxxyy_0, \
                             ta1_y_yz_xxxyz_0, \
                             ta1_y_yz_xxxzz_0, \
                             ta1_y_yz_xxyyy_0, \
                             ta1_y_yz_xxyyz_0, \
                             ta1_y_yz_xxyzz_0, \
                             ta1_y_yz_xxzzz_0, \
                             ta1_y_yz_xyyyy_0, \
                             ta1_y_yz_xyyyz_0, \
                             ta1_y_yz_xyyzz_0, \
                             ta1_y_yz_xyzzz_0, \
                             ta1_y_yz_xzzzz_0, \
                             ta1_y_yz_yyyyy_0, \
                             ta1_y_yz_yyyyz_0, \
                             ta1_y_yz_yyyzz_0, \
                             ta1_y_yz_yyzzz_0, \
                             ta1_y_yz_yzzzz_0, \
                             ta1_y_yz_zzzzz_0, \
                             ta1_y_z_xxxxz_0,  \
                             ta1_y_z_xxxxz_1,  \
                             ta1_y_z_xxxzz_0,  \
                             ta1_y_z_xxxzz_1,  \
                             ta1_y_z_xxzzz_0,  \
                             ta1_y_z_xxzzz_1,  \
                             ta1_y_z_xzzzz_0,  \
                             ta1_y_z_xzzzz_1,  \
                             ta1_y_z_zzzzz_0,  \
                             ta1_y_z_zzzzz_1,  \
                             ta_z_xxxxz_1,     \
                             ta_z_xxxzz_1,     \
                             ta_z_xxzzz_1,     \
                             ta_z_xzzzz_1,     \
                             ta_z_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yz_xxxxx_0[i] = ta1_y_y_xxxxx_0[i] * pa_z[i] - ta1_y_y_xxxxx_1[i] * pc_z[i];

        ta1_y_yz_xxxxy_0[i] = ta1_y_y_xxxxy_0[i] * pa_z[i] - ta1_y_y_xxxxy_1[i] * pc_z[i];

        ta1_y_yz_xxxxz_0[i] = ta_z_xxxxz_1[i] + ta1_y_z_xxxxz_0[i] * pa_y[i] - ta1_y_z_xxxxz_1[i] * pc_y[i];

        ta1_y_yz_xxxyy_0[i] = ta1_y_y_xxxyy_0[i] * pa_z[i] - ta1_y_y_xxxyy_1[i] * pc_z[i];

        ta1_y_yz_xxxyz_0[i] = ta1_y_y_xxxy_0[i] * fe_0 - ta1_y_y_xxxy_1[i] * fe_0 + ta1_y_y_xxxyz_0[i] * pa_z[i] - ta1_y_y_xxxyz_1[i] * pc_z[i];

        ta1_y_yz_xxxzz_0[i] = ta_z_xxxzz_1[i] + ta1_y_z_xxxzz_0[i] * pa_y[i] - ta1_y_z_xxxzz_1[i] * pc_y[i];

        ta1_y_yz_xxyyy_0[i] = ta1_y_y_xxyyy_0[i] * pa_z[i] - ta1_y_y_xxyyy_1[i] * pc_z[i];

        ta1_y_yz_xxyyz_0[i] = ta1_y_y_xxyy_0[i] * fe_0 - ta1_y_y_xxyy_1[i] * fe_0 + ta1_y_y_xxyyz_0[i] * pa_z[i] - ta1_y_y_xxyyz_1[i] * pc_z[i];

        ta1_y_yz_xxyzz_0[i] =
            2.0 * ta1_y_y_xxyz_0[i] * fe_0 - 2.0 * ta1_y_y_xxyz_1[i] * fe_0 + ta1_y_y_xxyzz_0[i] * pa_z[i] - ta1_y_y_xxyzz_1[i] * pc_z[i];

        ta1_y_yz_xxzzz_0[i] = ta_z_xxzzz_1[i] + ta1_y_z_xxzzz_0[i] * pa_y[i] - ta1_y_z_xxzzz_1[i] * pc_y[i];

        ta1_y_yz_xyyyy_0[i] = ta1_y_y_xyyyy_0[i] * pa_z[i] - ta1_y_y_xyyyy_1[i] * pc_z[i];

        ta1_y_yz_xyyyz_0[i] = ta1_y_y_xyyy_0[i] * fe_0 - ta1_y_y_xyyy_1[i] * fe_0 + ta1_y_y_xyyyz_0[i] * pa_z[i] - ta1_y_y_xyyyz_1[i] * pc_z[i];

        ta1_y_yz_xyyzz_0[i] =
            2.0 * ta1_y_y_xyyz_0[i] * fe_0 - 2.0 * ta1_y_y_xyyz_1[i] * fe_0 + ta1_y_y_xyyzz_0[i] * pa_z[i] - ta1_y_y_xyyzz_1[i] * pc_z[i];

        ta1_y_yz_xyzzz_0[i] =
            3.0 * ta1_y_y_xyzz_0[i] * fe_0 - 3.0 * ta1_y_y_xyzz_1[i] * fe_0 + ta1_y_y_xyzzz_0[i] * pa_z[i] - ta1_y_y_xyzzz_1[i] * pc_z[i];

        ta1_y_yz_xzzzz_0[i] = ta_z_xzzzz_1[i] + ta1_y_z_xzzzz_0[i] * pa_y[i] - ta1_y_z_xzzzz_1[i] * pc_y[i];

        ta1_y_yz_yyyyy_0[i] = ta1_y_y_yyyyy_0[i] * pa_z[i] - ta1_y_y_yyyyy_1[i] * pc_z[i];

        ta1_y_yz_yyyyz_0[i] = ta1_y_y_yyyy_0[i] * fe_0 - ta1_y_y_yyyy_1[i] * fe_0 + ta1_y_y_yyyyz_0[i] * pa_z[i] - ta1_y_y_yyyyz_1[i] * pc_z[i];

        ta1_y_yz_yyyzz_0[i] =
            2.0 * ta1_y_y_yyyz_0[i] * fe_0 - 2.0 * ta1_y_y_yyyz_1[i] * fe_0 + ta1_y_y_yyyzz_0[i] * pa_z[i] - ta1_y_y_yyyzz_1[i] * pc_z[i];

        ta1_y_yz_yyzzz_0[i] =
            3.0 * ta1_y_y_yyzz_0[i] * fe_0 - 3.0 * ta1_y_y_yyzz_1[i] * fe_0 + ta1_y_y_yyzzz_0[i] * pa_z[i] - ta1_y_y_yyzzz_1[i] * pc_z[i];

        ta1_y_yz_yzzzz_0[i] =
            4.0 * ta1_y_y_yzzz_0[i] * fe_0 - 4.0 * ta1_y_y_yzzz_1[i] * fe_0 + ta1_y_y_yzzzz_0[i] * pa_z[i] - ta1_y_y_yzzzz_1[i] * pc_z[i];

        ta1_y_yz_zzzzz_0[i] = ta_z_zzzzz_1[i] + ta1_y_z_zzzzz_0[i] * pa_y[i] - ta1_y_z_zzzzz_1[i] * pc_y[i];
    }

    // Set up 231-252 components of targeted buffer : DH

    auto ta1_y_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 231);

    auto ta1_y_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 232);

    auto ta1_y_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 233);

    auto ta1_y_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 234);

    auto ta1_y_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 235);

    auto ta1_y_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 236);

    auto ta1_y_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 237);

    auto ta1_y_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 238);

    auto ta1_y_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 239);

    auto ta1_y_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 240);

    auto ta1_y_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 241);

    auto ta1_y_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 242);

    auto ta1_y_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 243);

    auto ta1_y_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 244);

    auto ta1_y_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 245);

    auto ta1_y_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 246);

    auto ta1_y_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 247);

    auto ta1_y_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 248);

    auto ta1_y_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 249);

    auto ta1_y_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 250);

    auto ta1_y_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 251);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxz_0,  \
                             ta1_y_0_xxxxz_1,  \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyz_0,  \
                             ta1_y_0_xxxyz_1,  \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyz_0,  \
                             ta1_y_0_xxyyz_1,  \
                             ta1_y_0_xxyzz_0,  \
                             ta1_y_0_xxyzz_1,  \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyz_0,  \
                             ta1_y_0_xyyyz_1,  \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyzzz_0,  \
                             ta1_y_0_xyzzz_1,  \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_z_xxxx_0,   \
                             ta1_y_z_xxxx_1,   \
                             ta1_y_z_xxxxx_0,  \
                             ta1_y_z_xxxxx_1,  \
                             ta1_y_z_xxxxy_0,  \
                             ta1_y_z_xxxxy_1,  \
                             ta1_y_z_xxxxz_0,  \
                             ta1_y_z_xxxxz_1,  \
                             ta1_y_z_xxxy_0,   \
                             ta1_y_z_xxxy_1,   \
                             ta1_y_z_xxxyy_0,  \
                             ta1_y_z_xxxyy_1,  \
                             ta1_y_z_xxxyz_0,  \
                             ta1_y_z_xxxyz_1,  \
                             ta1_y_z_xxxz_0,   \
                             ta1_y_z_xxxz_1,   \
                             ta1_y_z_xxxzz_0,  \
                             ta1_y_z_xxxzz_1,  \
                             ta1_y_z_xxyy_0,   \
                             ta1_y_z_xxyy_1,   \
                             ta1_y_z_xxyyy_0,  \
                             ta1_y_z_xxyyy_1,  \
                             ta1_y_z_xxyyz_0,  \
                             ta1_y_z_xxyyz_1,  \
                             ta1_y_z_xxyz_0,   \
                             ta1_y_z_xxyz_1,   \
                             ta1_y_z_xxyzz_0,  \
                             ta1_y_z_xxyzz_1,  \
                             ta1_y_z_xxzz_0,   \
                             ta1_y_z_xxzz_1,   \
                             ta1_y_z_xxzzz_0,  \
                             ta1_y_z_xxzzz_1,  \
                             ta1_y_z_xyyy_0,   \
                             ta1_y_z_xyyy_1,   \
                             ta1_y_z_xyyyy_0,  \
                             ta1_y_z_xyyyy_1,  \
                             ta1_y_z_xyyyz_0,  \
                             ta1_y_z_xyyyz_1,  \
                             ta1_y_z_xyyz_0,   \
                             ta1_y_z_xyyz_1,   \
                             ta1_y_z_xyyzz_0,  \
                             ta1_y_z_xyyzz_1,  \
                             ta1_y_z_xyzz_0,   \
                             ta1_y_z_xyzz_1,   \
                             ta1_y_z_xyzzz_0,  \
                             ta1_y_z_xyzzz_1,  \
                             ta1_y_z_xzzz_0,   \
                             ta1_y_z_xzzz_1,   \
                             ta1_y_z_xzzzz_0,  \
                             ta1_y_z_xzzzz_1,  \
                             ta1_y_z_yyyy_0,   \
                             ta1_y_z_yyyy_1,   \
                             ta1_y_z_yyyyy_0,  \
                             ta1_y_z_yyyyy_1,  \
                             ta1_y_z_yyyyz_0,  \
                             ta1_y_z_yyyyz_1,  \
                             ta1_y_z_yyyz_0,   \
                             ta1_y_z_yyyz_1,   \
                             ta1_y_z_yyyzz_0,  \
                             ta1_y_z_yyyzz_1,  \
                             ta1_y_z_yyzz_0,   \
                             ta1_y_z_yyzz_1,   \
                             ta1_y_z_yyzzz_0,  \
                             ta1_y_z_yyzzz_1,  \
                             ta1_y_z_yzzz_0,   \
                             ta1_y_z_yzzz_1,   \
                             ta1_y_z_yzzzz_0,  \
                             ta1_y_z_yzzzz_1,  \
                             ta1_y_z_zzzz_0,   \
                             ta1_y_z_zzzz_1,   \
                             ta1_y_z_zzzzz_0,  \
                             ta1_y_z_zzzzz_1,  \
                             ta1_y_zz_xxxxx_0, \
                             ta1_y_zz_xxxxy_0, \
                             ta1_y_zz_xxxxz_0, \
                             ta1_y_zz_xxxyy_0, \
                             ta1_y_zz_xxxyz_0, \
                             ta1_y_zz_xxxzz_0, \
                             ta1_y_zz_xxyyy_0, \
                             ta1_y_zz_xxyyz_0, \
                             ta1_y_zz_xxyzz_0, \
                             ta1_y_zz_xxzzz_0, \
                             ta1_y_zz_xyyyy_0, \
                             ta1_y_zz_xyyyz_0, \
                             ta1_y_zz_xyyzz_0, \
                             ta1_y_zz_xyzzz_0, \
                             ta1_y_zz_xzzzz_0, \
                             ta1_y_zz_yyyyy_0, \
                             ta1_y_zz_yyyyz_0, \
                             ta1_y_zz_yyyzz_0, \
                             ta1_y_zz_yyzzz_0, \
                             ta1_y_zz_yzzzz_0, \
                             ta1_y_zz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zz_xxxxx_0[i] = ta1_y_0_xxxxx_0[i] * fe_0 - ta1_y_0_xxxxx_1[i] * fe_0 + ta1_y_z_xxxxx_0[i] * pa_z[i] - ta1_y_z_xxxxx_1[i] * pc_z[i];

        ta1_y_zz_xxxxy_0[i] = ta1_y_0_xxxxy_0[i] * fe_0 - ta1_y_0_xxxxy_1[i] * fe_0 + ta1_y_z_xxxxy_0[i] * pa_z[i] - ta1_y_z_xxxxy_1[i] * pc_z[i];

        ta1_y_zz_xxxxz_0[i] = ta1_y_0_xxxxz_0[i] * fe_0 - ta1_y_0_xxxxz_1[i] * fe_0 + ta1_y_z_xxxx_0[i] * fe_0 - ta1_y_z_xxxx_1[i] * fe_0 +
                              ta1_y_z_xxxxz_0[i] * pa_z[i] - ta1_y_z_xxxxz_1[i] * pc_z[i];

        ta1_y_zz_xxxyy_0[i] = ta1_y_0_xxxyy_0[i] * fe_0 - ta1_y_0_xxxyy_1[i] * fe_0 + ta1_y_z_xxxyy_0[i] * pa_z[i] - ta1_y_z_xxxyy_1[i] * pc_z[i];

        ta1_y_zz_xxxyz_0[i] = ta1_y_0_xxxyz_0[i] * fe_0 - ta1_y_0_xxxyz_1[i] * fe_0 + ta1_y_z_xxxy_0[i] * fe_0 - ta1_y_z_xxxy_1[i] * fe_0 +
                              ta1_y_z_xxxyz_0[i] * pa_z[i] - ta1_y_z_xxxyz_1[i] * pc_z[i];

        ta1_y_zz_xxxzz_0[i] = ta1_y_0_xxxzz_0[i] * fe_0 - ta1_y_0_xxxzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxxz_0[i] * fe_0 -
                              2.0 * ta1_y_z_xxxz_1[i] * fe_0 + ta1_y_z_xxxzz_0[i] * pa_z[i] - ta1_y_z_xxxzz_1[i] * pc_z[i];

        ta1_y_zz_xxyyy_0[i] = ta1_y_0_xxyyy_0[i] * fe_0 - ta1_y_0_xxyyy_1[i] * fe_0 + ta1_y_z_xxyyy_0[i] * pa_z[i] - ta1_y_z_xxyyy_1[i] * pc_z[i];

        ta1_y_zz_xxyyz_0[i] = ta1_y_0_xxyyz_0[i] * fe_0 - ta1_y_0_xxyyz_1[i] * fe_0 + ta1_y_z_xxyy_0[i] * fe_0 - ta1_y_z_xxyy_1[i] * fe_0 +
                              ta1_y_z_xxyyz_0[i] * pa_z[i] - ta1_y_z_xxyyz_1[i] * pc_z[i];

        ta1_y_zz_xxyzz_0[i] = ta1_y_0_xxyzz_0[i] * fe_0 - ta1_y_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xxyz_0[i] * fe_0 -
                              2.0 * ta1_y_z_xxyz_1[i] * fe_0 + ta1_y_z_xxyzz_0[i] * pa_z[i] - ta1_y_z_xxyzz_1[i] * pc_z[i];

        ta1_y_zz_xxzzz_0[i] = ta1_y_0_xxzzz_0[i] * fe_0 - ta1_y_0_xxzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xxzz_0[i] * fe_0 -
                              3.0 * ta1_y_z_xxzz_1[i] * fe_0 + ta1_y_z_xxzzz_0[i] * pa_z[i] - ta1_y_z_xxzzz_1[i] * pc_z[i];

        ta1_y_zz_xyyyy_0[i] = ta1_y_0_xyyyy_0[i] * fe_0 - ta1_y_0_xyyyy_1[i] * fe_0 + ta1_y_z_xyyyy_0[i] * pa_z[i] - ta1_y_z_xyyyy_1[i] * pc_z[i];

        ta1_y_zz_xyyyz_0[i] = ta1_y_0_xyyyz_0[i] * fe_0 - ta1_y_0_xyyyz_1[i] * fe_0 + ta1_y_z_xyyy_0[i] * fe_0 - ta1_y_z_xyyy_1[i] * fe_0 +
                              ta1_y_z_xyyyz_0[i] * pa_z[i] - ta1_y_z_xyyyz_1[i] * pc_z[i];

        ta1_y_zz_xyyzz_0[i] = ta1_y_0_xyyzz_0[i] * fe_0 - ta1_y_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_z_xyyz_0[i] * fe_0 -
                              2.0 * ta1_y_z_xyyz_1[i] * fe_0 + ta1_y_z_xyyzz_0[i] * pa_z[i] - ta1_y_z_xyyzz_1[i] * pc_z[i];

        ta1_y_zz_xyzzz_0[i] = ta1_y_0_xyzzz_0[i] * fe_0 - ta1_y_0_xyzzz_1[i] * fe_0 + 3.0 * ta1_y_z_xyzz_0[i] * fe_0 -
                              3.0 * ta1_y_z_xyzz_1[i] * fe_0 + ta1_y_z_xyzzz_0[i] * pa_z[i] - ta1_y_z_xyzzz_1[i] * pc_z[i];

        ta1_y_zz_xzzzz_0[i] = ta1_y_0_xzzzz_0[i] * fe_0 - ta1_y_0_xzzzz_1[i] * fe_0 + 4.0 * ta1_y_z_xzzz_0[i] * fe_0 -
                              4.0 * ta1_y_z_xzzz_1[i] * fe_0 + ta1_y_z_xzzzz_0[i] * pa_z[i] - ta1_y_z_xzzzz_1[i] * pc_z[i];

        ta1_y_zz_yyyyy_0[i] = ta1_y_0_yyyyy_0[i] * fe_0 - ta1_y_0_yyyyy_1[i] * fe_0 + ta1_y_z_yyyyy_0[i] * pa_z[i] - ta1_y_z_yyyyy_1[i] * pc_z[i];

        ta1_y_zz_yyyyz_0[i] = ta1_y_0_yyyyz_0[i] * fe_0 - ta1_y_0_yyyyz_1[i] * fe_0 + ta1_y_z_yyyy_0[i] * fe_0 - ta1_y_z_yyyy_1[i] * fe_0 +
                              ta1_y_z_yyyyz_0[i] * pa_z[i] - ta1_y_z_yyyyz_1[i] * pc_z[i];

        ta1_y_zz_yyyzz_0[i] = ta1_y_0_yyyzz_0[i] * fe_0 - ta1_y_0_yyyzz_1[i] * fe_0 + 2.0 * ta1_y_z_yyyz_0[i] * fe_0 -
                              2.0 * ta1_y_z_yyyz_1[i] * fe_0 + ta1_y_z_yyyzz_0[i] * pa_z[i] - ta1_y_z_yyyzz_1[i] * pc_z[i];

        ta1_y_zz_yyzzz_0[i] = ta1_y_0_yyzzz_0[i] * fe_0 - ta1_y_0_yyzzz_1[i] * fe_0 + 3.0 * ta1_y_z_yyzz_0[i] * fe_0 -
                              3.0 * ta1_y_z_yyzz_1[i] * fe_0 + ta1_y_z_yyzzz_0[i] * pa_z[i] - ta1_y_z_yyzzz_1[i] * pc_z[i];

        ta1_y_zz_yzzzz_0[i] = ta1_y_0_yzzzz_0[i] * fe_0 - ta1_y_0_yzzzz_1[i] * fe_0 + 4.0 * ta1_y_z_yzzz_0[i] * fe_0 -
                              4.0 * ta1_y_z_yzzz_1[i] * fe_0 + ta1_y_z_yzzzz_0[i] * pa_z[i] - ta1_y_z_yzzzz_1[i] * pc_z[i];

        ta1_y_zz_zzzzz_0[i] = ta1_y_0_zzzzz_0[i] * fe_0 - ta1_y_0_zzzzz_1[i] * fe_0 + 5.0 * ta1_y_z_zzzz_0[i] * fe_0 -
                              5.0 * ta1_y_z_zzzz_1[i] * fe_0 + ta1_y_z_zzzzz_0[i] * pa_z[i] - ta1_y_z_zzzzz_1[i] * pc_z[i];
    }

    // Set up 252-273 components of targeted buffer : DH

    auto ta1_z_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 252);

    auto ta1_z_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 253);

    auto ta1_z_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 254);

    auto ta1_z_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 255);

    auto ta1_z_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 256);

    auto ta1_z_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 257);

    auto ta1_z_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 258);

    auto ta1_z_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 259);

    auto ta1_z_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 260);

    auto ta1_z_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 261);

    auto ta1_z_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 262);

    auto ta1_z_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 263);

    auto ta1_z_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 264);

    auto ta1_z_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 265);

    auto ta1_z_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 266);

    auto ta1_z_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 267);

    auto ta1_z_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 268);

    auto ta1_z_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 269);

    auto ta1_z_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 270);

    auto ta1_z_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 271);

    auto ta1_z_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 272);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_x_xxxx_0,   \
                             ta1_z_x_xxxx_1,   \
                             ta1_z_x_xxxxx_0,  \
                             ta1_z_x_xxxxx_1,  \
                             ta1_z_x_xxxxy_0,  \
                             ta1_z_x_xxxxy_1,  \
                             ta1_z_x_xxxxz_0,  \
                             ta1_z_x_xxxxz_1,  \
                             ta1_z_x_xxxy_0,   \
                             ta1_z_x_xxxy_1,   \
                             ta1_z_x_xxxyy_0,  \
                             ta1_z_x_xxxyy_1,  \
                             ta1_z_x_xxxyz_0,  \
                             ta1_z_x_xxxyz_1,  \
                             ta1_z_x_xxxz_0,   \
                             ta1_z_x_xxxz_1,   \
                             ta1_z_x_xxxzz_0,  \
                             ta1_z_x_xxxzz_1,  \
                             ta1_z_x_xxyy_0,   \
                             ta1_z_x_xxyy_1,   \
                             ta1_z_x_xxyyy_0,  \
                             ta1_z_x_xxyyy_1,  \
                             ta1_z_x_xxyyz_0,  \
                             ta1_z_x_xxyyz_1,  \
                             ta1_z_x_xxyz_0,   \
                             ta1_z_x_xxyz_1,   \
                             ta1_z_x_xxyzz_0,  \
                             ta1_z_x_xxyzz_1,  \
                             ta1_z_x_xxzz_0,   \
                             ta1_z_x_xxzz_1,   \
                             ta1_z_x_xxzzz_0,  \
                             ta1_z_x_xxzzz_1,  \
                             ta1_z_x_xyyy_0,   \
                             ta1_z_x_xyyy_1,   \
                             ta1_z_x_xyyyy_0,  \
                             ta1_z_x_xyyyy_1,  \
                             ta1_z_x_xyyyz_0,  \
                             ta1_z_x_xyyyz_1,  \
                             ta1_z_x_xyyz_0,   \
                             ta1_z_x_xyyz_1,   \
                             ta1_z_x_xyyzz_0,  \
                             ta1_z_x_xyyzz_1,  \
                             ta1_z_x_xyzz_0,   \
                             ta1_z_x_xyzz_1,   \
                             ta1_z_x_xyzzz_0,  \
                             ta1_z_x_xyzzz_1,  \
                             ta1_z_x_xzzz_0,   \
                             ta1_z_x_xzzz_1,   \
                             ta1_z_x_xzzzz_0,  \
                             ta1_z_x_xzzzz_1,  \
                             ta1_z_x_yyyy_0,   \
                             ta1_z_x_yyyy_1,   \
                             ta1_z_x_yyyyy_0,  \
                             ta1_z_x_yyyyy_1,  \
                             ta1_z_x_yyyyz_0,  \
                             ta1_z_x_yyyyz_1,  \
                             ta1_z_x_yyyz_0,   \
                             ta1_z_x_yyyz_1,   \
                             ta1_z_x_yyyzz_0,  \
                             ta1_z_x_yyyzz_1,  \
                             ta1_z_x_yyzz_0,   \
                             ta1_z_x_yyzz_1,   \
                             ta1_z_x_yyzzz_0,  \
                             ta1_z_x_yyzzz_1,  \
                             ta1_z_x_yzzz_0,   \
                             ta1_z_x_yzzz_1,   \
                             ta1_z_x_yzzzz_0,  \
                             ta1_z_x_yzzzz_1,  \
                             ta1_z_x_zzzz_0,   \
                             ta1_z_x_zzzz_1,   \
                             ta1_z_x_zzzzz_0,  \
                             ta1_z_x_zzzzz_1,  \
                             ta1_z_xx_xxxxx_0, \
                             ta1_z_xx_xxxxy_0, \
                             ta1_z_xx_xxxxz_0, \
                             ta1_z_xx_xxxyy_0, \
                             ta1_z_xx_xxxyz_0, \
                             ta1_z_xx_xxxzz_0, \
                             ta1_z_xx_xxyyy_0, \
                             ta1_z_xx_xxyyz_0, \
                             ta1_z_xx_xxyzz_0, \
                             ta1_z_xx_xxzzz_0, \
                             ta1_z_xx_xyyyy_0, \
                             ta1_z_xx_xyyyz_0, \
                             ta1_z_xx_xyyzz_0, \
                             ta1_z_xx_xyzzz_0, \
                             ta1_z_xx_xzzzz_0, \
                             ta1_z_xx_yyyyy_0, \
                             ta1_z_xx_yyyyz_0, \
                             ta1_z_xx_yyyzz_0, \
                             ta1_z_xx_yyzzz_0, \
                             ta1_z_xx_yzzzz_0, \
                             ta1_z_xx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xx_xxxxx_0[i] = ta1_z_0_xxxxx_0[i] * fe_0 - ta1_z_0_xxxxx_1[i] * fe_0 + 5.0 * ta1_z_x_xxxx_0[i] * fe_0 -
                              5.0 * ta1_z_x_xxxx_1[i] * fe_0 + ta1_z_x_xxxxx_0[i] * pa_x[i] - ta1_z_x_xxxxx_1[i] * pc_x[i];

        ta1_z_xx_xxxxy_0[i] = ta1_z_0_xxxxy_0[i] * fe_0 - ta1_z_0_xxxxy_1[i] * fe_0 + 4.0 * ta1_z_x_xxxy_0[i] * fe_0 -
                              4.0 * ta1_z_x_xxxy_1[i] * fe_0 + ta1_z_x_xxxxy_0[i] * pa_x[i] - ta1_z_x_xxxxy_1[i] * pc_x[i];

        ta1_z_xx_xxxxz_0[i] = ta1_z_0_xxxxz_0[i] * fe_0 - ta1_z_0_xxxxz_1[i] * fe_0 + 4.0 * ta1_z_x_xxxz_0[i] * fe_0 -
                              4.0 * ta1_z_x_xxxz_1[i] * fe_0 + ta1_z_x_xxxxz_0[i] * pa_x[i] - ta1_z_x_xxxxz_1[i] * pc_x[i];

        ta1_z_xx_xxxyy_0[i] = ta1_z_0_xxxyy_0[i] * fe_0 - ta1_z_0_xxxyy_1[i] * fe_0 + 3.0 * ta1_z_x_xxyy_0[i] * fe_0 -
                              3.0 * ta1_z_x_xxyy_1[i] * fe_0 + ta1_z_x_xxxyy_0[i] * pa_x[i] - ta1_z_x_xxxyy_1[i] * pc_x[i];

        ta1_z_xx_xxxyz_0[i] = ta1_z_0_xxxyz_0[i] * fe_0 - ta1_z_0_xxxyz_1[i] * fe_0 + 3.0 * ta1_z_x_xxyz_0[i] * fe_0 -
                              3.0 * ta1_z_x_xxyz_1[i] * fe_0 + ta1_z_x_xxxyz_0[i] * pa_x[i] - ta1_z_x_xxxyz_1[i] * pc_x[i];

        ta1_z_xx_xxxzz_0[i] = ta1_z_0_xxxzz_0[i] * fe_0 - ta1_z_0_xxxzz_1[i] * fe_0 + 3.0 * ta1_z_x_xxzz_0[i] * fe_0 -
                              3.0 * ta1_z_x_xxzz_1[i] * fe_0 + ta1_z_x_xxxzz_0[i] * pa_x[i] - ta1_z_x_xxxzz_1[i] * pc_x[i];

        ta1_z_xx_xxyyy_0[i] = ta1_z_0_xxyyy_0[i] * fe_0 - ta1_z_0_xxyyy_1[i] * fe_0 + 2.0 * ta1_z_x_xyyy_0[i] * fe_0 -
                              2.0 * ta1_z_x_xyyy_1[i] * fe_0 + ta1_z_x_xxyyy_0[i] * pa_x[i] - ta1_z_x_xxyyy_1[i] * pc_x[i];

        ta1_z_xx_xxyyz_0[i] = ta1_z_0_xxyyz_0[i] * fe_0 - ta1_z_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_x_xyyz_0[i] * fe_0 -
                              2.0 * ta1_z_x_xyyz_1[i] * fe_0 + ta1_z_x_xxyyz_0[i] * pa_x[i] - ta1_z_x_xxyyz_1[i] * pc_x[i];

        ta1_z_xx_xxyzz_0[i] = ta1_z_0_xxyzz_0[i] * fe_0 - ta1_z_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_x_xyzz_0[i] * fe_0 -
                              2.0 * ta1_z_x_xyzz_1[i] * fe_0 + ta1_z_x_xxyzz_0[i] * pa_x[i] - ta1_z_x_xxyzz_1[i] * pc_x[i];

        ta1_z_xx_xxzzz_0[i] = ta1_z_0_xxzzz_0[i] * fe_0 - ta1_z_0_xxzzz_1[i] * fe_0 + 2.0 * ta1_z_x_xzzz_0[i] * fe_0 -
                              2.0 * ta1_z_x_xzzz_1[i] * fe_0 + ta1_z_x_xxzzz_0[i] * pa_x[i] - ta1_z_x_xxzzz_1[i] * pc_x[i];

        ta1_z_xx_xyyyy_0[i] = ta1_z_0_xyyyy_0[i] * fe_0 - ta1_z_0_xyyyy_1[i] * fe_0 + ta1_z_x_yyyy_0[i] * fe_0 - ta1_z_x_yyyy_1[i] * fe_0 +
                              ta1_z_x_xyyyy_0[i] * pa_x[i] - ta1_z_x_xyyyy_1[i] * pc_x[i];

        ta1_z_xx_xyyyz_0[i] = ta1_z_0_xyyyz_0[i] * fe_0 - ta1_z_0_xyyyz_1[i] * fe_0 + ta1_z_x_yyyz_0[i] * fe_0 - ta1_z_x_yyyz_1[i] * fe_0 +
                              ta1_z_x_xyyyz_0[i] * pa_x[i] - ta1_z_x_xyyyz_1[i] * pc_x[i];

        ta1_z_xx_xyyzz_0[i] = ta1_z_0_xyyzz_0[i] * fe_0 - ta1_z_0_xyyzz_1[i] * fe_0 + ta1_z_x_yyzz_0[i] * fe_0 - ta1_z_x_yyzz_1[i] * fe_0 +
                              ta1_z_x_xyyzz_0[i] * pa_x[i] - ta1_z_x_xyyzz_1[i] * pc_x[i];

        ta1_z_xx_xyzzz_0[i] = ta1_z_0_xyzzz_0[i] * fe_0 - ta1_z_0_xyzzz_1[i] * fe_0 + ta1_z_x_yzzz_0[i] * fe_0 - ta1_z_x_yzzz_1[i] * fe_0 +
                              ta1_z_x_xyzzz_0[i] * pa_x[i] - ta1_z_x_xyzzz_1[i] * pc_x[i];

        ta1_z_xx_xzzzz_0[i] = ta1_z_0_xzzzz_0[i] * fe_0 - ta1_z_0_xzzzz_1[i] * fe_0 + ta1_z_x_zzzz_0[i] * fe_0 - ta1_z_x_zzzz_1[i] * fe_0 +
                              ta1_z_x_xzzzz_0[i] * pa_x[i] - ta1_z_x_xzzzz_1[i] * pc_x[i];

        ta1_z_xx_yyyyy_0[i] = ta1_z_0_yyyyy_0[i] * fe_0 - ta1_z_0_yyyyy_1[i] * fe_0 + ta1_z_x_yyyyy_0[i] * pa_x[i] - ta1_z_x_yyyyy_1[i] * pc_x[i];

        ta1_z_xx_yyyyz_0[i] = ta1_z_0_yyyyz_0[i] * fe_0 - ta1_z_0_yyyyz_1[i] * fe_0 + ta1_z_x_yyyyz_0[i] * pa_x[i] - ta1_z_x_yyyyz_1[i] * pc_x[i];

        ta1_z_xx_yyyzz_0[i] = ta1_z_0_yyyzz_0[i] * fe_0 - ta1_z_0_yyyzz_1[i] * fe_0 + ta1_z_x_yyyzz_0[i] * pa_x[i] - ta1_z_x_yyyzz_1[i] * pc_x[i];

        ta1_z_xx_yyzzz_0[i] = ta1_z_0_yyzzz_0[i] * fe_0 - ta1_z_0_yyzzz_1[i] * fe_0 + ta1_z_x_yyzzz_0[i] * pa_x[i] - ta1_z_x_yyzzz_1[i] * pc_x[i];

        ta1_z_xx_yzzzz_0[i] = ta1_z_0_yzzzz_0[i] * fe_0 - ta1_z_0_yzzzz_1[i] * fe_0 + ta1_z_x_yzzzz_0[i] * pa_x[i] - ta1_z_x_yzzzz_1[i] * pc_x[i];

        ta1_z_xx_zzzzz_0[i] = ta1_z_0_zzzzz_0[i] * fe_0 - ta1_z_0_zzzzz_1[i] * fe_0 + ta1_z_x_zzzzz_0[i] * pa_x[i] - ta1_z_x_zzzzz_1[i] * pc_x[i];
    }

    // Set up 273-294 components of targeted buffer : DH

    auto ta1_z_xy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 273);

    auto ta1_z_xy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 274);

    auto ta1_z_xy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 275);

    auto ta1_z_xy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 276);

    auto ta1_z_xy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 277);

    auto ta1_z_xy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 278);

    auto ta1_z_xy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 279);

    auto ta1_z_xy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 280);

    auto ta1_z_xy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 281);

    auto ta1_z_xy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 282);

    auto ta1_z_xy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 283);

    auto ta1_z_xy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 284);

    auto ta1_z_xy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 285);

    auto ta1_z_xy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 286);

    auto ta1_z_xy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 287);

    auto ta1_z_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 288);

    auto ta1_z_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 289);

    auto ta1_z_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 290);

    auto ta1_z_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 291);

    auto ta1_z_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 292);

    auto ta1_z_xy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 293);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_x_xxxxx_0,  \
                             ta1_z_x_xxxxx_1,  \
                             ta1_z_x_xxxxz_0,  \
                             ta1_z_x_xxxxz_1,  \
                             ta1_z_x_xxxzz_0,  \
                             ta1_z_x_xxxzz_1,  \
                             ta1_z_x_xxzzz_0,  \
                             ta1_z_x_xxzzz_1,  \
                             ta1_z_x_xzzzz_0,  \
                             ta1_z_x_xzzzz_1,  \
                             ta1_z_xy_xxxxx_0, \
                             ta1_z_xy_xxxxy_0, \
                             ta1_z_xy_xxxxz_0, \
                             ta1_z_xy_xxxyy_0, \
                             ta1_z_xy_xxxyz_0, \
                             ta1_z_xy_xxxzz_0, \
                             ta1_z_xy_xxyyy_0, \
                             ta1_z_xy_xxyyz_0, \
                             ta1_z_xy_xxyzz_0, \
                             ta1_z_xy_xxzzz_0, \
                             ta1_z_xy_xyyyy_0, \
                             ta1_z_xy_xyyyz_0, \
                             ta1_z_xy_xyyzz_0, \
                             ta1_z_xy_xyzzz_0, \
                             ta1_z_xy_xzzzz_0, \
                             ta1_z_xy_yyyyy_0, \
                             ta1_z_xy_yyyyz_0, \
                             ta1_z_xy_yyyzz_0, \
                             ta1_z_xy_yyzzz_0, \
                             ta1_z_xy_yzzzz_0, \
                             ta1_z_xy_zzzzz_0, \
                             ta1_z_y_xxxxy_0,  \
                             ta1_z_y_xxxxy_1,  \
                             ta1_z_y_xxxy_0,   \
                             ta1_z_y_xxxy_1,   \
                             ta1_z_y_xxxyy_0,  \
                             ta1_z_y_xxxyy_1,  \
                             ta1_z_y_xxxyz_0,  \
                             ta1_z_y_xxxyz_1,  \
                             ta1_z_y_xxyy_0,   \
                             ta1_z_y_xxyy_1,   \
                             ta1_z_y_xxyyy_0,  \
                             ta1_z_y_xxyyy_1,  \
                             ta1_z_y_xxyyz_0,  \
                             ta1_z_y_xxyyz_1,  \
                             ta1_z_y_xxyz_0,   \
                             ta1_z_y_xxyz_1,   \
                             ta1_z_y_xxyzz_0,  \
                             ta1_z_y_xxyzz_1,  \
                             ta1_z_y_xyyy_0,   \
                             ta1_z_y_xyyy_1,   \
                             ta1_z_y_xyyyy_0,  \
                             ta1_z_y_xyyyy_1,  \
                             ta1_z_y_xyyyz_0,  \
                             ta1_z_y_xyyyz_1,  \
                             ta1_z_y_xyyz_0,   \
                             ta1_z_y_xyyz_1,   \
                             ta1_z_y_xyyzz_0,  \
                             ta1_z_y_xyyzz_1,  \
                             ta1_z_y_xyzz_0,   \
                             ta1_z_y_xyzz_1,   \
                             ta1_z_y_xyzzz_0,  \
                             ta1_z_y_xyzzz_1,  \
                             ta1_z_y_yyyy_0,   \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyyy_0,  \
                             ta1_z_y_yyyyy_1,  \
                             ta1_z_y_yyyyz_0,  \
                             ta1_z_y_yyyyz_1,  \
                             ta1_z_y_yyyz_0,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyyzz_0,  \
                             ta1_z_y_yyyzz_1,  \
                             ta1_z_y_yyzz_0,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yyzzz_0,  \
                             ta1_z_y_yyzzz_1,  \
                             ta1_z_y_yzzz_0,   \
                             ta1_z_y_yzzz_1,   \
                             ta1_z_y_yzzzz_0,  \
                             ta1_z_y_yzzzz_1,  \
                             ta1_z_y_zzzzz_0,  \
                             ta1_z_y_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xy_xxxxx_0[i] = ta1_z_x_xxxxx_0[i] * pa_y[i] - ta1_z_x_xxxxx_1[i] * pc_y[i];

        ta1_z_xy_xxxxy_0[i] =
            4.0 * ta1_z_y_xxxy_0[i] * fe_0 - 4.0 * ta1_z_y_xxxy_1[i] * fe_0 + ta1_z_y_xxxxy_0[i] * pa_x[i] - ta1_z_y_xxxxy_1[i] * pc_x[i];

        ta1_z_xy_xxxxz_0[i] = ta1_z_x_xxxxz_0[i] * pa_y[i] - ta1_z_x_xxxxz_1[i] * pc_y[i];

        ta1_z_xy_xxxyy_0[i] =
            3.0 * ta1_z_y_xxyy_0[i] * fe_0 - 3.0 * ta1_z_y_xxyy_1[i] * fe_0 + ta1_z_y_xxxyy_0[i] * pa_x[i] - ta1_z_y_xxxyy_1[i] * pc_x[i];

        ta1_z_xy_xxxyz_0[i] =
            3.0 * ta1_z_y_xxyz_0[i] * fe_0 - 3.0 * ta1_z_y_xxyz_1[i] * fe_0 + ta1_z_y_xxxyz_0[i] * pa_x[i] - ta1_z_y_xxxyz_1[i] * pc_x[i];

        ta1_z_xy_xxxzz_0[i] = ta1_z_x_xxxzz_0[i] * pa_y[i] - ta1_z_x_xxxzz_1[i] * pc_y[i];

        ta1_z_xy_xxyyy_0[i] =
            2.0 * ta1_z_y_xyyy_0[i] * fe_0 - 2.0 * ta1_z_y_xyyy_1[i] * fe_0 + ta1_z_y_xxyyy_0[i] * pa_x[i] - ta1_z_y_xxyyy_1[i] * pc_x[i];

        ta1_z_xy_xxyyz_0[i] =
            2.0 * ta1_z_y_xyyz_0[i] * fe_0 - 2.0 * ta1_z_y_xyyz_1[i] * fe_0 + ta1_z_y_xxyyz_0[i] * pa_x[i] - ta1_z_y_xxyyz_1[i] * pc_x[i];

        ta1_z_xy_xxyzz_0[i] =
            2.0 * ta1_z_y_xyzz_0[i] * fe_0 - 2.0 * ta1_z_y_xyzz_1[i] * fe_0 + ta1_z_y_xxyzz_0[i] * pa_x[i] - ta1_z_y_xxyzz_1[i] * pc_x[i];

        ta1_z_xy_xxzzz_0[i] = ta1_z_x_xxzzz_0[i] * pa_y[i] - ta1_z_x_xxzzz_1[i] * pc_y[i];

        ta1_z_xy_xyyyy_0[i] = ta1_z_y_yyyy_0[i] * fe_0 - ta1_z_y_yyyy_1[i] * fe_0 + ta1_z_y_xyyyy_0[i] * pa_x[i] - ta1_z_y_xyyyy_1[i] * pc_x[i];

        ta1_z_xy_xyyyz_0[i] = ta1_z_y_yyyz_0[i] * fe_0 - ta1_z_y_yyyz_1[i] * fe_0 + ta1_z_y_xyyyz_0[i] * pa_x[i] - ta1_z_y_xyyyz_1[i] * pc_x[i];

        ta1_z_xy_xyyzz_0[i] = ta1_z_y_yyzz_0[i] * fe_0 - ta1_z_y_yyzz_1[i] * fe_0 + ta1_z_y_xyyzz_0[i] * pa_x[i] - ta1_z_y_xyyzz_1[i] * pc_x[i];

        ta1_z_xy_xyzzz_0[i] = ta1_z_y_yzzz_0[i] * fe_0 - ta1_z_y_yzzz_1[i] * fe_0 + ta1_z_y_xyzzz_0[i] * pa_x[i] - ta1_z_y_xyzzz_1[i] * pc_x[i];

        ta1_z_xy_xzzzz_0[i] = ta1_z_x_xzzzz_0[i] * pa_y[i] - ta1_z_x_xzzzz_1[i] * pc_y[i];

        ta1_z_xy_yyyyy_0[i] = ta1_z_y_yyyyy_0[i] * pa_x[i] - ta1_z_y_yyyyy_1[i] * pc_x[i];

        ta1_z_xy_yyyyz_0[i] = ta1_z_y_yyyyz_0[i] * pa_x[i] - ta1_z_y_yyyyz_1[i] * pc_x[i];

        ta1_z_xy_yyyzz_0[i] = ta1_z_y_yyyzz_0[i] * pa_x[i] - ta1_z_y_yyyzz_1[i] * pc_x[i];

        ta1_z_xy_yyzzz_0[i] = ta1_z_y_yyzzz_0[i] * pa_x[i] - ta1_z_y_yyzzz_1[i] * pc_x[i];

        ta1_z_xy_yzzzz_0[i] = ta1_z_y_yzzzz_0[i] * pa_x[i] - ta1_z_y_yzzzz_1[i] * pc_x[i];

        ta1_z_xy_zzzzz_0[i] = ta1_z_y_zzzzz_0[i] * pa_x[i] - ta1_z_y_zzzzz_1[i] * pc_x[i];
    }

    // Set up 294-315 components of targeted buffer : DH

    auto ta1_z_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 294);

    auto ta1_z_xz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 295);

    auto ta1_z_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 296);

    auto ta1_z_xz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 297);

    auto ta1_z_xz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 298);

    auto ta1_z_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 299);

    auto ta1_z_xz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 300);

    auto ta1_z_xz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 301);

    auto ta1_z_xz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 302);

    auto ta1_z_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 303);

    auto ta1_z_xz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 304);

    auto ta1_z_xz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 305);

    auto ta1_z_xz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 306);

    auto ta1_z_xz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 307);

    auto ta1_z_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 308);

    auto ta1_z_xz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 309);

    auto ta1_z_xz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 310);

    auto ta1_z_xz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 311);

    auto ta1_z_xz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 312);

    auto ta1_z_xz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 313);

    auto ta1_z_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 314);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_x_xxxxx_0,  \
                             ta1_z_x_xxxxx_1,  \
                             ta1_z_x_xxxxy_0,  \
                             ta1_z_x_xxxxy_1,  \
                             ta1_z_x_xxxyy_0,  \
                             ta1_z_x_xxxyy_1,  \
                             ta1_z_x_xxyyy_0,  \
                             ta1_z_x_xxyyy_1,  \
                             ta1_z_x_xyyyy_0,  \
                             ta1_z_x_xyyyy_1,  \
                             ta1_z_xz_xxxxx_0, \
                             ta1_z_xz_xxxxy_0, \
                             ta1_z_xz_xxxxz_0, \
                             ta1_z_xz_xxxyy_0, \
                             ta1_z_xz_xxxyz_0, \
                             ta1_z_xz_xxxzz_0, \
                             ta1_z_xz_xxyyy_0, \
                             ta1_z_xz_xxyyz_0, \
                             ta1_z_xz_xxyzz_0, \
                             ta1_z_xz_xxzzz_0, \
                             ta1_z_xz_xyyyy_0, \
                             ta1_z_xz_xyyyz_0, \
                             ta1_z_xz_xyyzz_0, \
                             ta1_z_xz_xyzzz_0, \
                             ta1_z_xz_xzzzz_0, \
                             ta1_z_xz_yyyyy_0, \
                             ta1_z_xz_yyyyz_0, \
                             ta1_z_xz_yyyzz_0, \
                             ta1_z_xz_yyzzz_0, \
                             ta1_z_xz_yzzzz_0, \
                             ta1_z_xz_zzzzz_0, \
                             ta1_z_z_xxxxz_0,  \
                             ta1_z_z_xxxxz_1,  \
                             ta1_z_z_xxxyz_0,  \
                             ta1_z_z_xxxyz_1,  \
                             ta1_z_z_xxxz_0,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxxzz_0,  \
                             ta1_z_z_xxxzz_1,  \
                             ta1_z_z_xxyyz_0,  \
                             ta1_z_z_xxyyz_1,  \
                             ta1_z_z_xxyz_0,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxyzz_0,  \
                             ta1_z_z_xxyzz_1,  \
                             ta1_z_z_xxzz_0,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xxzzz_0,  \
                             ta1_z_z_xxzzz_1,  \
                             ta1_z_z_xyyyz_0,  \
                             ta1_z_z_xyyyz_1,  \
                             ta1_z_z_xyyz_0,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyyzz_0,  \
                             ta1_z_z_xyyzz_1,  \
                             ta1_z_z_xyzz_0,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xyzzz_0,  \
                             ta1_z_z_xyzzz_1,  \
                             ta1_z_z_xzzz_0,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_xzzzz_0,  \
                             ta1_z_z_xzzzz_1,  \
                             ta1_z_z_yyyyy_0,  \
                             ta1_z_z_yyyyy_1,  \
                             ta1_z_z_yyyyz_0,  \
                             ta1_z_z_yyyyz_1,  \
                             ta1_z_z_yyyz_0,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyyzz_0,  \
                             ta1_z_z_yyyzz_1,  \
                             ta1_z_z_yyzz_0,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yyzzz_0,  \
                             ta1_z_z_yyzzz_1,  \
                             ta1_z_z_yzzz_0,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_yzzzz_0,  \
                             ta1_z_z_yzzzz_1,  \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta1_z_z_zzzzz_0,  \
                             ta1_z_z_zzzzz_1,  \
                             ta_x_xxxxx_1,     \
                             ta_x_xxxxy_1,     \
                             ta_x_xxxyy_1,     \
                             ta_x_xxyyy_1,     \
                             ta_x_xyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xz_xxxxx_0[i] = ta_x_xxxxx_1[i] + ta1_z_x_xxxxx_0[i] * pa_z[i] - ta1_z_x_xxxxx_1[i] * pc_z[i];

        ta1_z_xz_xxxxy_0[i] = ta_x_xxxxy_1[i] + ta1_z_x_xxxxy_0[i] * pa_z[i] - ta1_z_x_xxxxy_1[i] * pc_z[i];

        ta1_z_xz_xxxxz_0[i] =
            4.0 * ta1_z_z_xxxz_0[i] * fe_0 - 4.0 * ta1_z_z_xxxz_1[i] * fe_0 + ta1_z_z_xxxxz_0[i] * pa_x[i] - ta1_z_z_xxxxz_1[i] * pc_x[i];

        ta1_z_xz_xxxyy_0[i] = ta_x_xxxyy_1[i] + ta1_z_x_xxxyy_0[i] * pa_z[i] - ta1_z_x_xxxyy_1[i] * pc_z[i];

        ta1_z_xz_xxxyz_0[i] =
            3.0 * ta1_z_z_xxyz_0[i] * fe_0 - 3.0 * ta1_z_z_xxyz_1[i] * fe_0 + ta1_z_z_xxxyz_0[i] * pa_x[i] - ta1_z_z_xxxyz_1[i] * pc_x[i];

        ta1_z_xz_xxxzz_0[i] =
            3.0 * ta1_z_z_xxzz_0[i] * fe_0 - 3.0 * ta1_z_z_xxzz_1[i] * fe_0 + ta1_z_z_xxxzz_0[i] * pa_x[i] - ta1_z_z_xxxzz_1[i] * pc_x[i];

        ta1_z_xz_xxyyy_0[i] = ta_x_xxyyy_1[i] + ta1_z_x_xxyyy_0[i] * pa_z[i] - ta1_z_x_xxyyy_1[i] * pc_z[i];

        ta1_z_xz_xxyyz_0[i] =
            2.0 * ta1_z_z_xyyz_0[i] * fe_0 - 2.0 * ta1_z_z_xyyz_1[i] * fe_0 + ta1_z_z_xxyyz_0[i] * pa_x[i] - ta1_z_z_xxyyz_1[i] * pc_x[i];

        ta1_z_xz_xxyzz_0[i] =
            2.0 * ta1_z_z_xyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzz_1[i] * fe_0 + ta1_z_z_xxyzz_0[i] * pa_x[i] - ta1_z_z_xxyzz_1[i] * pc_x[i];

        ta1_z_xz_xxzzz_0[i] =
            2.0 * ta1_z_z_xzzz_0[i] * fe_0 - 2.0 * ta1_z_z_xzzz_1[i] * fe_0 + ta1_z_z_xxzzz_0[i] * pa_x[i] - ta1_z_z_xxzzz_1[i] * pc_x[i];

        ta1_z_xz_xyyyy_0[i] = ta_x_xyyyy_1[i] + ta1_z_x_xyyyy_0[i] * pa_z[i] - ta1_z_x_xyyyy_1[i] * pc_z[i];

        ta1_z_xz_xyyyz_0[i] = ta1_z_z_yyyz_0[i] * fe_0 - ta1_z_z_yyyz_1[i] * fe_0 + ta1_z_z_xyyyz_0[i] * pa_x[i] - ta1_z_z_xyyyz_1[i] * pc_x[i];

        ta1_z_xz_xyyzz_0[i] = ta1_z_z_yyzz_0[i] * fe_0 - ta1_z_z_yyzz_1[i] * fe_0 + ta1_z_z_xyyzz_0[i] * pa_x[i] - ta1_z_z_xyyzz_1[i] * pc_x[i];

        ta1_z_xz_xyzzz_0[i] = ta1_z_z_yzzz_0[i] * fe_0 - ta1_z_z_yzzz_1[i] * fe_0 + ta1_z_z_xyzzz_0[i] * pa_x[i] - ta1_z_z_xyzzz_1[i] * pc_x[i];

        ta1_z_xz_xzzzz_0[i] = ta1_z_z_zzzz_0[i] * fe_0 - ta1_z_z_zzzz_1[i] * fe_0 + ta1_z_z_xzzzz_0[i] * pa_x[i] - ta1_z_z_xzzzz_1[i] * pc_x[i];

        ta1_z_xz_yyyyy_0[i] = ta1_z_z_yyyyy_0[i] * pa_x[i] - ta1_z_z_yyyyy_1[i] * pc_x[i];

        ta1_z_xz_yyyyz_0[i] = ta1_z_z_yyyyz_0[i] * pa_x[i] - ta1_z_z_yyyyz_1[i] * pc_x[i];

        ta1_z_xz_yyyzz_0[i] = ta1_z_z_yyyzz_0[i] * pa_x[i] - ta1_z_z_yyyzz_1[i] * pc_x[i];

        ta1_z_xz_yyzzz_0[i] = ta1_z_z_yyzzz_0[i] * pa_x[i] - ta1_z_z_yyzzz_1[i] * pc_x[i];

        ta1_z_xz_yzzzz_0[i] = ta1_z_z_yzzzz_0[i] * pa_x[i] - ta1_z_z_yzzzz_1[i] * pc_x[i];

        ta1_z_xz_zzzzz_0[i] = ta1_z_z_zzzzz_0[i] * pa_x[i] - ta1_z_z_zzzzz_1[i] * pc_x[i];
    }

    // Set up 315-336 components of targeted buffer : DH

    auto ta1_z_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 315);

    auto ta1_z_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 316);

    auto ta1_z_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 317);

    auto ta1_z_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 318);

    auto ta1_z_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 319);

    auto ta1_z_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 320);

    auto ta1_z_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 321);

    auto ta1_z_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 322);

    auto ta1_z_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 323);

    auto ta1_z_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 324);

    auto ta1_z_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 325);

    auto ta1_z_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 326);

    auto ta1_z_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 327);

    auto ta1_z_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 328);

    auto ta1_z_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 329);

    auto ta1_z_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 330);

    auto ta1_z_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 331);

    auto ta1_z_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 332);

    auto ta1_z_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 333);

    auto ta1_z_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 334);

    auto ta1_z_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 335);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_y_xxxx_0,   \
                             ta1_z_y_xxxx_1,   \
                             ta1_z_y_xxxxx_0,  \
                             ta1_z_y_xxxxx_1,  \
                             ta1_z_y_xxxxy_0,  \
                             ta1_z_y_xxxxy_1,  \
                             ta1_z_y_xxxxz_0,  \
                             ta1_z_y_xxxxz_1,  \
                             ta1_z_y_xxxy_0,   \
                             ta1_z_y_xxxy_1,   \
                             ta1_z_y_xxxyy_0,  \
                             ta1_z_y_xxxyy_1,  \
                             ta1_z_y_xxxyz_0,  \
                             ta1_z_y_xxxyz_1,  \
                             ta1_z_y_xxxz_0,   \
                             ta1_z_y_xxxz_1,   \
                             ta1_z_y_xxxzz_0,  \
                             ta1_z_y_xxxzz_1,  \
                             ta1_z_y_xxyy_0,   \
                             ta1_z_y_xxyy_1,   \
                             ta1_z_y_xxyyy_0,  \
                             ta1_z_y_xxyyy_1,  \
                             ta1_z_y_xxyyz_0,  \
                             ta1_z_y_xxyyz_1,  \
                             ta1_z_y_xxyz_0,   \
                             ta1_z_y_xxyz_1,   \
                             ta1_z_y_xxyzz_0,  \
                             ta1_z_y_xxyzz_1,  \
                             ta1_z_y_xxzz_0,   \
                             ta1_z_y_xxzz_1,   \
                             ta1_z_y_xxzzz_0,  \
                             ta1_z_y_xxzzz_1,  \
                             ta1_z_y_xyyy_0,   \
                             ta1_z_y_xyyy_1,   \
                             ta1_z_y_xyyyy_0,  \
                             ta1_z_y_xyyyy_1,  \
                             ta1_z_y_xyyyz_0,  \
                             ta1_z_y_xyyyz_1,  \
                             ta1_z_y_xyyz_0,   \
                             ta1_z_y_xyyz_1,   \
                             ta1_z_y_xyyzz_0,  \
                             ta1_z_y_xyyzz_1,  \
                             ta1_z_y_xyzz_0,   \
                             ta1_z_y_xyzz_1,   \
                             ta1_z_y_xyzzz_0,  \
                             ta1_z_y_xyzzz_1,  \
                             ta1_z_y_xzzz_0,   \
                             ta1_z_y_xzzz_1,   \
                             ta1_z_y_xzzzz_0,  \
                             ta1_z_y_xzzzz_1,  \
                             ta1_z_y_yyyy_0,   \
                             ta1_z_y_yyyy_1,   \
                             ta1_z_y_yyyyy_0,  \
                             ta1_z_y_yyyyy_1,  \
                             ta1_z_y_yyyyz_0,  \
                             ta1_z_y_yyyyz_1,  \
                             ta1_z_y_yyyz_0,   \
                             ta1_z_y_yyyz_1,   \
                             ta1_z_y_yyyzz_0,  \
                             ta1_z_y_yyyzz_1,  \
                             ta1_z_y_yyzz_0,   \
                             ta1_z_y_yyzz_1,   \
                             ta1_z_y_yyzzz_0,  \
                             ta1_z_y_yyzzz_1,  \
                             ta1_z_y_yzzz_0,   \
                             ta1_z_y_yzzz_1,   \
                             ta1_z_y_yzzzz_0,  \
                             ta1_z_y_yzzzz_1,  \
                             ta1_z_y_zzzz_0,   \
                             ta1_z_y_zzzz_1,   \
                             ta1_z_y_zzzzz_0,  \
                             ta1_z_y_zzzzz_1,  \
                             ta1_z_yy_xxxxx_0, \
                             ta1_z_yy_xxxxy_0, \
                             ta1_z_yy_xxxxz_0, \
                             ta1_z_yy_xxxyy_0, \
                             ta1_z_yy_xxxyz_0, \
                             ta1_z_yy_xxxzz_0, \
                             ta1_z_yy_xxyyy_0, \
                             ta1_z_yy_xxyyz_0, \
                             ta1_z_yy_xxyzz_0, \
                             ta1_z_yy_xxzzz_0, \
                             ta1_z_yy_xyyyy_0, \
                             ta1_z_yy_xyyyz_0, \
                             ta1_z_yy_xyyzz_0, \
                             ta1_z_yy_xyzzz_0, \
                             ta1_z_yy_xzzzz_0, \
                             ta1_z_yy_yyyyy_0, \
                             ta1_z_yy_yyyyz_0, \
                             ta1_z_yy_yyyzz_0, \
                             ta1_z_yy_yyzzz_0, \
                             ta1_z_yy_yzzzz_0, \
                             ta1_z_yy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yy_xxxxx_0[i] = ta1_z_0_xxxxx_0[i] * fe_0 - ta1_z_0_xxxxx_1[i] * fe_0 + ta1_z_y_xxxxx_0[i] * pa_y[i] - ta1_z_y_xxxxx_1[i] * pc_y[i];

        ta1_z_yy_xxxxy_0[i] = ta1_z_0_xxxxy_0[i] * fe_0 - ta1_z_0_xxxxy_1[i] * fe_0 + ta1_z_y_xxxx_0[i] * fe_0 - ta1_z_y_xxxx_1[i] * fe_0 +
                              ta1_z_y_xxxxy_0[i] * pa_y[i] - ta1_z_y_xxxxy_1[i] * pc_y[i];

        ta1_z_yy_xxxxz_0[i] = ta1_z_0_xxxxz_0[i] * fe_0 - ta1_z_0_xxxxz_1[i] * fe_0 + ta1_z_y_xxxxz_0[i] * pa_y[i] - ta1_z_y_xxxxz_1[i] * pc_y[i];

        ta1_z_yy_xxxyy_0[i] = ta1_z_0_xxxyy_0[i] * fe_0 - ta1_z_0_xxxyy_1[i] * fe_0 + 2.0 * ta1_z_y_xxxy_0[i] * fe_0 -
                              2.0 * ta1_z_y_xxxy_1[i] * fe_0 + ta1_z_y_xxxyy_0[i] * pa_y[i] - ta1_z_y_xxxyy_1[i] * pc_y[i];

        ta1_z_yy_xxxyz_0[i] = ta1_z_0_xxxyz_0[i] * fe_0 - ta1_z_0_xxxyz_1[i] * fe_0 + ta1_z_y_xxxz_0[i] * fe_0 - ta1_z_y_xxxz_1[i] * fe_0 +
                              ta1_z_y_xxxyz_0[i] * pa_y[i] - ta1_z_y_xxxyz_1[i] * pc_y[i];

        ta1_z_yy_xxxzz_0[i] = ta1_z_0_xxxzz_0[i] * fe_0 - ta1_z_0_xxxzz_1[i] * fe_0 + ta1_z_y_xxxzz_0[i] * pa_y[i] - ta1_z_y_xxxzz_1[i] * pc_y[i];

        ta1_z_yy_xxyyy_0[i] = ta1_z_0_xxyyy_0[i] * fe_0 - ta1_z_0_xxyyy_1[i] * fe_0 + 3.0 * ta1_z_y_xxyy_0[i] * fe_0 -
                              3.0 * ta1_z_y_xxyy_1[i] * fe_0 + ta1_z_y_xxyyy_0[i] * pa_y[i] - ta1_z_y_xxyyy_1[i] * pc_y[i];

        ta1_z_yy_xxyyz_0[i] = ta1_z_0_xxyyz_0[i] * fe_0 - ta1_z_0_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_y_xxyz_0[i] * fe_0 -
                              2.0 * ta1_z_y_xxyz_1[i] * fe_0 + ta1_z_y_xxyyz_0[i] * pa_y[i] - ta1_z_y_xxyyz_1[i] * pc_y[i];

        ta1_z_yy_xxyzz_0[i] = ta1_z_0_xxyzz_0[i] * fe_0 - ta1_z_0_xxyzz_1[i] * fe_0 + ta1_z_y_xxzz_0[i] * fe_0 - ta1_z_y_xxzz_1[i] * fe_0 +
                              ta1_z_y_xxyzz_0[i] * pa_y[i] - ta1_z_y_xxyzz_1[i] * pc_y[i];

        ta1_z_yy_xxzzz_0[i] = ta1_z_0_xxzzz_0[i] * fe_0 - ta1_z_0_xxzzz_1[i] * fe_0 + ta1_z_y_xxzzz_0[i] * pa_y[i] - ta1_z_y_xxzzz_1[i] * pc_y[i];

        ta1_z_yy_xyyyy_0[i] = ta1_z_0_xyyyy_0[i] * fe_0 - ta1_z_0_xyyyy_1[i] * fe_0 + 4.0 * ta1_z_y_xyyy_0[i] * fe_0 -
                              4.0 * ta1_z_y_xyyy_1[i] * fe_0 + ta1_z_y_xyyyy_0[i] * pa_y[i] - ta1_z_y_xyyyy_1[i] * pc_y[i];

        ta1_z_yy_xyyyz_0[i] = ta1_z_0_xyyyz_0[i] * fe_0 - ta1_z_0_xyyyz_1[i] * fe_0 + 3.0 * ta1_z_y_xyyz_0[i] * fe_0 -
                              3.0 * ta1_z_y_xyyz_1[i] * fe_0 + ta1_z_y_xyyyz_0[i] * pa_y[i] - ta1_z_y_xyyyz_1[i] * pc_y[i];

        ta1_z_yy_xyyzz_0[i] = ta1_z_0_xyyzz_0[i] * fe_0 - ta1_z_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_y_xyzz_0[i] * fe_0 -
                              2.0 * ta1_z_y_xyzz_1[i] * fe_0 + ta1_z_y_xyyzz_0[i] * pa_y[i] - ta1_z_y_xyyzz_1[i] * pc_y[i];

        ta1_z_yy_xyzzz_0[i] = ta1_z_0_xyzzz_0[i] * fe_0 - ta1_z_0_xyzzz_1[i] * fe_0 + ta1_z_y_xzzz_0[i] * fe_0 - ta1_z_y_xzzz_1[i] * fe_0 +
                              ta1_z_y_xyzzz_0[i] * pa_y[i] - ta1_z_y_xyzzz_1[i] * pc_y[i];

        ta1_z_yy_xzzzz_0[i] = ta1_z_0_xzzzz_0[i] * fe_0 - ta1_z_0_xzzzz_1[i] * fe_0 + ta1_z_y_xzzzz_0[i] * pa_y[i] - ta1_z_y_xzzzz_1[i] * pc_y[i];

        ta1_z_yy_yyyyy_0[i] = ta1_z_0_yyyyy_0[i] * fe_0 - ta1_z_0_yyyyy_1[i] * fe_0 + 5.0 * ta1_z_y_yyyy_0[i] * fe_0 -
                              5.0 * ta1_z_y_yyyy_1[i] * fe_0 + ta1_z_y_yyyyy_0[i] * pa_y[i] - ta1_z_y_yyyyy_1[i] * pc_y[i];

        ta1_z_yy_yyyyz_0[i] = ta1_z_0_yyyyz_0[i] * fe_0 - ta1_z_0_yyyyz_1[i] * fe_0 + 4.0 * ta1_z_y_yyyz_0[i] * fe_0 -
                              4.0 * ta1_z_y_yyyz_1[i] * fe_0 + ta1_z_y_yyyyz_0[i] * pa_y[i] - ta1_z_y_yyyyz_1[i] * pc_y[i];

        ta1_z_yy_yyyzz_0[i] = ta1_z_0_yyyzz_0[i] * fe_0 - ta1_z_0_yyyzz_1[i] * fe_0 + 3.0 * ta1_z_y_yyzz_0[i] * fe_0 -
                              3.0 * ta1_z_y_yyzz_1[i] * fe_0 + ta1_z_y_yyyzz_0[i] * pa_y[i] - ta1_z_y_yyyzz_1[i] * pc_y[i];

        ta1_z_yy_yyzzz_0[i] = ta1_z_0_yyzzz_0[i] * fe_0 - ta1_z_0_yyzzz_1[i] * fe_0 + 2.0 * ta1_z_y_yzzz_0[i] * fe_0 -
                              2.0 * ta1_z_y_yzzz_1[i] * fe_0 + ta1_z_y_yyzzz_0[i] * pa_y[i] - ta1_z_y_yyzzz_1[i] * pc_y[i];

        ta1_z_yy_yzzzz_0[i] = ta1_z_0_yzzzz_0[i] * fe_0 - ta1_z_0_yzzzz_1[i] * fe_0 + ta1_z_y_zzzz_0[i] * fe_0 - ta1_z_y_zzzz_1[i] * fe_0 +
                              ta1_z_y_yzzzz_0[i] * pa_y[i] - ta1_z_y_yzzzz_1[i] * pc_y[i];

        ta1_z_yy_zzzzz_0[i] = ta1_z_0_zzzzz_0[i] * fe_0 - ta1_z_0_zzzzz_1[i] * fe_0 + ta1_z_y_zzzzz_0[i] * pa_y[i] - ta1_z_y_zzzzz_1[i] * pc_y[i];
    }

    // Set up 336-357 components of targeted buffer : DH

    auto ta1_z_yz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 336);

    auto ta1_z_yz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 337);

    auto ta1_z_yz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 338);

    auto ta1_z_yz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 339);

    auto ta1_z_yz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 340);

    auto ta1_z_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 341);

    auto ta1_z_yz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 342);

    auto ta1_z_yz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 343);

    auto ta1_z_yz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 344);

    auto ta1_z_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 345);

    auto ta1_z_yz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 346);

    auto ta1_z_yz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 347);

    auto ta1_z_yz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 348);

    auto ta1_z_yz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 349);

    auto ta1_z_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 350);

    auto ta1_z_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 351);

    auto ta1_z_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 352);

    auto ta1_z_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 353);

    auto ta1_z_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 354);

    auto ta1_z_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 355);

    auto ta1_z_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 356);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_y_xxxxy_0,  \
                             ta1_z_y_xxxxy_1,  \
                             ta1_z_y_xxxyy_0,  \
                             ta1_z_y_xxxyy_1,  \
                             ta1_z_y_xxyyy_0,  \
                             ta1_z_y_xxyyy_1,  \
                             ta1_z_y_xyyyy_0,  \
                             ta1_z_y_xyyyy_1,  \
                             ta1_z_y_yyyyy_0,  \
                             ta1_z_y_yyyyy_1,  \
                             ta1_z_yz_xxxxx_0, \
                             ta1_z_yz_xxxxy_0, \
                             ta1_z_yz_xxxxz_0, \
                             ta1_z_yz_xxxyy_0, \
                             ta1_z_yz_xxxyz_0, \
                             ta1_z_yz_xxxzz_0, \
                             ta1_z_yz_xxyyy_0, \
                             ta1_z_yz_xxyyz_0, \
                             ta1_z_yz_xxyzz_0, \
                             ta1_z_yz_xxzzz_0, \
                             ta1_z_yz_xyyyy_0, \
                             ta1_z_yz_xyyyz_0, \
                             ta1_z_yz_xyyzz_0, \
                             ta1_z_yz_xyzzz_0, \
                             ta1_z_yz_xzzzz_0, \
                             ta1_z_yz_yyyyy_0, \
                             ta1_z_yz_yyyyz_0, \
                             ta1_z_yz_yyyzz_0, \
                             ta1_z_yz_yyzzz_0, \
                             ta1_z_yz_yzzzz_0, \
                             ta1_z_yz_zzzzz_0, \
                             ta1_z_z_xxxxx_0,  \
                             ta1_z_z_xxxxx_1,  \
                             ta1_z_z_xxxxz_0,  \
                             ta1_z_z_xxxxz_1,  \
                             ta1_z_z_xxxyz_0,  \
                             ta1_z_z_xxxyz_1,  \
                             ta1_z_z_xxxz_0,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxxzz_0,  \
                             ta1_z_z_xxxzz_1,  \
                             ta1_z_z_xxyyz_0,  \
                             ta1_z_z_xxyyz_1,  \
                             ta1_z_z_xxyz_0,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxyzz_0,  \
                             ta1_z_z_xxyzz_1,  \
                             ta1_z_z_xxzz_0,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xxzzz_0,  \
                             ta1_z_z_xxzzz_1,  \
                             ta1_z_z_xyyyz_0,  \
                             ta1_z_z_xyyyz_1,  \
                             ta1_z_z_xyyz_0,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyyzz_0,  \
                             ta1_z_z_xyyzz_1,  \
                             ta1_z_z_xyzz_0,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xyzzz_0,  \
                             ta1_z_z_xyzzz_1,  \
                             ta1_z_z_xzzz_0,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_xzzzz_0,  \
                             ta1_z_z_xzzzz_1,  \
                             ta1_z_z_yyyyz_0,  \
                             ta1_z_z_yyyyz_1,  \
                             ta1_z_z_yyyz_0,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyyzz_0,  \
                             ta1_z_z_yyyzz_1,  \
                             ta1_z_z_yyzz_0,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yyzzz_0,  \
                             ta1_z_z_yyzzz_1,  \
                             ta1_z_z_yzzz_0,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_yzzzz_0,  \
                             ta1_z_z_yzzzz_1,  \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta1_z_z_zzzzz_0,  \
                             ta1_z_z_zzzzz_1,  \
                             ta_y_xxxxy_1,     \
                             ta_y_xxxyy_1,     \
                             ta_y_xxyyy_1,     \
                             ta_y_xyyyy_1,     \
                             ta_y_yyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yz_xxxxx_0[i] = ta1_z_z_xxxxx_0[i] * pa_y[i] - ta1_z_z_xxxxx_1[i] * pc_y[i];

        ta1_z_yz_xxxxy_0[i] = ta_y_xxxxy_1[i] + ta1_z_y_xxxxy_0[i] * pa_z[i] - ta1_z_y_xxxxy_1[i] * pc_z[i];

        ta1_z_yz_xxxxz_0[i] = ta1_z_z_xxxxz_0[i] * pa_y[i] - ta1_z_z_xxxxz_1[i] * pc_y[i];

        ta1_z_yz_xxxyy_0[i] = ta_y_xxxyy_1[i] + ta1_z_y_xxxyy_0[i] * pa_z[i] - ta1_z_y_xxxyy_1[i] * pc_z[i];

        ta1_z_yz_xxxyz_0[i] = ta1_z_z_xxxz_0[i] * fe_0 - ta1_z_z_xxxz_1[i] * fe_0 + ta1_z_z_xxxyz_0[i] * pa_y[i] - ta1_z_z_xxxyz_1[i] * pc_y[i];

        ta1_z_yz_xxxzz_0[i] = ta1_z_z_xxxzz_0[i] * pa_y[i] - ta1_z_z_xxxzz_1[i] * pc_y[i];

        ta1_z_yz_xxyyy_0[i] = ta_y_xxyyy_1[i] + ta1_z_y_xxyyy_0[i] * pa_z[i] - ta1_z_y_xxyyy_1[i] * pc_z[i];

        ta1_z_yz_xxyyz_0[i] =
            2.0 * ta1_z_z_xxyz_0[i] * fe_0 - 2.0 * ta1_z_z_xxyz_1[i] * fe_0 + ta1_z_z_xxyyz_0[i] * pa_y[i] - ta1_z_z_xxyyz_1[i] * pc_y[i];

        ta1_z_yz_xxyzz_0[i] = ta1_z_z_xxzz_0[i] * fe_0 - ta1_z_z_xxzz_1[i] * fe_0 + ta1_z_z_xxyzz_0[i] * pa_y[i] - ta1_z_z_xxyzz_1[i] * pc_y[i];

        ta1_z_yz_xxzzz_0[i] = ta1_z_z_xxzzz_0[i] * pa_y[i] - ta1_z_z_xxzzz_1[i] * pc_y[i];

        ta1_z_yz_xyyyy_0[i] = ta_y_xyyyy_1[i] + ta1_z_y_xyyyy_0[i] * pa_z[i] - ta1_z_y_xyyyy_1[i] * pc_z[i];

        ta1_z_yz_xyyyz_0[i] =
            3.0 * ta1_z_z_xyyz_0[i] * fe_0 - 3.0 * ta1_z_z_xyyz_1[i] * fe_0 + ta1_z_z_xyyyz_0[i] * pa_y[i] - ta1_z_z_xyyyz_1[i] * pc_y[i];

        ta1_z_yz_xyyzz_0[i] =
            2.0 * ta1_z_z_xyzz_0[i] * fe_0 - 2.0 * ta1_z_z_xyzz_1[i] * fe_0 + ta1_z_z_xyyzz_0[i] * pa_y[i] - ta1_z_z_xyyzz_1[i] * pc_y[i];

        ta1_z_yz_xyzzz_0[i] = ta1_z_z_xzzz_0[i] * fe_0 - ta1_z_z_xzzz_1[i] * fe_0 + ta1_z_z_xyzzz_0[i] * pa_y[i] - ta1_z_z_xyzzz_1[i] * pc_y[i];

        ta1_z_yz_xzzzz_0[i] = ta1_z_z_xzzzz_0[i] * pa_y[i] - ta1_z_z_xzzzz_1[i] * pc_y[i];

        ta1_z_yz_yyyyy_0[i] = ta_y_yyyyy_1[i] + ta1_z_y_yyyyy_0[i] * pa_z[i] - ta1_z_y_yyyyy_1[i] * pc_z[i];

        ta1_z_yz_yyyyz_0[i] =
            4.0 * ta1_z_z_yyyz_0[i] * fe_0 - 4.0 * ta1_z_z_yyyz_1[i] * fe_0 + ta1_z_z_yyyyz_0[i] * pa_y[i] - ta1_z_z_yyyyz_1[i] * pc_y[i];

        ta1_z_yz_yyyzz_0[i] =
            3.0 * ta1_z_z_yyzz_0[i] * fe_0 - 3.0 * ta1_z_z_yyzz_1[i] * fe_0 + ta1_z_z_yyyzz_0[i] * pa_y[i] - ta1_z_z_yyyzz_1[i] * pc_y[i];

        ta1_z_yz_yyzzz_0[i] =
            2.0 * ta1_z_z_yzzz_0[i] * fe_0 - 2.0 * ta1_z_z_yzzz_1[i] * fe_0 + ta1_z_z_yyzzz_0[i] * pa_y[i] - ta1_z_z_yyzzz_1[i] * pc_y[i];

        ta1_z_yz_yzzzz_0[i] = ta1_z_z_zzzz_0[i] * fe_0 - ta1_z_z_zzzz_1[i] * fe_0 + ta1_z_z_yzzzz_0[i] * pa_y[i] - ta1_z_z_yzzzz_1[i] * pc_y[i];

        ta1_z_yz_zzzzz_0[i] = ta1_z_z_zzzzz_0[i] * pa_y[i] - ta1_z_z_zzzzz_1[i] * pc_y[i];
    }

    // Set up 357-378 components of targeted buffer : DH

    auto ta1_z_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 357);

    auto ta1_z_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 358);

    auto ta1_z_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 359);

    auto ta1_z_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 360);

    auto ta1_z_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 361);

    auto ta1_z_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 362);

    auto ta1_z_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 363);

    auto ta1_z_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 364);

    auto ta1_z_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 365);

    auto ta1_z_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 366);

    auto ta1_z_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 367);

    auto ta1_z_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 368);

    auto ta1_z_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 369);

    auto ta1_z_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 370);

    auto ta1_z_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 371);

    auto ta1_z_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 372);

    auto ta1_z_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 373);

    auto ta1_z_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 374);

    auto ta1_z_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 375);

    auto ta1_z_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 376);

    auto ta1_z_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 377);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxy_0,  \
                             ta1_z_0_xxxxy_1,  \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyz_0,  \
                             ta1_z_0_xxxyz_1,  \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyz_0,  \
                             ta1_z_0_xxyyz_1,  \
                             ta1_z_0_xxyzz_0,  \
                             ta1_z_0_xxyzz_1,  \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyz_0,  \
                             ta1_z_0_xyyyz_1,  \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyzzz_0,  \
                             ta1_z_0_xyzzz_1,  \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_z_xxxx_0,   \
                             ta1_z_z_xxxx_1,   \
                             ta1_z_z_xxxxx_0,  \
                             ta1_z_z_xxxxx_1,  \
                             ta1_z_z_xxxxy_0,  \
                             ta1_z_z_xxxxy_1,  \
                             ta1_z_z_xxxxz_0,  \
                             ta1_z_z_xxxxz_1,  \
                             ta1_z_z_xxxy_0,   \
                             ta1_z_z_xxxy_1,   \
                             ta1_z_z_xxxyy_0,  \
                             ta1_z_z_xxxyy_1,  \
                             ta1_z_z_xxxyz_0,  \
                             ta1_z_z_xxxyz_1,  \
                             ta1_z_z_xxxz_0,   \
                             ta1_z_z_xxxz_1,   \
                             ta1_z_z_xxxzz_0,  \
                             ta1_z_z_xxxzz_1,  \
                             ta1_z_z_xxyy_0,   \
                             ta1_z_z_xxyy_1,   \
                             ta1_z_z_xxyyy_0,  \
                             ta1_z_z_xxyyy_1,  \
                             ta1_z_z_xxyyz_0,  \
                             ta1_z_z_xxyyz_1,  \
                             ta1_z_z_xxyz_0,   \
                             ta1_z_z_xxyz_1,   \
                             ta1_z_z_xxyzz_0,  \
                             ta1_z_z_xxyzz_1,  \
                             ta1_z_z_xxzz_0,   \
                             ta1_z_z_xxzz_1,   \
                             ta1_z_z_xxzzz_0,  \
                             ta1_z_z_xxzzz_1,  \
                             ta1_z_z_xyyy_0,   \
                             ta1_z_z_xyyy_1,   \
                             ta1_z_z_xyyyy_0,  \
                             ta1_z_z_xyyyy_1,  \
                             ta1_z_z_xyyyz_0,  \
                             ta1_z_z_xyyyz_1,  \
                             ta1_z_z_xyyz_0,   \
                             ta1_z_z_xyyz_1,   \
                             ta1_z_z_xyyzz_0,  \
                             ta1_z_z_xyyzz_1,  \
                             ta1_z_z_xyzz_0,   \
                             ta1_z_z_xyzz_1,   \
                             ta1_z_z_xyzzz_0,  \
                             ta1_z_z_xyzzz_1,  \
                             ta1_z_z_xzzz_0,   \
                             ta1_z_z_xzzz_1,   \
                             ta1_z_z_xzzzz_0,  \
                             ta1_z_z_xzzzz_1,  \
                             ta1_z_z_yyyy_0,   \
                             ta1_z_z_yyyy_1,   \
                             ta1_z_z_yyyyy_0,  \
                             ta1_z_z_yyyyy_1,  \
                             ta1_z_z_yyyyz_0,  \
                             ta1_z_z_yyyyz_1,  \
                             ta1_z_z_yyyz_0,   \
                             ta1_z_z_yyyz_1,   \
                             ta1_z_z_yyyzz_0,  \
                             ta1_z_z_yyyzz_1,  \
                             ta1_z_z_yyzz_0,   \
                             ta1_z_z_yyzz_1,   \
                             ta1_z_z_yyzzz_0,  \
                             ta1_z_z_yyzzz_1,  \
                             ta1_z_z_yzzz_0,   \
                             ta1_z_z_yzzz_1,   \
                             ta1_z_z_yzzzz_0,  \
                             ta1_z_z_yzzzz_1,  \
                             ta1_z_z_zzzz_0,   \
                             ta1_z_z_zzzz_1,   \
                             ta1_z_z_zzzzz_0,  \
                             ta1_z_z_zzzzz_1,  \
                             ta1_z_zz_xxxxx_0, \
                             ta1_z_zz_xxxxy_0, \
                             ta1_z_zz_xxxxz_0, \
                             ta1_z_zz_xxxyy_0, \
                             ta1_z_zz_xxxyz_0, \
                             ta1_z_zz_xxxzz_0, \
                             ta1_z_zz_xxyyy_0, \
                             ta1_z_zz_xxyyz_0, \
                             ta1_z_zz_xxyzz_0, \
                             ta1_z_zz_xxzzz_0, \
                             ta1_z_zz_xyyyy_0, \
                             ta1_z_zz_xyyyz_0, \
                             ta1_z_zz_xyyzz_0, \
                             ta1_z_zz_xyzzz_0, \
                             ta1_z_zz_xzzzz_0, \
                             ta1_z_zz_yyyyy_0, \
                             ta1_z_zz_yyyyz_0, \
                             ta1_z_zz_yyyzz_0, \
                             ta1_z_zz_yyzzz_0, \
                             ta1_z_zz_yzzzz_0, \
                             ta1_z_zz_zzzzz_0, \
                             ta_z_xxxxx_1,     \
                             ta_z_xxxxy_1,     \
                             ta_z_xxxxz_1,     \
                             ta_z_xxxyy_1,     \
                             ta_z_xxxyz_1,     \
                             ta_z_xxxzz_1,     \
                             ta_z_xxyyy_1,     \
                             ta_z_xxyyz_1,     \
                             ta_z_xxyzz_1,     \
                             ta_z_xxzzz_1,     \
                             ta_z_xyyyy_1,     \
                             ta_z_xyyyz_1,     \
                             ta_z_xyyzz_1,     \
                             ta_z_xyzzz_1,     \
                             ta_z_xzzzz_1,     \
                             ta_z_yyyyy_1,     \
                             ta_z_yyyyz_1,     \
                             ta_z_yyyzz_1,     \
                             ta_z_yyzzz_1,     \
                             ta_z_yzzzz_1,     \
                             ta_z_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zz_xxxxx_0[i] =
            ta1_z_0_xxxxx_0[i] * fe_0 - ta1_z_0_xxxxx_1[i] * fe_0 + ta_z_xxxxx_1[i] + ta1_z_z_xxxxx_0[i] * pa_z[i] - ta1_z_z_xxxxx_1[i] * pc_z[i];

        ta1_z_zz_xxxxy_0[i] =
            ta1_z_0_xxxxy_0[i] * fe_0 - ta1_z_0_xxxxy_1[i] * fe_0 + ta_z_xxxxy_1[i] + ta1_z_z_xxxxy_0[i] * pa_z[i] - ta1_z_z_xxxxy_1[i] * pc_z[i];

        ta1_z_zz_xxxxz_0[i] = ta1_z_0_xxxxz_0[i] * fe_0 - ta1_z_0_xxxxz_1[i] * fe_0 + ta1_z_z_xxxx_0[i] * fe_0 - ta1_z_z_xxxx_1[i] * fe_0 +
                              ta_z_xxxxz_1[i] + ta1_z_z_xxxxz_0[i] * pa_z[i] - ta1_z_z_xxxxz_1[i] * pc_z[i];

        ta1_z_zz_xxxyy_0[i] =
            ta1_z_0_xxxyy_0[i] * fe_0 - ta1_z_0_xxxyy_1[i] * fe_0 + ta_z_xxxyy_1[i] + ta1_z_z_xxxyy_0[i] * pa_z[i] - ta1_z_z_xxxyy_1[i] * pc_z[i];

        ta1_z_zz_xxxyz_0[i] = ta1_z_0_xxxyz_0[i] * fe_0 - ta1_z_0_xxxyz_1[i] * fe_0 + ta1_z_z_xxxy_0[i] * fe_0 - ta1_z_z_xxxy_1[i] * fe_0 +
                              ta_z_xxxyz_1[i] + ta1_z_z_xxxyz_0[i] * pa_z[i] - ta1_z_z_xxxyz_1[i] * pc_z[i];

        ta1_z_zz_xxxzz_0[i] = ta1_z_0_xxxzz_0[i] * fe_0 - ta1_z_0_xxxzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxxz_0[i] * fe_0 -
                              2.0 * ta1_z_z_xxxz_1[i] * fe_0 + ta_z_xxxzz_1[i] + ta1_z_z_xxxzz_0[i] * pa_z[i] - ta1_z_z_xxxzz_1[i] * pc_z[i];

        ta1_z_zz_xxyyy_0[i] =
            ta1_z_0_xxyyy_0[i] * fe_0 - ta1_z_0_xxyyy_1[i] * fe_0 + ta_z_xxyyy_1[i] + ta1_z_z_xxyyy_0[i] * pa_z[i] - ta1_z_z_xxyyy_1[i] * pc_z[i];

        ta1_z_zz_xxyyz_0[i] = ta1_z_0_xxyyz_0[i] * fe_0 - ta1_z_0_xxyyz_1[i] * fe_0 + ta1_z_z_xxyy_0[i] * fe_0 - ta1_z_z_xxyy_1[i] * fe_0 +
                              ta_z_xxyyz_1[i] + ta1_z_z_xxyyz_0[i] * pa_z[i] - ta1_z_z_xxyyz_1[i] * pc_z[i];

        ta1_z_zz_xxyzz_0[i] = ta1_z_0_xxyzz_0[i] * fe_0 - ta1_z_0_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xxyz_0[i] * fe_0 -
                              2.0 * ta1_z_z_xxyz_1[i] * fe_0 + ta_z_xxyzz_1[i] + ta1_z_z_xxyzz_0[i] * pa_z[i] - ta1_z_z_xxyzz_1[i] * pc_z[i];

        ta1_z_zz_xxzzz_0[i] = ta1_z_0_xxzzz_0[i] * fe_0 - ta1_z_0_xxzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xxzz_0[i] * fe_0 -
                              3.0 * ta1_z_z_xxzz_1[i] * fe_0 + ta_z_xxzzz_1[i] + ta1_z_z_xxzzz_0[i] * pa_z[i] - ta1_z_z_xxzzz_1[i] * pc_z[i];

        ta1_z_zz_xyyyy_0[i] =
            ta1_z_0_xyyyy_0[i] * fe_0 - ta1_z_0_xyyyy_1[i] * fe_0 + ta_z_xyyyy_1[i] + ta1_z_z_xyyyy_0[i] * pa_z[i] - ta1_z_z_xyyyy_1[i] * pc_z[i];

        ta1_z_zz_xyyyz_0[i] = ta1_z_0_xyyyz_0[i] * fe_0 - ta1_z_0_xyyyz_1[i] * fe_0 + ta1_z_z_xyyy_0[i] * fe_0 - ta1_z_z_xyyy_1[i] * fe_0 +
                              ta_z_xyyyz_1[i] + ta1_z_z_xyyyz_0[i] * pa_z[i] - ta1_z_z_xyyyz_1[i] * pc_z[i];

        ta1_z_zz_xyyzz_0[i] = ta1_z_0_xyyzz_0[i] * fe_0 - ta1_z_0_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_z_xyyz_0[i] * fe_0 -
                              2.0 * ta1_z_z_xyyz_1[i] * fe_0 + ta_z_xyyzz_1[i] + ta1_z_z_xyyzz_0[i] * pa_z[i] - ta1_z_z_xyyzz_1[i] * pc_z[i];

        ta1_z_zz_xyzzz_0[i] = ta1_z_0_xyzzz_0[i] * fe_0 - ta1_z_0_xyzzz_1[i] * fe_0 + 3.0 * ta1_z_z_xyzz_0[i] * fe_0 -
                              3.0 * ta1_z_z_xyzz_1[i] * fe_0 + ta_z_xyzzz_1[i] + ta1_z_z_xyzzz_0[i] * pa_z[i] - ta1_z_z_xyzzz_1[i] * pc_z[i];

        ta1_z_zz_xzzzz_0[i] = ta1_z_0_xzzzz_0[i] * fe_0 - ta1_z_0_xzzzz_1[i] * fe_0 + 4.0 * ta1_z_z_xzzz_0[i] * fe_0 -
                              4.0 * ta1_z_z_xzzz_1[i] * fe_0 + ta_z_xzzzz_1[i] + ta1_z_z_xzzzz_0[i] * pa_z[i] - ta1_z_z_xzzzz_1[i] * pc_z[i];

        ta1_z_zz_yyyyy_0[i] =
            ta1_z_0_yyyyy_0[i] * fe_0 - ta1_z_0_yyyyy_1[i] * fe_0 + ta_z_yyyyy_1[i] + ta1_z_z_yyyyy_0[i] * pa_z[i] - ta1_z_z_yyyyy_1[i] * pc_z[i];

        ta1_z_zz_yyyyz_0[i] = ta1_z_0_yyyyz_0[i] * fe_0 - ta1_z_0_yyyyz_1[i] * fe_0 + ta1_z_z_yyyy_0[i] * fe_0 - ta1_z_z_yyyy_1[i] * fe_0 +
                              ta_z_yyyyz_1[i] + ta1_z_z_yyyyz_0[i] * pa_z[i] - ta1_z_z_yyyyz_1[i] * pc_z[i];

        ta1_z_zz_yyyzz_0[i] = ta1_z_0_yyyzz_0[i] * fe_0 - ta1_z_0_yyyzz_1[i] * fe_0 + 2.0 * ta1_z_z_yyyz_0[i] * fe_0 -
                              2.0 * ta1_z_z_yyyz_1[i] * fe_0 + ta_z_yyyzz_1[i] + ta1_z_z_yyyzz_0[i] * pa_z[i] - ta1_z_z_yyyzz_1[i] * pc_z[i];

        ta1_z_zz_yyzzz_0[i] = ta1_z_0_yyzzz_0[i] * fe_0 - ta1_z_0_yyzzz_1[i] * fe_0 + 3.0 * ta1_z_z_yyzz_0[i] * fe_0 -
                              3.0 * ta1_z_z_yyzz_1[i] * fe_0 + ta_z_yyzzz_1[i] + ta1_z_z_yyzzz_0[i] * pa_z[i] - ta1_z_z_yyzzz_1[i] * pc_z[i];

        ta1_z_zz_yzzzz_0[i] = ta1_z_0_yzzzz_0[i] * fe_0 - ta1_z_0_yzzzz_1[i] * fe_0 + 4.0 * ta1_z_z_yzzz_0[i] * fe_0 -
                              4.0 * ta1_z_z_yzzz_1[i] * fe_0 + ta_z_yzzzz_1[i] + ta1_z_z_yzzzz_0[i] * pa_z[i] - ta1_z_z_yzzzz_1[i] * pc_z[i];

        ta1_z_zz_zzzzz_0[i] = ta1_z_0_zzzzz_0[i] * fe_0 - ta1_z_0_zzzzz_1[i] * fe_0 + 5.0 * ta1_z_z_zzzz_0[i] * fe_0 -
                              5.0 * ta1_z_z_zzzz_1[i] * fe_0 + ta_z_zzzzz_1[i] + ta1_z_z_zzzzz_0[i] * pa_z[i] - ta1_z_z_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
