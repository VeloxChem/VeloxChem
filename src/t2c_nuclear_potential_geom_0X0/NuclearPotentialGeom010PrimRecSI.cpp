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

#include "NuclearPotentialGeom010PrimRecSI.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_si(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_si,
                                        const size_t              idx_npot_geom_010_0_sg,
                                        const size_t              idx_npot_geom_010_1_sg,
                                        const size_t              idx_npot_1_sh,
                                        const size_t              idx_npot_geom_010_0_sh,
                                        const size_t              idx_npot_geom_010_1_sh,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpb,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg);

    auto ta1_x_0_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 1);

    auto ta1_x_0_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 2);

    auto ta1_x_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 3);

    auto ta1_x_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 5);

    auto ta1_x_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 10);

    auto ta1_x_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 12);

    auto ta1_x_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 13);

    auto ta1_x_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 14);

    auto ta1_y_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 15);

    auto ta1_y_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 18);

    auto ta1_y_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 20);

    auto ta1_y_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 21);

    auto ta1_y_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 24);

    auto ta1_y_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 25);

    auto ta1_y_0_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 26);

    auto ta1_y_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 27);

    auto ta1_y_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 29);

    auto ta1_z_0_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_sg + 30);

    auto ta1_z_0_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 33);

    auto ta1_z_0_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 35);

    auto ta1_z_0_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 36);

    auto ta1_z_0_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 39);

    auto ta1_z_0_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_sg + 40);

    auto ta1_z_0_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 42);

    auto ta1_z_0_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 43);

    auto ta1_z_0_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_sg + 44);

    // Set up components of auxiliary buffer : SG

    auto ta1_x_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg);

    auto ta1_x_0_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 1);

    auto ta1_x_0_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 2);

    auto ta1_x_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 3);

    auto ta1_x_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 5);

    auto ta1_x_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 10);

    auto ta1_x_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 12);

    auto ta1_x_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 13);

    auto ta1_x_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 14);

    auto ta1_y_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 15);

    auto ta1_y_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 18);

    auto ta1_y_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 20);

    auto ta1_y_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 21);

    auto ta1_y_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 24);

    auto ta1_y_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 25);

    auto ta1_y_0_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 26);

    auto ta1_y_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 27);

    auto ta1_y_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 29);

    auto ta1_z_0_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_sg + 30);

    auto ta1_z_0_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 33);

    auto ta1_z_0_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 35);

    auto ta1_z_0_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 36);

    auto ta1_z_0_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 39);

    auto ta1_z_0_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_sg + 40);

    auto ta1_z_0_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 42);

    auto ta1_z_0_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 43);

    auto ta1_z_0_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_sg + 44);

    // Set up components of auxiliary buffer : SH

    auto ta_0_xxxxx_1 = pbuffer.data(idx_npot_1_sh);

    auto ta_0_xxxyy_1 = pbuffer.data(idx_npot_1_sh + 3);

    auto ta_0_xxxzz_1 = pbuffer.data(idx_npot_1_sh + 5);

    auto ta_0_xxyyy_1 = pbuffer.data(idx_npot_1_sh + 6);

    auto ta_0_xxzzz_1 = pbuffer.data(idx_npot_1_sh + 9);

    auto ta_0_yyyyy_1 = pbuffer.data(idx_npot_1_sh + 15);

    auto ta_0_yyyzz_1 = pbuffer.data(idx_npot_1_sh + 17);

    auto ta_0_yyzzz_1 = pbuffer.data(idx_npot_1_sh + 18);

    auto ta_0_zzzzz_1 = pbuffer.data(idx_npot_1_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto ta1_x_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh);

    auto ta1_x_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 1);

    auto ta1_x_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 2);

    auto ta1_x_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 3);

    auto ta1_x_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 5);

    auto ta1_x_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 6);

    auto ta1_x_0_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 8);

    auto ta1_x_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 9);

    auto ta1_x_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 10);

    auto ta1_x_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 14);

    auto ta1_x_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 15);

    auto ta1_x_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 17);

    auto ta1_x_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 18);

    auto ta1_x_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 19);

    auto ta1_x_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 20);

    auto ta1_y_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 21);

    auto ta1_y_0_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 22);

    auto ta1_y_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 24);

    auto ta1_y_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 26);

    auto ta1_y_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 27);

    auto ta1_y_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 30);

    auto ta1_y_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 31);

    auto ta1_y_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 33);

    auto ta1_y_0_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 35);

    auto ta1_y_0_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 36);

    auto ta1_y_0_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 37);

    auto ta1_y_0_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 38);

    auto ta1_y_0_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 39);

    auto ta1_y_0_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 40);

    auto ta1_y_0_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 41);

    auto ta1_z_0_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_sh + 42);

    auto ta1_z_0_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 44);

    auto ta1_z_0_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 45);

    auto ta1_z_0_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 47);

    auto ta1_z_0_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 48);

    auto ta1_z_0_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 51);

    auto ta1_z_0_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_sh + 52);

    auto ta1_z_0_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_sh + 54);

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

    auto ta1_x_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 5);

    auto ta1_x_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 6);

    auto ta1_x_0_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 8);

    auto ta1_x_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 9);

    auto ta1_x_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 10);

    auto ta1_x_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 14);

    auto ta1_x_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 15);

    auto ta1_x_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 17);

    auto ta1_x_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 18);

    auto ta1_x_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 19);

    auto ta1_x_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 20);

    auto ta1_y_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 21);

    auto ta1_y_0_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 22);

    auto ta1_y_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 24);

    auto ta1_y_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 26);

    auto ta1_y_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 27);

    auto ta1_y_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 30);

    auto ta1_y_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 31);

    auto ta1_y_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 33);

    auto ta1_y_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 35);

    auto ta1_y_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 36);

    auto ta1_y_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 37);

    auto ta1_y_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 38);

    auto ta1_y_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 39);

    auto ta1_y_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 40);

    auto ta1_y_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 41);

    auto ta1_z_0_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_sh + 42);

    auto ta1_z_0_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 44);

    auto ta1_z_0_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 45);

    auto ta1_z_0_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 47);

    auto ta1_z_0_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 48);

    auto ta1_z_0_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 51);

    auto ta1_z_0_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 52);

    auto ta1_z_0_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 54);

    auto ta1_z_0_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 56);

    auto ta1_z_0_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_sh + 57);

    auto ta1_z_0_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 58);

    auto ta1_z_0_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 59);

    auto ta1_z_0_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 60);

    auto ta1_z_0_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 61);

    auto ta1_z_0_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_sh + 62);

    // Set up components of targeted buffer : SI

    auto ta1_x_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si);

    auto ta1_x_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 1);

    auto ta1_x_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 2);

    auto ta1_x_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 3);

    auto ta1_x_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 4);

    auto ta1_x_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 5);

    auto ta1_x_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 6);

    auto ta1_x_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 7);

    auto ta1_x_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 8);

    auto ta1_x_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 9);

    auto ta1_x_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 10);

    auto ta1_x_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 11);

    auto ta1_x_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 12);

    auto ta1_x_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 13);

    auto ta1_x_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 14);

    auto ta1_x_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 15);

    auto ta1_x_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 16);

    auto ta1_x_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 17);

    auto ta1_x_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 18);

    auto ta1_x_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 19);

    auto ta1_x_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 20);

    auto ta1_x_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 21);

    auto ta1_x_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 22);

    auto ta1_x_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 23);

    auto ta1_x_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 24);

    auto ta1_x_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 25);

    auto ta1_x_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 26);

    auto ta1_x_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 27);

    auto ta1_y_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si + 28);

    auto ta1_y_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 29);

    auto ta1_y_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 30);

    auto ta1_y_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 31);

    auto ta1_y_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 32);

    auto ta1_y_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 33);

    auto ta1_y_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 34);

    auto ta1_y_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 35);

    auto ta1_y_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 36);

    auto ta1_y_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 37);

    auto ta1_y_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 38);

    auto ta1_y_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 39);

    auto ta1_y_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 40);

    auto ta1_y_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 41);

    auto ta1_y_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 42);

    auto ta1_y_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 43);

    auto ta1_y_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 44);

    auto ta1_y_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 45);

    auto ta1_y_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 46);

    auto ta1_y_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 47);

    auto ta1_y_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 48);

    auto ta1_y_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 49);

    auto ta1_y_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 50);

    auto ta1_y_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 51);

    auto ta1_y_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 52);

    auto ta1_y_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 53);

    auto ta1_y_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 54);

    auto ta1_y_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 55);

    auto ta1_z_0_xxxxxx_0 = pbuffer.data(idx_npot_geom_010_0_si + 56);

    auto ta1_z_0_xxxxxy_0 = pbuffer.data(idx_npot_geom_010_0_si + 57);

    auto ta1_z_0_xxxxxz_0 = pbuffer.data(idx_npot_geom_010_0_si + 58);

    auto ta1_z_0_xxxxyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 59);

    auto ta1_z_0_xxxxyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 60);

    auto ta1_z_0_xxxxzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 61);

    auto ta1_z_0_xxxyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 62);

    auto ta1_z_0_xxxyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 63);

    auto ta1_z_0_xxxyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 64);

    auto ta1_z_0_xxxzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 65);

    auto ta1_z_0_xxyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 66);

    auto ta1_z_0_xxyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 67);

    auto ta1_z_0_xxyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 68);

    auto ta1_z_0_xxyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 69);

    auto ta1_z_0_xxzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 70);

    auto ta1_z_0_xyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 71);

    auto ta1_z_0_xyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 72);

    auto ta1_z_0_xyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 73);

    auto ta1_z_0_xyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 74);

    auto ta1_z_0_xyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 75);

    auto ta1_z_0_xzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 76);

    auto ta1_z_0_yyyyyy_0 = pbuffer.data(idx_npot_geom_010_0_si + 77);

    auto ta1_z_0_yyyyyz_0 = pbuffer.data(idx_npot_geom_010_0_si + 78);

    auto ta1_z_0_yyyyzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 79);

    auto ta1_z_0_yyyzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 80);

    auto ta1_z_0_yyzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 81);

    auto ta1_z_0_yzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 82);

    auto ta1_z_0_zzzzzz_0 = pbuffer.data(idx_npot_geom_010_0_si + 83);

#pragma omp simd aligned(pb_x,                 \
                             pb_y,             \
                             pb_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_0_xxxx_0,   \
                             ta1_x_0_xxxx_1,   \
                             ta1_x_0_xxxxx_0,  \
                             ta1_x_0_xxxxx_1,  \
                             ta1_x_0_xxxxxx_0, \
                             ta1_x_0_xxxxxy_0, \
                             ta1_x_0_xxxxxz_0, \
                             ta1_x_0_xxxxy_0,  \
                             ta1_x_0_xxxxy_1,  \
                             ta1_x_0_xxxxyy_0, \
                             ta1_x_0_xxxxyz_0, \
                             ta1_x_0_xxxxz_0,  \
                             ta1_x_0_xxxxz_1,  \
                             ta1_x_0_xxxxzz_0, \
                             ta1_x_0_xxxy_0,   \
                             ta1_x_0_xxxy_1,   \
                             ta1_x_0_xxxyy_0,  \
                             ta1_x_0_xxxyy_1,  \
                             ta1_x_0_xxxyyy_0, \
                             ta1_x_0_xxxyyz_0, \
                             ta1_x_0_xxxyzz_0, \
                             ta1_x_0_xxxz_0,   \
                             ta1_x_0_xxxz_1,   \
                             ta1_x_0_xxxzz_0,  \
                             ta1_x_0_xxxzz_1,  \
                             ta1_x_0_xxxzzz_0, \
                             ta1_x_0_xxyy_0,   \
                             ta1_x_0_xxyy_1,   \
                             ta1_x_0_xxyyy_0,  \
                             ta1_x_0_xxyyy_1,  \
                             ta1_x_0_xxyyyy_0, \
                             ta1_x_0_xxyyyz_0, \
                             ta1_x_0_xxyyzz_0, \
                             ta1_x_0_xxyzz_0,  \
                             ta1_x_0_xxyzz_1,  \
                             ta1_x_0_xxyzzz_0, \
                             ta1_x_0_xxzz_0,   \
                             ta1_x_0_xxzz_1,   \
                             ta1_x_0_xxzzz_0,  \
                             ta1_x_0_xxzzz_1,  \
                             ta1_x_0_xxzzzz_0, \
                             ta1_x_0_xyyyy_0,  \
                             ta1_x_0_xyyyy_1,  \
                             ta1_x_0_xyyyyy_0, \
                             ta1_x_0_xyyyyz_0, \
                             ta1_x_0_xyyyzz_0, \
                             ta1_x_0_xyyzzz_0, \
                             ta1_x_0_xyzzzz_0, \
                             ta1_x_0_xzzzz_0,  \
                             ta1_x_0_xzzzz_1,  \
                             ta1_x_0_xzzzzz_0, \
                             ta1_x_0_yyyy_0,   \
                             ta1_x_0_yyyy_1,   \
                             ta1_x_0_yyyyy_0,  \
                             ta1_x_0_yyyyy_1,  \
                             ta1_x_0_yyyyyy_0, \
                             ta1_x_0_yyyyyz_0, \
                             ta1_x_0_yyyyzz_0, \
                             ta1_x_0_yyyzz_0,  \
                             ta1_x_0_yyyzz_1,  \
                             ta1_x_0_yyyzzz_0, \
                             ta1_x_0_yyzz_0,   \
                             ta1_x_0_yyzz_1,   \
                             ta1_x_0_yyzzz_0,  \
                             ta1_x_0_yyzzz_1,  \
                             ta1_x_0_yyzzzz_0, \
                             ta1_x_0_yzzz_0,   \
                             ta1_x_0_yzzz_1,   \
                             ta1_x_0_yzzzz_0,  \
                             ta1_x_0_yzzzz_1,  \
                             ta1_x_0_yzzzzz_0, \
                             ta1_x_0_zzzz_0,   \
                             ta1_x_0_zzzz_1,   \
                             ta1_x_0_zzzzz_0,  \
                             ta1_x_0_zzzzz_1,  \
                             ta1_x_0_zzzzzz_0, \
                             ta1_y_0_xxxx_0,   \
                             ta1_y_0_xxxx_1,   \
                             ta1_y_0_xxxxx_0,  \
                             ta1_y_0_xxxxx_1,  \
                             ta1_y_0_xxxxxx_0, \
                             ta1_y_0_xxxxxy_0, \
                             ta1_y_0_xxxxxz_0, \
                             ta1_y_0_xxxxy_0,  \
                             ta1_y_0_xxxxy_1,  \
                             ta1_y_0_xxxxyy_0, \
                             ta1_y_0_xxxxyz_0, \
                             ta1_y_0_xxxxzz_0, \
                             ta1_y_0_xxxyy_0,  \
                             ta1_y_0_xxxyy_1,  \
                             ta1_y_0_xxxyyy_0, \
                             ta1_y_0_xxxyyz_0, \
                             ta1_y_0_xxxyzz_0, \
                             ta1_y_0_xxxzz_0,  \
                             ta1_y_0_xxxzz_1,  \
                             ta1_y_0_xxxzzz_0, \
                             ta1_y_0_xxyy_0,   \
                             ta1_y_0_xxyy_1,   \
                             ta1_y_0_xxyyy_0,  \
                             ta1_y_0_xxyyy_1,  \
                             ta1_y_0_xxyyyy_0, \
                             ta1_y_0_xxyyyz_0, \
                             ta1_y_0_xxyyzz_0, \
                             ta1_y_0_xxyzzz_0, \
                             ta1_y_0_xxzz_0,   \
                             ta1_y_0_xxzz_1,   \
                             ta1_y_0_xxzzz_0,  \
                             ta1_y_0_xxzzz_1,  \
                             ta1_y_0_xxzzzz_0, \
                             ta1_y_0_xyyy_0,   \
                             ta1_y_0_xyyy_1,   \
                             ta1_y_0_xyyyy_0,  \
                             ta1_y_0_xyyyy_1,  \
                             ta1_y_0_xyyyyy_0, \
                             ta1_y_0_xyyyyz_0, \
                             ta1_y_0_xyyyzz_0, \
                             ta1_y_0_xyyzz_0,  \
                             ta1_y_0_xyyzz_1,  \
                             ta1_y_0_xyyzzz_0, \
                             ta1_y_0_xyzzzz_0, \
                             ta1_y_0_xzzz_0,   \
                             ta1_y_0_xzzz_1,   \
                             ta1_y_0_xzzzz_0,  \
                             ta1_y_0_xzzzz_1,  \
                             ta1_y_0_xzzzzz_0, \
                             ta1_y_0_yyyy_0,   \
                             ta1_y_0_yyyy_1,   \
                             ta1_y_0_yyyyy_0,  \
                             ta1_y_0_yyyyy_1,  \
                             ta1_y_0_yyyyyy_0, \
                             ta1_y_0_yyyyyz_0, \
                             ta1_y_0_yyyyz_0,  \
                             ta1_y_0_yyyyz_1,  \
                             ta1_y_0_yyyyzz_0, \
                             ta1_y_0_yyyz_0,   \
                             ta1_y_0_yyyz_1,   \
                             ta1_y_0_yyyzz_0,  \
                             ta1_y_0_yyyzz_1,  \
                             ta1_y_0_yyyzzz_0, \
                             ta1_y_0_yyzz_0,   \
                             ta1_y_0_yyzz_1,   \
                             ta1_y_0_yyzzz_0,  \
                             ta1_y_0_yyzzz_1,  \
                             ta1_y_0_yyzzzz_0, \
                             ta1_y_0_yzzzz_0,  \
                             ta1_y_0_yzzzz_1,  \
                             ta1_y_0_yzzzzz_0, \
                             ta1_y_0_zzzz_0,   \
                             ta1_y_0_zzzz_1,   \
                             ta1_y_0_zzzzz_0,  \
                             ta1_y_0_zzzzz_1,  \
                             ta1_y_0_zzzzzz_0, \
                             ta1_z_0_xxxx_0,   \
                             ta1_z_0_xxxx_1,   \
                             ta1_z_0_xxxxx_0,  \
                             ta1_z_0_xxxxx_1,  \
                             ta1_z_0_xxxxxx_0, \
                             ta1_z_0_xxxxxy_0, \
                             ta1_z_0_xxxxxz_0, \
                             ta1_z_0_xxxxyy_0, \
                             ta1_z_0_xxxxyz_0, \
                             ta1_z_0_xxxxz_0,  \
                             ta1_z_0_xxxxz_1,  \
                             ta1_z_0_xxxxzz_0, \
                             ta1_z_0_xxxyy_0,  \
                             ta1_z_0_xxxyy_1,  \
                             ta1_z_0_xxxyyy_0, \
                             ta1_z_0_xxxyyz_0, \
                             ta1_z_0_xxxyzz_0, \
                             ta1_z_0_xxxzz_0,  \
                             ta1_z_0_xxxzz_1,  \
                             ta1_z_0_xxxzzz_0, \
                             ta1_z_0_xxyy_0,   \
                             ta1_z_0_xxyy_1,   \
                             ta1_z_0_xxyyy_0,  \
                             ta1_z_0_xxyyy_1,  \
                             ta1_z_0_xxyyyy_0, \
                             ta1_z_0_xxyyyz_0, \
                             ta1_z_0_xxyyzz_0, \
                             ta1_z_0_xxyzzz_0, \
                             ta1_z_0_xxzz_0,   \
                             ta1_z_0_xxzz_1,   \
                             ta1_z_0_xxzzz_0,  \
                             ta1_z_0_xxzzz_1,  \
                             ta1_z_0_xxzzzz_0, \
                             ta1_z_0_xyyy_0,   \
                             ta1_z_0_xyyy_1,   \
                             ta1_z_0_xyyyy_0,  \
                             ta1_z_0_xyyyy_1,  \
                             ta1_z_0_xyyyyy_0, \
                             ta1_z_0_xyyyyz_0, \
                             ta1_z_0_xyyyzz_0, \
                             ta1_z_0_xyyzz_0,  \
                             ta1_z_0_xyyzz_1,  \
                             ta1_z_0_xyyzzz_0, \
                             ta1_z_0_xyzzzz_0, \
                             ta1_z_0_xzzz_0,   \
                             ta1_z_0_xzzz_1,   \
                             ta1_z_0_xzzzz_0,  \
                             ta1_z_0_xzzzz_1,  \
                             ta1_z_0_xzzzzz_0, \
                             ta1_z_0_yyyy_0,   \
                             ta1_z_0_yyyy_1,   \
                             ta1_z_0_yyyyy_0,  \
                             ta1_z_0_yyyyy_1,  \
                             ta1_z_0_yyyyyy_0, \
                             ta1_z_0_yyyyyz_0, \
                             ta1_z_0_yyyyz_0,  \
                             ta1_z_0_yyyyz_1,  \
                             ta1_z_0_yyyyzz_0, \
                             ta1_z_0_yyyzz_0,  \
                             ta1_z_0_yyyzz_1,  \
                             ta1_z_0_yyyzzz_0, \
                             ta1_z_0_yyzz_0,   \
                             ta1_z_0_yyzz_1,   \
                             ta1_z_0_yyzzz_0,  \
                             ta1_z_0_yyzzz_1,  \
                             ta1_z_0_yyzzzz_0, \
                             ta1_z_0_yzzz_0,   \
                             ta1_z_0_yzzz_1,   \
                             ta1_z_0_yzzzz_0,  \
                             ta1_z_0_yzzzz_1,  \
                             ta1_z_0_yzzzzz_0, \
                             ta1_z_0_zzzz_0,   \
                             ta1_z_0_zzzz_1,   \
                             ta1_z_0_zzzzz_0,  \
                             ta1_z_0_zzzzz_1,  \
                             ta1_z_0_zzzzzz_0, \
                             ta_0_xxxxx_1,     \
                             ta_0_xxxyy_1,     \
                             ta_0_xxxzz_1,     \
                             ta_0_xxyyy_1,     \
                             ta_0_xxzzz_1,     \
                             ta_0_yyyyy_1,     \
                             ta_0_yyyzz_1,     \
                             ta_0_yyzzz_1,     \
                             ta_0_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_0_xxxxxx_0[i] = 5.0 * ta1_x_0_xxxx_0[i] * fe_0 - 5.0 * ta1_x_0_xxxx_1[i] * fe_0 + ta_0_xxxxx_1[i] + ta1_x_0_xxxxx_0[i] * pb_x[i] -
                              ta1_x_0_xxxxx_1[i] * pc_x[i];

        ta1_x_0_xxxxxy_0[i] = ta1_x_0_xxxxx_0[i] * pb_y[i] - ta1_x_0_xxxxx_1[i] * pc_y[i];

        ta1_x_0_xxxxxz_0[i] = ta1_x_0_xxxxx_0[i] * pb_z[i] - ta1_x_0_xxxxx_1[i] * pc_z[i];

        ta1_x_0_xxxxyy_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_0_xxxxy_0[i] * pb_y[i] - ta1_x_0_xxxxy_1[i] * pc_y[i];

        ta1_x_0_xxxxyz_0[i] = ta1_x_0_xxxxz_0[i] * pb_y[i] - ta1_x_0_xxxxz_1[i] * pc_y[i];

        ta1_x_0_xxxxzz_0[i] = ta1_x_0_xxxx_0[i] * fe_0 - ta1_x_0_xxxx_1[i] * fe_0 + ta1_x_0_xxxxz_0[i] * pb_z[i] - ta1_x_0_xxxxz_1[i] * pc_z[i];

        ta1_x_0_xxxyyy_0[i] =
            2.0 * ta1_x_0_xxxy_0[i] * fe_0 - 2.0 * ta1_x_0_xxxy_1[i] * fe_0 + ta1_x_0_xxxyy_0[i] * pb_y[i] - ta1_x_0_xxxyy_1[i] * pc_y[i];

        ta1_x_0_xxxyyz_0[i] = ta1_x_0_xxxyy_0[i] * pb_z[i] - ta1_x_0_xxxyy_1[i] * pc_z[i];

        ta1_x_0_xxxyzz_0[i] = ta1_x_0_xxxzz_0[i] * pb_y[i] - ta1_x_0_xxxzz_1[i] * pc_y[i];

        ta1_x_0_xxxzzz_0[i] =
            2.0 * ta1_x_0_xxxz_0[i] * fe_0 - 2.0 * ta1_x_0_xxxz_1[i] * fe_0 + ta1_x_0_xxxzz_0[i] * pb_z[i] - ta1_x_0_xxxzz_1[i] * pc_z[i];

        ta1_x_0_xxyyyy_0[i] =
            3.0 * ta1_x_0_xxyy_0[i] * fe_0 - 3.0 * ta1_x_0_xxyy_1[i] * fe_0 + ta1_x_0_xxyyy_0[i] * pb_y[i] - ta1_x_0_xxyyy_1[i] * pc_y[i];

        ta1_x_0_xxyyyz_0[i] = ta1_x_0_xxyyy_0[i] * pb_z[i] - ta1_x_0_xxyyy_1[i] * pc_z[i];

        ta1_x_0_xxyyzz_0[i] = ta1_x_0_xxzz_0[i] * fe_0 - ta1_x_0_xxzz_1[i] * fe_0 + ta1_x_0_xxyzz_0[i] * pb_y[i] - ta1_x_0_xxyzz_1[i] * pc_y[i];

        ta1_x_0_xxyzzz_0[i] = ta1_x_0_xxzzz_0[i] * pb_y[i] - ta1_x_0_xxzzz_1[i] * pc_y[i];

        ta1_x_0_xxzzzz_0[i] =
            3.0 * ta1_x_0_xxzz_0[i] * fe_0 - 3.0 * ta1_x_0_xxzz_1[i] * fe_0 + ta1_x_0_xxzzz_0[i] * pb_z[i] - ta1_x_0_xxzzz_1[i] * pc_z[i];

        ta1_x_0_xyyyyy_0[i] = ta_0_yyyyy_1[i] + ta1_x_0_yyyyy_0[i] * pb_x[i] - ta1_x_0_yyyyy_1[i] * pc_x[i];

        ta1_x_0_xyyyyz_0[i] = ta1_x_0_xyyyy_0[i] * pb_z[i] - ta1_x_0_xyyyy_1[i] * pc_z[i];

        ta1_x_0_xyyyzz_0[i] = ta_0_yyyzz_1[i] + ta1_x_0_yyyzz_0[i] * pb_x[i] - ta1_x_0_yyyzz_1[i] * pc_x[i];

        ta1_x_0_xyyzzz_0[i] = ta_0_yyzzz_1[i] + ta1_x_0_yyzzz_0[i] * pb_x[i] - ta1_x_0_yyzzz_1[i] * pc_x[i];

        ta1_x_0_xyzzzz_0[i] = ta1_x_0_xzzzz_0[i] * pb_y[i] - ta1_x_0_xzzzz_1[i] * pc_y[i];

        ta1_x_0_xzzzzz_0[i] = ta_0_zzzzz_1[i] + ta1_x_0_zzzzz_0[i] * pb_x[i] - ta1_x_0_zzzzz_1[i] * pc_x[i];

        ta1_x_0_yyyyyy_0[i] =
            5.0 * ta1_x_0_yyyy_0[i] * fe_0 - 5.0 * ta1_x_0_yyyy_1[i] * fe_0 + ta1_x_0_yyyyy_0[i] * pb_y[i] - ta1_x_0_yyyyy_1[i] * pc_y[i];

        ta1_x_0_yyyyyz_0[i] = ta1_x_0_yyyyy_0[i] * pb_z[i] - ta1_x_0_yyyyy_1[i] * pc_z[i];

        ta1_x_0_yyyyzz_0[i] =
            3.0 * ta1_x_0_yyzz_0[i] * fe_0 - 3.0 * ta1_x_0_yyzz_1[i] * fe_0 + ta1_x_0_yyyzz_0[i] * pb_y[i] - ta1_x_0_yyyzz_1[i] * pc_y[i];

        ta1_x_0_yyyzzz_0[i] =
            2.0 * ta1_x_0_yzzz_0[i] * fe_0 - 2.0 * ta1_x_0_yzzz_1[i] * fe_0 + ta1_x_0_yyzzz_0[i] * pb_y[i] - ta1_x_0_yyzzz_1[i] * pc_y[i];

        ta1_x_0_yyzzzz_0[i] = ta1_x_0_zzzz_0[i] * fe_0 - ta1_x_0_zzzz_1[i] * fe_0 + ta1_x_0_yzzzz_0[i] * pb_y[i] - ta1_x_0_yzzzz_1[i] * pc_y[i];

        ta1_x_0_yzzzzz_0[i] = ta1_x_0_zzzzz_0[i] * pb_y[i] - ta1_x_0_zzzzz_1[i] * pc_y[i];

        ta1_x_0_zzzzzz_0[i] =
            5.0 * ta1_x_0_zzzz_0[i] * fe_0 - 5.0 * ta1_x_0_zzzz_1[i] * fe_0 + ta1_x_0_zzzzz_0[i] * pb_z[i] - ta1_x_0_zzzzz_1[i] * pc_z[i];

        ta1_y_0_xxxxxx_0[i] =
            5.0 * ta1_y_0_xxxx_0[i] * fe_0 - 5.0 * ta1_y_0_xxxx_1[i] * fe_0 + ta1_y_0_xxxxx_0[i] * pb_x[i] - ta1_y_0_xxxxx_1[i] * pc_x[i];

        ta1_y_0_xxxxxy_0[i] = ta_0_xxxxx_1[i] + ta1_y_0_xxxxx_0[i] * pb_y[i] - ta1_y_0_xxxxx_1[i] * pc_y[i];

        ta1_y_0_xxxxxz_0[i] = ta1_y_0_xxxxx_0[i] * pb_z[i] - ta1_y_0_xxxxx_1[i] * pc_z[i];

        ta1_y_0_xxxxyy_0[i] =
            3.0 * ta1_y_0_xxyy_0[i] * fe_0 - 3.0 * ta1_y_0_xxyy_1[i] * fe_0 + ta1_y_0_xxxyy_0[i] * pb_x[i] - ta1_y_0_xxxyy_1[i] * pc_x[i];

        ta1_y_0_xxxxyz_0[i] = ta1_y_0_xxxxy_0[i] * pb_z[i] - ta1_y_0_xxxxy_1[i] * pc_z[i];

        ta1_y_0_xxxxzz_0[i] =
            3.0 * ta1_y_0_xxzz_0[i] * fe_0 - 3.0 * ta1_y_0_xxzz_1[i] * fe_0 + ta1_y_0_xxxzz_0[i] * pb_x[i] - ta1_y_0_xxxzz_1[i] * pc_x[i];

        ta1_y_0_xxxyyy_0[i] =
            2.0 * ta1_y_0_xyyy_0[i] * fe_0 - 2.0 * ta1_y_0_xyyy_1[i] * fe_0 + ta1_y_0_xxyyy_0[i] * pb_x[i] - ta1_y_0_xxyyy_1[i] * pc_x[i];

        ta1_y_0_xxxyyz_0[i] = ta1_y_0_xxxyy_0[i] * pb_z[i] - ta1_y_0_xxxyy_1[i] * pc_z[i];

        ta1_y_0_xxxyzz_0[i] = ta_0_xxxzz_1[i] + ta1_y_0_xxxzz_0[i] * pb_y[i] - ta1_y_0_xxxzz_1[i] * pc_y[i];

        ta1_y_0_xxxzzz_0[i] =
            2.0 * ta1_y_0_xzzz_0[i] * fe_0 - 2.0 * ta1_y_0_xzzz_1[i] * fe_0 + ta1_y_0_xxzzz_0[i] * pb_x[i] - ta1_y_0_xxzzz_1[i] * pc_x[i];

        ta1_y_0_xxyyyy_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_0_xyyyy_0[i] * pb_x[i] - ta1_y_0_xyyyy_1[i] * pc_x[i];

        ta1_y_0_xxyyyz_0[i] = ta1_y_0_xxyyy_0[i] * pb_z[i] - ta1_y_0_xxyyy_1[i] * pc_z[i];

        ta1_y_0_xxyyzz_0[i] = ta1_y_0_yyzz_0[i] * fe_0 - ta1_y_0_yyzz_1[i] * fe_0 + ta1_y_0_xyyzz_0[i] * pb_x[i] - ta1_y_0_xyyzz_1[i] * pc_x[i];

        ta1_y_0_xxyzzz_0[i] = ta_0_xxzzz_1[i] + ta1_y_0_xxzzz_0[i] * pb_y[i] - ta1_y_0_xxzzz_1[i] * pc_y[i];

        ta1_y_0_xxzzzz_0[i] = ta1_y_0_zzzz_0[i] * fe_0 - ta1_y_0_zzzz_1[i] * fe_0 + ta1_y_0_xzzzz_0[i] * pb_x[i] - ta1_y_0_xzzzz_1[i] * pc_x[i];

        ta1_y_0_xyyyyy_0[i] = ta1_y_0_yyyyy_0[i] * pb_x[i] - ta1_y_0_yyyyy_1[i] * pc_x[i];

        ta1_y_0_xyyyyz_0[i] = ta1_y_0_yyyyz_0[i] * pb_x[i] - ta1_y_0_yyyyz_1[i] * pc_x[i];

        ta1_y_0_xyyyzz_0[i] = ta1_y_0_yyyzz_0[i] * pb_x[i] - ta1_y_0_yyyzz_1[i] * pc_x[i];

        ta1_y_0_xyyzzz_0[i] = ta1_y_0_yyzzz_0[i] * pb_x[i] - ta1_y_0_yyzzz_1[i] * pc_x[i];

        ta1_y_0_xyzzzz_0[i] = ta1_y_0_yzzzz_0[i] * pb_x[i] - ta1_y_0_yzzzz_1[i] * pc_x[i];

        ta1_y_0_xzzzzz_0[i] = ta1_y_0_zzzzz_0[i] * pb_x[i] - ta1_y_0_zzzzz_1[i] * pc_x[i];

        ta1_y_0_yyyyyy_0[i] = 5.0 * ta1_y_0_yyyy_0[i] * fe_0 - 5.0 * ta1_y_0_yyyy_1[i] * fe_0 + ta_0_yyyyy_1[i] + ta1_y_0_yyyyy_0[i] * pb_y[i] -
                              ta1_y_0_yyyyy_1[i] * pc_y[i];

        ta1_y_0_yyyyyz_0[i] = ta1_y_0_yyyyy_0[i] * pb_z[i] - ta1_y_0_yyyyy_1[i] * pc_z[i];

        ta1_y_0_yyyyzz_0[i] = ta1_y_0_yyyy_0[i] * fe_0 - ta1_y_0_yyyy_1[i] * fe_0 + ta1_y_0_yyyyz_0[i] * pb_z[i] - ta1_y_0_yyyyz_1[i] * pc_z[i];

        ta1_y_0_yyyzzz_0[i] =
            2.0 * ta1_y_0_yyyz_0[i] * fe_0 - 2.0 * ta1_y_0_yyyz_1[i] * fe_0 + ta1_y_0_yyyzz_0[i] * pb_z[i] - ta1_y_0_yyyzz_1[i] * pc_z[i];

        ta1_y_0_yyzzzz_0[i] =
            3.0 * ta1_y_0_yyzz_0[i] * fe_0 - 3.0 * ta1_y_0_yyzz_1[i] * fe_0 + ta1_y_0_yyzzz_0[i] * pb_z[i] - ta1_y_0_yyzzz_1[i] * pc_z[i];

        ta1_y_0_yzzzzz_0[i] = ta_0_zzzzz_1[i] + ta1_y_0_zzzzz_0[i] * pb_y[i] - ta1_y_0_zzzzz_1[i] * pc_y[i];

        ta1_y_0_zzzzzz_0[i] =
            5.0 * ta1_y_0_zzzz_0[i] * fe_0 - 5.0 * ta1_y_0_zzzz_1[i] * fe_0 + ta1_y_0_zzzzz_0[i] * pb_z[i] - ta1_y_0_zzzzz_1[i] * pc_z[i];

        ta1_z_0_xxxxxx_0[i] =
            5.0 * ta1_z_0_xxxx_0[i] * fe_0 - 5.0 * ta1_z_0_xxxx_1[i] * fe_0 + ta1_z_0_xxxxx_0[i] * pb_x[i] - ta1_z_0_xxxxx_1[i] * pc_x[i];

        ta1_z_0_xxxxxy_0[i] = ta1_z_0_xxxxx_0[i] * pb_y[i] - ta1_z_0_xxxxx_1[i] * pc_y[i];

        ta1_z_0_xxxxxz_0[i] = ta_0_xxxxx_1[i] + ta1_z_0_xxxxx_0[i] * pb_z[i] - ta1_z_0_xxxxx_1[i] * pc_z[i];

        ta1_z_0_xxxxyy_0[i] =
            3.0 * ta1_z_0_xxyy_0[i] * fe_0 - 3.0 * ta1_z_0_xxyy_1[i] * fe_0 + ta1_z_0_xxxyy_0[i] * pb_x[i] - ta1_z_0_xxxyy_1[i] * pc_x[i];

        ta1_z_0_xxxxyz_0[i] = ta1_z_0_xxxxz_0[i] * pb_y[i] - ta1_z_0_xxxxz_1[i] * pc_y[i];

        ta1_z_0_xxxxzz_0[i] =
            3.0 * ta1_z_0_xxzz_0[i] * fe_0 - 3.0 * ta1_z_0_xxzz_1[i] * fe_0 + ta1_z_0_xxxzz_0[i] * pb_x[i] - ta1_z_0_xxxzz_1[i] * pc_x[i];

        ta1_z_0_xxxyyy_0[i] =
            2.0 * ta1_z_0_xyyy_0[i] * fe_0 - 2.0 * ta1_z_0_xyyy_1[i] * fe_0 + ta1_z_0_xxyyy_0[i] * pb_x[i] - ta1_z_0_xxyyy_1[i] * pc_x[i];

        ta1_z_0_xxxyyz_0[i] = ta_0_xxxyy_1[i] + ta1_z_0_xxxyy_0[i] * pb_z[i] - ta1_z_0_xxxyy_1[i] * pc_z[i];

        ta1_z_0_xxxyzz_0[i] = ta1_z_0_xxxzz_0[i] * pb_y[i] - ta1_z_0_xxxzz_1[i] * pc_y[i];

        ta1_z_0_xxxzzz_0[i] =
            2.0 * ta1_z_0_xzzz_0[i] * fe_0 - 2.0 * ta1_z_0_xzzz_1[i] * fe_0 + ta1_z_0_xxzzz_0[i] * pb_x[i] - ta1_z_0_xxzzz_1[i] * pc_x[i];

        ta1_z_0_xxyyyy_0[i] = ta1_z_0_yyyy_0[i] * fe_0 - ta1_z_0_yyyy_1[i] * fe_0 + ta1_z_0_xyyyy_0[i] * pb_x[i] - ta1_z_0_xyyyy_1[i] * pc_x[i];

        ta1_z_0_xxyyyz_0[i] = ta_0_xxyyy_1[i] + ta1_z_0_xxyyy_0[i] * pb_z[i] - ta1_z_0_xxyyy_1[i] * pc_z[i];

        ta1_z_0_xxyyzz_0[i] = ta1_z_0_yyzz_0[i] * fe_0 - ta1_z_0_yyzz_1[i] * fe_0 + ta1_z_0_xyyzz_0[i] * pb_x[i] - ta1_z_0_xyyzz_1[i] * pc_x[i];

        ta1_z_0_xxyzzz_0[i] = ta1_z_0_xxzzz_0[i] * pb_y[i] - ta1_z_0_xxzzz_1[i] * pc_y[i];

        ta1_z_0_xxzzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_0_xzzzz_0[i] * pb_x[i] - ta1_z_0_xzzzz_1[i] * pc_x[i];

        ta1_z_0_xyyyyy_0[i] = ta1_z_0_yyyyy_0[i] * pb_x[i] - ta1_z_0_yyyyy_1[i] * pc_x[i];

        ta1_z_0_xyyyyz_0[i] = ta1_z_0_yyyyz_0[i] * pb_x[i] - ta1_z_0_yyyyz_1[i] * pc_x[i];

        ta1_z_0_xyyyzz_0[i] = ta1_z_0_yyyzz_0[i] * pb_x[i] - ta1_z_0_yyyzz_1[i] * pc_x[i];

        ta1_z_0_xyyzzz_0[i] = ta1_z_0_yyzzz_0[i] * pb_x[i] - ta1_z_0_yyzzz_1[i] * pc_x[i];

        ta1_z_0_xyzzzz_0[i] = ta1_z_0_yzzzz_0[i] * pb_x[i] - ta1_z_0_yzzzz_1[i] * pc_x[i];

        ta1_z_0_xzzzzz_0[i] = ta1_z_0_zzzzz_0[i] * pb_x[i] - ta1_z_0_zzzzz_1[i] * pc_x[i];

        ta1_z_0_yyyyyy_0[i] =
            5.0 * ta1_z_0_yyyy_0[i] * fe_0 - 5.0 * ta1_z_0_yyyy_1[i] * fe_0 + ta1_z_0_yyyyy_0[i] * pb_y[i] - ta1_z_0_yyyyy_1[i] * pc_y[i];

        ta1_z_0_yyyyyz_0[i] = ta_0_yyyyy_1[i] + ta1_z_0_yyyyy_0[i] * pb_z[i] - ta1_z_0_yyyyy_1[i] * pc_z[i];

        ta1_z_0_yyyyzz_0[i] =
            3.0 * ta1_z_0_yyzz_0[i] * fe_0 - 3.0 * ta1_z_0_yyzz_1[i] * fe_0 + ta1_z_0_yyyzz_0[i] * pb_y[i] - ta1_z_0_yyyzz_1[i] * pc_y[i];

        ta1_z_0_yyyzzz_0[i] =
            2.0 * ta1_z_0_yzzz_0[i] * fe_0 - 2.0 * ta1_z_0_yzzz_1[i] * fe_0 + ta1_z_0_yyzzz_0[i] * pb_y[i] - ta1_z_0_yyzzz_1[i] * pc_y[i];

        ta1_z_0_yyzzzz_0[i] = ta1_z_0_zzzz_0[i] * fe_0 - ta1_z_0_zzzz_1[i] * fe_0 + ta1_z_0_yzzzz_0[i] * pb_y[i] - ta1_z_0_yzzzz_1[i] * pc_y[i];

        ta1_z_0_yzzzzz_0[i] = ta1_z_0_zzzzz_0[i] * pb_y[i] - ta1_z_0_zzzzz_1[i] * pc_y[i];

        ta1_z_0_zzzzzz_0[i] = 5.0 * ta1_z_0_zzzz_0[i] * fe_0 - 5.0 * ta1_z_0_zzzz_1[i] * fe_0 + ta_0_zzzzz_1[i] + ta1_z_0_zzzzz_0[i] * pb_z[i] -
                              ta1_z_0_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
