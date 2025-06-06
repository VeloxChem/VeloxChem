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

#include "NuclearPotentialGeom010PrimRecIS.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_is(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_is,
                                        const size_t              idx_npot_geom_010_0_gs,
                                        const size_t              idx_npot_geom_010_1_gs,
                                        const size_t              idx_npot_1_hs,
                                        const size_t              idx_npot_geom_010_0_hs,
                                        const size_t              idx_npot_geom_010_1_hs,
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

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs);

    auto ta1_x_xxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 1);

    auto ta1_x_xxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 2);

    auto ta1_x_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 3);

    auto ta1_x_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 5);

    auto ta1_x_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 10);

    auto ta1_x_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 12);

    auto ta1_x_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 13);

    auto ta1_x_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 14);

    auto ta1_y_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 15);

    auto ta1_y_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 18);

    auto ta1_y_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 20);

    auto ta1_y_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 21);

    auto ta1_y_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 24);

    auto ta1_y_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 25);

    auto ta1_y_yyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 26);

    auto ta1_y_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 27);

    auto ta1_y_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 29);

    auto ta1_z_xxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 30);

    auto ta1_z_xxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 33);

    auto ta1_z_xxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 35);

    auto ta1_z_xyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 36);

    auto ta1_z_xzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 39);

    auto ta1_z_yyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 40);

    auto ta1_z_yyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 42);

    auto ta1_z_yzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 43);

    auto ta1_z_zzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_gs + 44);

    // Set up components of auxiliary buffer : GS

    auto ta1_x_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs);

    auto ta1_x_xxxy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 1);

    auto ta1_x_xxxz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 2);

    auto ta1_x_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 3);

    auto ta1_x_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 5);

    auto ta1_x_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 10);

    auto ta1_x_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 12);

    auto ta1_x_yzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 13);

    auto ta1_x_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 14);

    auto ta1_y_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 15);

    auto ta1_y_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 18);

    auto ta1_y_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 20);

    auto ta1_y_xyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 21);

    auto ta1_y_xzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 24);

    auto ta1_y_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 25);

    auto ta1_y_yyyz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 26);

    auto ta1_y_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 27);

    auto ta1_y_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 29);

    auto ta1_z_xxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 30);

    auto ta1_z_xxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 33);

    auto ta1_z_xxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 35);

    auto ta1_z_xyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 36);

    auto ta1_z_xzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 39);

    auto ta1_z_yyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 40);

    auto ta1_z_yyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 42);

    auto ta1_z_yzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 43);

    auto ta1_z_zzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_gs + 44);

    // Set up components of auxiliary buffer : HS

    auto ta_xxxxx_0_1 = pbuffer.data(idx_npot_1_hs);

    auto ta_xxxyy_0_1 = pbuffer.data(idx_npot_1_hs + 3);

    auto ta_xxxzz_0_1 = pbuffer.data(idx_npot_1_hs + 5);

    auto ta_xxyyy_0_1 = pbuffer.data(idx_npot_1_hs + 6);

    auto ta_xxzzz_0_1 = pbuffer.data(idx_npot_1_hs + 9);

    auto ta_yyyyy_0_1 = pbuffer.data(idx_npot_1_hs + 15);

    auto ta_yyyzz_0_1 = pbuffer.data(idx_npot_1_hs + 17);

    auto ta_yyzzz_0_1 = pbuffer.data(idx_npot_1_hs + 18);

    auto ta_zzzzz_0_1 = pbuffer.data(idx_npot_1_hs + 20);

    // Set up components of auxiliary buffer : HS

    auto ta1_x_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs);

    auto ta1_x_xxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 1);

    auto ta1_x_xxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 2);

    auto ta1_x_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 3);

    auto ta1_x_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 5);

    auto ta1_x_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 6);

    auto ta1_x_xxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 8);

    auto ta1_x_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 9);

    auto ta1_x_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 10);

    auto ta1_x_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 14);

    auto ta1_x_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 15);

    auto ta1_x_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 17);

    auto ta1_x_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 18);

    auto ta1_x_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 19);

    auto ta1_x_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 20);

    auto ta1_y_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 21);

    auto ta1_y_xxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 22);

    auto ta1_y_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 24);

    auto ta1_y_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 26);

    auto ta1_y_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 27);

    auto ta1_y_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 30);

    auto ta1_y_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 31);

    auto ta1_y_xyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 33);

    auto ta1_y_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 35);

    auto ta1_y_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 36);

    auto ta1_y_yyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 37);

    auto ta1_y_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 38);

    auto ta1_y_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 39);

    auto ta1_y_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 40);

    auto ta1_y_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 41);

    auto ta1_z_xxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 42);

    auto ta1_z_xxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 44);

    auto ta1_z_xxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 45);

    auto ta1_z_xxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 47);

    auto ta1_z_xxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 48);

    auto ta1_z_xxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 51);

    auto ta1_z_xyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 52);

    auto ta1_z_xyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 54);

    auto ta1_z_xzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 56);

    auto ta1_z_yyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 57);

    auto ta1_z_yyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 58);

    auto ta1_z_yyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 59);

    auto ta1_z_yyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 60);

    auto ta1_z_yzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 61);

    auto ta1_z_zzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_hs + 62);

    // Set up components of auxiliary buffer : HS

    auto ta1_x_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs);

    auto ta1_x_xxxxy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 1);

    auto ta1_x_xxxxz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 2);

    auto ta1_x_xxxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 3);

    auto ta1_x_xxxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 5);

    auto ta1_x_xxyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 6);

    auto ta1_x_xxyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 8);

    auto ta1_x_xxzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 9);

    auto ta1_x_xyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 10);

    auto ta1_x_xzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 14);

    auto ta1_x_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 15);

    auto ta1_x_yyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 17);

    auto ta1_x_yyzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 18);

    auto ta1_x_yzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 19);

    auto ta1_x_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 20);

    auto ta1_y_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 21);

    auto ta1_y_xxxxy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 22);

    auto ta1_y_xxxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 24);

    auto ta1_y_xxxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 26);

    auto ta1_y_xxyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 27);

    auto ta1_y_xxzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 30);

    auto ta1_y_xyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 31);

    auto ta1_y_xyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 33);

    auto ta1_y_xzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 35);

    auto ta1_y_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 36);

    auto ta1_y_yyyyz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 37);

    auto ta1_y_yyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 38);

    auto ta1_y_yyzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 39);

    auto ta1_y_yzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 40);

    auto ta1_y_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 41);

    auto ta1_z_xxxxx_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 42);

    auto ta1_z_xxxxz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 44);

    auto ta1_z_xxxyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 45);

    auto ta1_z_xxxzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 47);

    auto ta1_z_xxyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 48);

    auto ta1_z_xxzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 51);

    auto ta1_z_xyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 52);

    auto ta1_z_xyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 54);

    auto ta1_z_xzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 56);

    auto ta1_z_yyyyy_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 57);

    auto ta1_z_yyyyz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 58);

    auto ta1_z_yyyzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 59);

    auto ta1_z_yyzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 60);

    auto ta1_z_yzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 61);

    auto ta1_z_zzzzz_0_1 = pbuffer.data(idx_npot_geom_010_1_hs + 62);

    // Set up components of targeted buffer : IS

    auto ta1_x_xxxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_is);

    auto ta1_x_xxxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 1);

    auto ta1_x_xxxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 2);

    auto ta1_x_xxxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 3);

    auto ta1_x_xxxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 4);

    auto ta1_x_xxxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 5);

    auto ta1_x_xxxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 6);

    auto ta1_x_xxxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 7);

    auto ta1_x_xxxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 8);

    auto ta1_x_xxxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 9);

    auto ta1_x_xxyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 10);

    auto ta1_x_xxyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 11);

    auto ta1_x_xxyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 12);

    auto ta1_x_xxyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 13);

    auto ta1_x_xxzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 14);

    auto ta1_x_xyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 15);

    auto ta1_x_xyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 16);

    auto ta1_x_xyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 17);

    auto ta1_x_xyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 18);

    auto ta1_x_xyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 19);

    auto ta1_x_xzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 20);

    auto ta1_x_yyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 21);

    auto ta1_x_yyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 22);

    auto ta1_x_yyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 23);

    auto ta1_x_yyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 24);

    auto ta1_x_yyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 25);

    auto ta1_x_yzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 26);

    auto ta1_x_zzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 27);

    auto ta1_y_xxxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 28);

    auto ta1_y_xxxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 29);

    auto ta1_y_xxxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 30);

    auto ta1_y_xxxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 31);

    auto ta1_y_xxxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 32);

    auto ta1_y_xxxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 33);

    auto ta1_y_xxxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 34);

    auto ta1_y_xxxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 35);

    auto ta1_y_xxxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 36);

    auto ta1_y_xxxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 37);

    auto ta1_y_xxyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 38);

    auto ta1_y_xxyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 39);

    auto ta1_y_xxyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 40);

    auto ta1_y_xxyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 41);

    auto ta1_y_xxzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 42);

    auto ta1_y_xyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 43);

    auto ta1_y_xyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 44);

    auto ta1_y_xyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 45);

    auto ta1_y_xyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 46);

    auto ta1_y_xyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 47);

    auto ta1_y_xzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 48);

    auto ta1_y_yyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 49);

    auto ta1_y_yyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 50);

    auto ta1_y_yyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 51);

    auto ta1_y_yyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 52);

    auto ta1_y_yyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 53);

    auto ta1_y_yzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 54);

    auto ta1_y_zzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 55);

    auto ta1_z_xxxxxx_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 56);

    auto ta1_z_xxxxxy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 57);

    auto ta1_z_xxxxxz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 58);

    auto ta1_z_xxxxyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 59);

    auto ta1_z_xxxxyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 60);

    auto ta1_z_xxxxzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 61);

    auto ta1_z_xxxyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 62);

    auto ta1_z_xxxyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 63);

    auto ta1_z_xxxyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 64);

    auto ta1_z_xxxzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 65);

    auto ta1_z_xxyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 66);

    auto ta1_z_xxyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 67);

    auto ta1_z_xxyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 68);

    auto ta1_z_xxyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 69);

    auto ta1_z_xxzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 70);

    auto ta1_z_xyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 71);

    auto ta1_z_xyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 72);

    auto ta1_z_xyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 73);

    auto ta1_z_xyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 74);

    auto ta1_z_xyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 75);

    auto ta1_z_xzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 76);

    auto ta1_z_yyyyyy_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 77);

    auto ta1_z_yyyyyz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 78);

    auto ta1_z_yyyyzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 79);

    auto ta1_z_yyyzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 80);

    auto ta1_z_yyzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 81);

    auto ta1_z_yzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 82);

    auto ta1_z_zzzzzz_0_0 = pbuffer.data(idx_npot_geom_010_0_is + 83);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxxx_0_0,   \
                             ta1_x_xxxx_0_1,   \
                             ta1_x_xxxxx_0_0,  \
                             ta1_x_xxxxx_0_1,  \
                             ta1_x_xxxxxx_0_0, \
                             ta1_x_xxxxxy_0_0, \
                             ta1_x_xxxxxz_0_0, \
                             ta1_x_xxxxy_0_0,  \
                             ta1_x_xxxxy_0_1,  \
                             ta1_x_xxxxyy_0_0, \
                             ta1_x_xxxxyz_0_0, \
                             ta1_x_xxxxz_0_0,  \
                             ta1_x_xxxxz_0_1,  \
                             ta1_x_xxxxzz_0_0, \
                             ta1_x_xxxy_0_0,   \
                             ta1_x_xxxy_0_1,   \
                             ta1_x_xxxyy_0_0,  \
                             ta1_x_xxxyy_0_1,  \
                             ta1_x_xxxyyy_0_0, \
                             ta1_x_xxxyyz_0_0, \
                             ta1_x_xxxyzz_0_0, \
                             ta1_x_xxxz_0_0,   \
                             ta1_x_xxxz_0_1,   \
                             ta1_x_xxxzz_0_0,  \
                             ta1_x_xxxzz_0_1,  \
                             ta1_x_xxxzzz_0_0, \
                             ta1_x_xxyy_0_0,   \
                             ta1_x_xxyy_0_1,   \
                             ta1_x_xxyyy_0_0,  \
                             ta1_x_xxyyy_0_1,  \
                             ta1_x_xxyyyy_0_0, \
                             ta1_x_xxyyyz_0_0, \
                             ta1_x_xxyyzz_0_0, \
                             ta1_x_xxyzz_0_0,  \
                             ta1_x_xxyzz_0_1,  \
                             ta1_x_xxyzzz_0_0, \
                             ta1_x_xxzz_0_0,   \
                             ta1_x_xxzz_0_1,   \
                             ta1_x_xxzzz_0_0,  \
                             ta1_x_xxzzz_0_1,  \
                             ta1_x_xxzzzz_0_0, \
                             ta1_x_xyyyy_0_0,  \
                             ta1_x_xyyyy_0_1,  \
                             ta1_x_xyyyyy_0_0, \
                             ta1_x_xyyyyz_0_0, \
                             ta1_x_xyyyzz_0_0, \
                             ta1_x_xyyzzz_0_0, \
                             ta1_x_xyzzzz_0_0, \
                             ta1_x_xzzzz_0_0,  \
                             ta1_x_xzzzz_0_1,  \
                             ta1_x_xzzzzz_0_0, \
                             ta1_x_yyyy_0_0,   \
                             ta1_x_yyyy_0_1,   \
                             ta1_x_yyyyy_0_0,  \
                             ta1_x_yyyyy_0_1,  \
                             ta1_x_yyyyyy_0_0, \
                             ta1_x_yyyyyz_0_0, \
                             ta1_x_yyyyzz_0_0, \
                             ta1_x_yyyzz_0_0,  \
                             ta1_x_yyyzz_0_1,  \
                             ta1_x_yyyzzz_0_0, \
                             ta1_x_yyzz_0_0,   \
                             ta1_x_yyzz_0_1,   \
                             ta1_x_yyzzz_0_0,  \
                             ta1_x_yyzzz_0_1,  \
                             ta1_x_yyzzzz_0_0, \
                             ta1_x_yzzz_0_0,   \
                             ta1_x_yzzz_0_1,   \
                             ta1_x_yzzzz_0_0,  \
                             ta1_x_yzzzz_0_1,  \
                             ta1_x_yzzzzz_0_0, \
                             ta1_x_zzzz_0_0,   \
                             ta1_x_zzzz_0_1,   \
                             ta1_x_zzzzz_0_0,  \
                             ta1_x_zzzzz_0_1,  \
                             ta1_x_zzzzzz_0_0, \
                             ta1_y_xxxx_0_0,   \
                             ta1_y_xxxx_0_1,   \
                             ta1_y_xxxxx_0_0,  \
                             ta1_y_xxxxx_0_1,  \
                             ta1_y_xxxxxx_0_0, \
                             ta1_y_xxxxxy_0_0, \
                             ta1_y_xxxxxz_0_0, \
                             ta1_y_xxxxy_0_0,  \
                             ta1_y_xxxxy_0_1,  \
                             ta1_y_xxxxyy_0_0, \
                             ta1_y_xxxxyz_0_0, \
                             ta1_y_xxxxzz_0_0, \
                             ta1_y_xxxyy_0_0,  \
                             ta1_y_xxxyy_0_1,  \
                             ta1_y_xxxyyy_0_0, \
                             ta1_y_xxxyyz_0_0, \
                             ta1_y_xxxyzz_0_0, \
                             ta1_y_xxxzz_0_0,  \
                             ta1_y_xxxzz_0_1,  \
                             ta1_y_xxxzzz_0_0, \
                             ta1_y_xxyy_0_0,   \
                             ta1_y_xxyy_0_1,   \
                             ta1_y_xxyyy_0_0,  \
                             ta1_y_xxyyy_0_1,  \
                             ta1_y_xxyyyy_0_0, \
                             ta1_y_xxyyyz_0_0, \
                             ta1_y_xxyyzz_0_0, \
                             ta1_y_xxyzzz_0_0, \
                             ta1_y_xxzz_0_0,   \
                             ta1_y_xxzz_0_1,   \
                             ta1_y_xxzzz_0_0,  \
                             ta1_y_xxzzz_0_1,  \
                             ta1_y_xxzzzz_0_0, \
                             ta1_y_xyyy_0_0,   \
                             ta1_y_xyyy_0_1,   \
                             ta1_y_xyyyy_0_0,  \
                             ta1_y_xyyyy_0_1,  \
                             ta1_y_xyyyyy_0_0, \
                             ta1_y_xyyyyz_0_0, \
                             ta1_y_xyyyzz_0_0, \
                             ta1_y_xyyzz_0_0,  \
                             ta1_y_xyyzz_0_1,  \
                             ta1_y_xyyzzz_0_0, \
                             ta1_y_xyzzzz_0_0, \
                             ta1_y_xzzz_0_0,   \
                             ta1_y_xzzz_0_1,   \
                             ta1_y_xzzzz_0_0,  \
                             ta1_y_xzzzz_0_1,  \
                             ta1_y_xzzzzz_0_0, \
                             ta1_y_yyyy_0_0,   \
                             ta1_y_yyyy_0_1,   \
                             ta1_y_yyyyy_0_0,  \
                             ta1_y_yyyyy_0_1,  \
                             ta1_y_yyyyyy_0_0, \
                             ta1_y_yyyyyz_0_0, \
                             ta1_y_yyyyz_0_0,  \
                             ta1_y_yyyyz_0_1,  \
                             ta1_y_yyyyzz_0_0, \
                             ta1_y_yyyz_0_0,   \
                             ta1_y_yyyz_0_1,   \
                             ta1_y_yyyzz_0_0,  \
                             ta1_y_yyyzz_0_1,  \
                             ta1_y_yyyzzz_0_0, \
                             ta1_y_yyzz_0_0,   \
                             ta1_y_yyzz_0_1,   \
                             ta1_y_yyzzz_0_0,  \
                             ta1_y_yyzzz_0_1,  \
                             ta1_y_yyzzzz_0_0, \
                             ta1_y_yzzzz_0_0,  \
                             ta1_y_yzzzz_0_1,  \
                             ta1_y_yzzzzz_0_0, \
                             ta1_y_zzzz_0_0,   \
                             ta1_y_zzzz_0_1,   \
                             ta1_y_zzzzz_0_0,  \
                             ta1_y_zzzzz_0_1,  \
                             ta1_y_zzzzzz_0_0, \
                             ta1_z_xxxx_0_0,   \
                             ta1_z_xxxx_0_1,   \
                             ta1_z_xxxxx_0_0,  \
                             ta1_z_xxxxx_0_1,  \
                             ta1_z_xxxxxx_0_0, \
                             ta1_z_xxxxxy_0_0, \
                             ta1_z_xxxxxz_0_0, \
                             ta1_z_xxxxyy_0_0, \
                             ta1_z_xxxxyz_0_0, \
                             ta1_z_xxxxz_0_0,  \
                             ta1_z_xxxxz_0_1,  \
                             ta1_z_xxxxzz_0_0, \
                             ta1_z_xxxyy_0_0,  \
                             ta1_z_xxxyy_0_1,  \
                             ta1_z_xxxyyy_0_0, \
                             ta1_z_xxxyyz_0_0, \
                             ta1_z_xxxyzz_0_0, \
                             ta1_z_xxxzz_0_0,  \
                             ta1_z_xxxzz_0_1,  \
                             ta1_z_xxxzzz_0_0, \
                             ta1_z_xxyy_0_0,   \
                             ta1_z_xxyy_0_1,   \
                             ta1_z_xxyyy_0_0,  \
                             ta1_z_xxyyy_0_1,  \
                             ta1_z_xxyyyy_0_0, \
                             ta1_z_xxyyyz_0_0, \
                             ta1_z_xxyyzz_0_0, \
                             ta1_z_xxyzzz_0_0, \
                             ta1_z_xxzz_0_0,   \
                             ta1_z_xxzz_0_1,   \
                             ta1_z_xxzzz_0_0,  \
                             ta1_z_xxzzz_0_1,  \
                             ta1_z_xxzzzz_0_0, \
                             ta1_z_xyyy_0_0,   \
                             ta1_z_xyyy_0_1,   \
                             ta1_z_xyyyy_0_0,  \
                             ta1_z_xyyyy_0_1,  \
                             ta1_z_xyyyyy_0_0, \
                             ta1_z_xyyyyz_0_0, \
                             ta1_z_xyyyzz_0_0, \
                             ta1_z_xyyzz_0_0,  \
                             ta1_z_xyyzz_0_1,  \
                             ta1_z_xyyzzz_0_0, \
                             ta1_z_xyzzzz_0_0, \
                             ta1_z_xzzz_0_0,   \
                             ta1_z_xzzz_0_1,   \
                             ta1_z_xzzzz_0_0,  \
                             ta1_z_xzzzz_0_1,  \
                             ta1_z_xzzzzz_0_0, \
                             ta1_z_yyyy_0_0,   \
                             ta1_z_yyyy_0_1,   \
                             ta1_z_yyyyy_0_0,  \
                             ta1_z_yyyyy_0_1,  \
                             ta1_z_yyyyyy_0_0, \
                             ta1_z_yyyyyz_0_0, \
                             ta1_z_yyyyz_0_0,  \
                             ta1_z_yyyyz_0_1,  \
                             ta1_z_yyyyzz_0_0, \
                             ta1_z_yyyzz_0_0,  \
                             ta1_z_yyyzz_0_1,  \
                             ta1_z_yyyzzz_0_0, \
                             ta1_z_yyzz_0_0,   \
                             ta1_z_yyzz_0_1,   \
                             ta1_z_yyzzz_0_0,  \
                             ta1_z_yyzzz_0_1,  \
                             ta1_z_yyzzzz_0_0, \
                             ta1_z_yzzz_0_0,   \
                             ta1_z_yzzz_0_1,   \
                             ta1_z_yzzzz_0_0,  \
                             ta1_z_yzzzz_0_1,  \
                             ta1_z_yzzzzz_0_0, \
                             ta1_z_zzzz_0_0,   \
                             ta1_z_zzzz_0_1,   \
                             ta1_z_zzzzz_0_0,  \
                             ta1_z_zzzzz_0_1,  \
                             ta1_z_zzzzzz_0_0, \
                             ta_xxxxx_0_1,     \
                             ta_xxxyy_0_1,     \
                             ta_xxxzz_0_1,     \
                             ta_xxyyy_0_1,     \
                             ta_xxzzz_0_1,     \
                             ta_yyyyy_0_1,     \
                             ta_yyyzz_0_1,     \
                             ta_yyzzz_0_1,     \
                             ta_zzzzz_0_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxx_0_0[i] = 5.0 * ta1_x_xxxx_0_0[i] * fe_0 - 5.0 * ta1_x_xxxx_0_1[i] * fe_0 + ta_xxxxx_0_1[i] + ta1_x_xxxxx_0_0[i] * pa_x[i] -
                              ta1_x_xxxxx_0_1[i] * pc_x[i];

        ta1_x_xxxxxy_0_0[i] = ta1_x_xxxxx_0_0[i] * pa_y[i] - ta1_x_xxxxx_0_1[i] * pc_y[i];

        ta1_x_xxxxxz_0_0[i] = ta1_x_xxxxx_0_0[i] * pa_z[i] - ta1_x_xxxxx_0_1[i] * pc_z[i];

        ta1_x_xxxxyy_0_0[i] = ta1_x_xxxx_0_0[i] * fe_0 - ta1_x_xxxx_0_1[i] * fe_0 + ta1_x_xxxxy_0_0[i] * pa_y[i] - ta1_x_xxxxy_0_1[i] * pc_y[i];

        ta1_x_xxxxyz_0_0[i] = ta1_x_xxxxz_0_0[i] * pa_y[i] - ta1_x_xxxxz_0_1[i] * pc_y[i];

        ta1_x_xxxxzz_0_0[i] = ta1_x_xxxx_0_0[i] * fe_0 - ta1_x_xxxx_0_1[i] * fe_0 + ta1_x_xxxxz_0_0[i] * pa_z[i] - ta1_x_xxxxz_0_1[i] * pc_z[i];

        ta1_x_xxxyyy_0_0[i] =
            2.0 * ta1_x_xxxy_0_0[i] * fe_0 - 2.0 * ta1_x_xxxy_0_1[i] * fe_0 + ta1_x_xxxyy_0_0[i] * pa_y[i] - ta1_x_xxxyy_0_1[i] * pc_y[i];

        ta1_x_xxxyyz_0_0[i] = ta1_x_xxxyy_0_0[i] * pa_z[i] - ta1_x_xxxyy_0_1[i] * pc_z[i];

        ta1_x_xxxyzz_0_0[i] = ta1_x_xxxzz_0_0[i] * pa_y[i] - ta1_x_xxxzz_0_1[i] * pc_y[i];

        ta1_x_xxxzzz_0_0[i] =
            2.0 * ta1_x_xxxz_0_0[i] * fe_0 - 2.0 * ta1_x_xxxz_0_1[i] * fe_0 + ta1_x_xxxzz_0_0[i] * pa_z[i] - ta1_x_xxxzz_0_1[i] * pc_z[i];

        ta1_x_xxyyyy_0_0[i] =
            3.0 * ta1_x_xxyy_0_0[i] * fe_0 - 3.0 * ta1_x_xxyy_0_1[i] * fe_0 + ta1_x_xxyyy_0_0[i] * pa_y[i] - ta1_x_xxyyy_0_1[i] * pc_y[i];

        ta1_x_xxyyyz_0_0[i] = ta1_x_xxyyy_0_0[i] * pa_z[i] - ta1_x_xxyyy_0_1[i] * pc_z[i];

        ta1_x_xxyyzz_0_0[i] = ta1_x_xxzz_0_0[i] * fe_0 - ta1_x_xxzz_0_1[i] * fe_0 + ta1_x_xxyzz_0_0[i] * pa_y[i] - ta1_x_xxyzz_0_1[i] * pc_y[i];

        ta1_x_xxyzzz_0_0[i] = ta1_x_xxzzz_0_0[i] * pa_y[i] - ta1_x_xxzzz_0_1[i] * pc_y[i];

        ta1_x_xxzzzz_0_0[i] =
            3.0 * ta1_x_xxzz_0_0[i] * fe_0 - 3.0 * ta1_x_xxzz_0_1[i] * fe_0 + ta1_x_xxzzz_0_0[i] * pa_z[i] - ta1_x_xxzzz_0_1[i] * pc_z[i];

        ta1_x_xyyyyy_0_0[i] = ta_yyyyy_0_1[i] + ta1_x_yyyyy_0_0[i] * pa_x[i] - ta1_x_yyyyy_0_1[i] * pc_x[i];

        ta1_x_xyyyyz_0_0[i] = ta1_x_xyyyy_0_0[i] * pa_z[i] - ta1_x_xyyyy_0_1[i] * pc_z[i];

        ta1_x_xyyyzz_0_0[i] = ta_yyyzz_0_1[i] + ta1_x_yyyzz_0_0[i] * pa_x[i] - ta1_x_yyyzz_0_1[i] * pc_x[i];

        ta1_x_xyyzzz_0_0[i] = ta_yyzzz_0_1[i] + ta1_x_yyzzz_0_0[i] * pa_x[i] - ta1_x_yyzzz_0_1[i] * pc_x[i];

        ta1_x_xyzzzz_0_0[i] = ta1_x_xzzzz_0_0[i] * pa_y[i] - ta1_x_xzzzz_0_1[i] * pc_y[i];

        ta1_x_xzzzzz_0_0[i] = ta_zzzzz_0_1[i] + ta1_x_zzzzz_0_0[i] * pa_x[i] - ta1_x_zzzzz_0_1[i] * pc_x[i];

        ta1_x_yyyyyy_0_0[i] =
            5.0 * ta1_x_yyyy_0_0[i] * fe_0 - 5.0 * ta1_x_yyyy_0_1[i] * fe_0 + ta1_x_yyyyy_0_0[i] * pa_y[i] - ta1_x_yyyyy_0_1[i] * pc_y[i];

        ta1_x_yyyyyz_0_0[i] = ta1_x_yyyyy_0_0[i] * pa_z[i] - ta1_x_yyyyy_0_1[i] * pc_z[i];

        ta1_x_yyyyzz_0_0[i] =
            3.0 * ta1_x_yyzz_0_0[i] * fe_0 - 3.0 * ta1_x_yyzz_0_1[i] * fe_0 + ta1_x_yyyzz_0_0[i] * pa_y[i] - ta1_x_yyyzz_0_1[i] * pc_y[i];

        ta1_x_yyyzzz_0_0[i] =
            2.0 * ta1_x_yzzz_0_0[i] * fe_0 - 2.0 * ta1_x_yzzz_0_1[i] * fe_0 + ta1_x_yyzzz_0_0[i] * pa_y[i] - ta1_x_yyzzz_0_1[i] * pc_y[i];

        ta1_x_yyzzzz_0_0[i] = ta1_x_zzzz_0_0[i] * fe_0 - ta1_x_zzzz_0_1[i] * fe_0 + ta1_x_yzzzz_0_0[i] * pa_y[i] - ta1_x_yzzzz_0_1[i] * pc_y[i];

        ta1_x_yzzzzz_0_0[i] = ta1_x_zzzzz_0_0[i] * pa_y[i] - ta1_x_zzzzz_0_1[i] * pc_y[i];

        ta1_x_zzzzzz_0_0[i] =
            5.0 * ta1_x_zzzz_0_0[i] * fe_0 - 5.0 * ta1_x_zzzz_0_1[i] * fe_0 + ta1_x_zzzzz_0_0[i] * pa_z[i] - ta1_x_zzzzz_0_1[i] * pc_z[i];

        ta1_y_xxxxxx_0_0[i] =
            5.0 * ta1_y_xxxx_0_0[i] * fe_0 - 5.0 * ta1_y_xxxx_0_1[i] * fe_0 + ta1_y_xxxxx_0_0[i] * pa_x[i] - ta1_y_xxxxx_0_1[i] * pc_x[i];

        ta1_y_xxxxxy_0_0[i] = ta_xxxxx_0_1[i] + ta1_y_xxxxx_0_0[i] * pa_y[i] - ta1_y_xxxxx_0_1[i] * pc_y[i];

        ta1_y_xxxxxz_0_0[i] = ta1_y_xxxxx_0_0[i] * pa_z[i] - ta1_y_xxxxx_0_1[i] * pc_z[i];

        ta1_y_xxxxyy_0_0[i] =
            3.0 * ta1_y_xxyy_0_0[i] * fe_0 - 3.0 * ta1_y_xxyy_0_1[i] * fe_0 + ta1_y_xxxyy_0_0[i] * pa_x[i] - ta1_y_xxxyy_0_1[i] * pc_x[i];

        ta1_y_xxxxyz_0_0[i] = ta1_y_xxxxy_0_0[i] * pa_z[i] - ta1_y_xxxxy_0_1[i] * pc_z[i];

        ta1_y_xxxxzz_0_0[i] =
            3.0 * ta1_y_xxzz_0_0[i] * fe_0 - 3.0 * ta1_y_xxzz_0_1[i] * fe_0 + ta1_y_xxxzz_0_0[i] * pa_x[i] - ta1_y_xxxzz_0_1[i] * pc_x[i];

        ta1_y_xxxyyy_0_0[i] =
            2.0 * ta1_y_xyyy_0_0[i] * fe_0 - 2.0 * ta1_y_xyyy_0_1[i] * fe_0 + ta1_y_xxyyy_0_0[i] * pa_x[i] - ta1_y_xxyyy_0_1[i] * pc_x[i];

        ta1_y_xxxyyz_0_0[i] = ta1_y_xxxyy_0_0[i] * pa_z[i] - ta1_y_xxxyy_0_1[i] * pc_z[i];

        ta1_y_xxxyzz_0_0[i] = ta_xxxzz_0_1[i] + ta1_y_xxxzz_0_0[i] * pa_y[i] - ta1_y_xxxzz_0_1[i] * pc_y[i];

        ta1_y_xxxzzz_0_0[i] =
            2.0 * ta1_y_xzzz_0_0[i] * fe_0 - 2.0 * ta1_y_xzzz_0_1[i] * fe_0 + ta1_y_xxzzz_0_0[i] * pa_x[i] - ta1_y_xxzzz_0_1[i] * pc_x[i];

        ta1_y_xxyyyy_0_0[i] = ta1_y_yyyy_0_0[i] * fe_0 - ta1_y_yyyy_0_1[i] * fe_0 + ta1_y_xyyyy_0_0[i] * pa_x[i] - ta1_y_xyyyy_0_1[i] * pc_x[i];

        ta1_y_xxyyyz_0_0[i] = ta1_y_xxyyy_0_0[i] * pa_z[i] - ta1_y_xxyyy_0_1[i] * pc_z[i];

        ta1_y_xxyyzz_0_0[i] = ta1_y_yyzz_0_0[i] * fe_0 - ta1_y_yyzz_0_1[i] * fe_0 + ta1_y_xyyzz_0_0[i] * pa_x[i] - ta1_y_xyyzz_0_1[i] * pc_x[i];

        ta1_y_xxyzzz_0_0[i] = ta_xxzzz_0_1[i] + ta1_y_xxzzz_0_0[i] * pa_y[i] - ta1_y_xxzzz_0_1[i] * pc_y[i];

        ta1_y_xxzzzz_0_0[i] = ta1_y_zzzz_0_0[i] * fe_0 - ta1_y_zzzz_0_1[i] * fe_0 + ta1_y_xzzzz_0_0[i] * pa_x[i] - ta1_y_xzzzz_0_1[i] * pc_x[i];

        ta1_y_xyyyyy_0_0[i] = ta1_y_yyyyy_0_0[i] * pa_x[i] - ta1_y_yyyyy_0_1[i] * pc_x[i];

        ta1_y_xyyyyz_0_0[i] = ta1_y_yyyyz_0_0[i] * pa_x[i] - ta1_y_yyyyz_0_1[i] * pc_x[i];

        ta1_y_xyyyzz_0_0[i] = ta1_y_yyyzz_0_0[i] * pa_x[i] - ta1_y_yyyzz_0_1[i] * pc_x[i];

        ta1_y_xyyzzz_0_0[i] = ta1_y_yyzzz_0_0[i] * pa_x[i] - ta1_y_yyzzz_0_1[i] * pc_x[i];

        ta1_y_xyzzzz_0_0[i] = ta1_y_yzzzz_0_0[i] * pa_x[i] - ta1_y_yzzzz_0_1[i] * pc_x[i];

        ta1_y_xzzzzz_0_0[i] = ta1_y_zzzzz_0_0[i] * pa_x[i] - ta1_y_zzzzz_0_1[i] * pc_x[i];

        ta1_y_yyyyyy_0_0[i] = 5.0 * ta1_y_yyyy_0_0[i] * fe_0 - 5.0 * ta1_y_yyyy_0_1[i] * fe_0 + ta_yyyyy_0_1[i] + ta1_y_yyyyy_0_0[i] * pa_y[i] -
                              ta1_y_yyyyy_0_1[i] * pc_y[i];

        ta1_y_yyyyyz_0_0[i] = ta1_y_yyyyy_0_0[i] * pa_z[i] - ta1_y_yyyyy_0_1[i] * pc_z[i];

        ta1_y_yyyyzz_0_0[i] = ta1_y_yyyy_0_0[i] * fe_0 - ta1_y_yyyy_0_1[i] * fe_0 + ta1_y_yyyyz_0_0[i] * pa_z[i] - ta1_y_yyyyz_0_1[i] * pc_z[i];

        ta1_y_yyyzzz_0_0[i] =
            2.0 * ta1_y_yyyz_0_0[i] * fe_0 - 2.0 * ta1_y_yyyz_0_1[i] * fe_0 + ta1_y_yyyzz_0_0[i] * pa_z[i] - ta1_y_yyyzz_0_1[i] * pc_z[i];

        ta1_y_yyzzzz_0_0[i] =
            3.0 * ta1_y_yyzz_0_0[i] * fe_0 - 3.0 * ta1_y_yyzz_0_1[i] * fe_0 + ta1_y_yyzzz_0_0[i] * pa_z[i] - ta1_y_yyzzz_0_1[i] * pc_z[i];

        ta1_y_yzzzzz_0_0[i] = ta_zzzzz_0_1[i] + ta1_y_zzzzz_0_0[i] * pa_y[i] - ta1_y_zzzzz_0_1[i] * pc_y[i];

        ta1_y_zzzzzz_0_0[i] =
            5.0 * ta1_y_zzzz_0_0[i] * fe_0 - 5.0 * ta1_y_zzzz_0_1[i] * fe_0 + ta1_y_zzzzz_0_0[i] * pa_z[i] - ta1_y_zzzzz_0_1[i] * pc_z[i];

        ta1_z_xxxxxx_0_0[i] =
            5.0 * ta1_z_xxxx_0_0[i] * fe_0 - 5.0 * ta1_z_xxxx_0_1[i] * fe_0 + ta1_z_xxxxx_0_0[i] * pa_x[i] - ta1_z_xxxxx_0_1[i] * pc_x[i];

        ta1_z_xxxxxy_0_0[i] = ta1_z_xxxxx_0_0[i] * pa_y[i] - ta1_z_xxxxx_0_1[i] * pc_y[i];

        ta1_z_xxxxxz_0_0[i] = ta_xxxxx_0_1[i] + ta1_z_xxxxx_0_0[i] * pa_z[i] - ta1_z_xxxxx_0_1[i] * pc_z[i];

        ta1_z_xxxxyy_0_0[i] =
            3.0 * ta1_z_xxyy_0_0[i] * fe_0 - 3.0 * ta1_z_xxyy_0_1[i] * fe_0 + ta1_z_xxxyy_0_0[i] * pa_x[i] - ta1_z_xxxyy_0_1[i] * pc_x[i];

        ta1_z_xxxxyz_0_0[i] = ta1_z_xxxxz_0_0[i] * pa_y[i] - ta1_z_xxxxz_0_1[i] * pc_y[i];

        ta1_z_xxxxzz_0_0[i] =
            3.0 * ta1_z_xxzz_0_0[i] * fe_0 - 3.0 * ta1_z_xxzz_0_1[i] * fe_0 + ta1_z_xxxzz_0_0[i] * pa_x[i] - ta1_z_xxxzz_0_1[i] * pc_x[i];

        ta1_z_xxxyyy_0_0[i] =
            2.0 * ta1_z_xyyy_0_0[i] * fe_0 - 2.0 * ta1_z_xyyy_0_1[i] * fe_0 + ta1_z_xxyyy_0_0[i] * pa_x[i] - ta1_z_xxyyy_0_1[i] * pc_x[i];

        ta1_z_xxxyyz_0_0[i] = ta_xxxyy_0_1[i] + ta1_z_xxxyy_0_0[i] * pa_z[i] - ta1_z_xxxyy_0_1[i] * pc_z[i];

        ta1_z_xxxyzz_0_0[i] = ta1_z_xxxzz_0_0[i] * pa_y[i] - ta1_z_xxxzz_0_1[i] * pc_y[i];

        ta1_z_xxxzzz_0_0[i] =
            2.0 * ta1_z_xzzz_0_0[i] * fe_0 - 2.0 * ta1_z_xzzz_0_1[i] * fe_0 + ta1_z_xxzzz_0_0[i] * pa_x[i] - ta1_z_xxzzz_0_1[i] * pc_x[i];

        ta1_z_xxyyyy_0_0[i] = ta1_z_yyyy_0_0[i] * fe_0 - ta1_z_yyyy_0_1[i] * fe_0 + ta1_z_xyyyy_0_0[i] * pa_x[i] - ta1_z_xyyyy_0_1[i] * pc_x[i];

        ta1_z_xxyyyz_0_0[i] = ta_xxyyy_0_1[i] + ta1_z_xxyyy_0_0[i] * pa_z[i] - ta1_z_xxyyy_0_1[i] * pc_z[i];

        ta1_z_xxyyzz_0_0[i] = ta1_z_yyzz_0_0[i] * fe_0 - ta1_z_yyzz_0_1[i] * fe_0 + ta1_z_xyyzz_0_0[i] * pa_x[i] - ta1_z_xyyzz_0_1[i] * pc_x[i];

        ta1_z_xxyzzz_0_0[i] = ta1_z_xxzzz_0_0[i] * pa_y[i] - ta1_z_xxzzz_0_1[i] * pc_y[i];

        ta1_z_xxzzzz_0_0[i] = ta1_z_zzzz_0_0[i] * fe_0 - ta1_z_zzzz_0_1[i] * fe_0 + ta1_z_xzzzz_0_0[i] * pa_x[i] - ta1_z_xzzzz_0_1[i] * pc_x[i];

        ta1_z_xyyyyy_0_0[i] = ta1_z_yyyyy_0_0[i] * pa_x[i] - ta1_z_yyyyy_0_1[i] * pc_x[i];

        ta1_z_xyyyyz_0_0[i] = ta1_z_yyyyz_0_0[i] * pa_x[i] - ta1_z_yyyyz_0_1[i] * pc_x[i];

        ta1_z_xyyyzz_0_0[i] = ta1_z_yyyzz_0_0[i] * pa_x[i] - ta1_z_yyyzz_0_1[i] * pc_x[i];

        ta1_z_xyyzzz_0_0[i] = ta1_z_yyzzz_0_0[i] * pa_x[i] - ta1_z_yyzzz_0_1[i] * pc_x[i];

        ta1_z_xyzzzz_0_0[i] = ta1_z_yzzzz_0_0[i] * pa_x[i] - ta1_z_yzzzz_0_1[i] * pc_x[i];

        ta1_z_xzzzzz_0_0[i] = ta1_z_zzzzz_0_0[i] * pa_x[i] - ta1_z_zzzzz_0_1[i] * pc_x[i];

        ta1_z_yyyyyy_0_0[i] =
            5.0 * ta1_z_yyyy_0_0[i] * fe_0 - 5.0 * ta1_z_yyyy_0_1[i] * fe_0 + ta1_z_yyyyy_0_0[i] * pa_y[i] - ta1_z_yyyyy_0_1[i] * pc_y[i];

        ta1_z_yyyyyz_0_0[i] = ta_yyyyy_0_1[i] + ta1_z_yyyyy_0_0[i] * pa_z[i] - ta1_z_yyyyy_0_1[i] * pc_z[i];

        ta1_z_yyyyzz_0_0[i] =
            3.0 * ta1_z_yyzz_0_0[i] * fe_0 - 3.0 * ta1_z_yyzz_0_1[i] * fe_0 + ta1_z_yyyzz_0_0[i] * pa_y[i] - ta1_z_yyyzz_0_1[i] * pc_y[i];

        ta1_z_yyyzzz_0_0[i] =
            2.0 * ta1_z_yzzz_0_0[i] * fe_0 - 2.0 * ta1_z_yzzz_0_1[i] * fe_0 + ta1_z_yyzzz_0_0[i] * pa_y[i] - ta1_z_yyzzz_0_1[i] * pc_y[i];

        ta1_z_yyzzzz_0_0[i] = ta1_z_zzzz_0_0[i] * fe_0 - ta1_z_zzzz_0_1[i] * fe_0 + ta1_z_yzzzz_0_0[i] * pa_y[i] - ta1_z_yzzzz_0_1[i] * pc_y[i];

        ta1_z_yzzzzz_0_0[i] = ta1_z_zzzzz_0_0[i] * pa_y[i] - ta1_z_zzzzz_0_1[i] * pc_y[i];

        ta1_z_zzzzzz_0_0[i] = 5.0 * ta1_z_zzzz_0_0[i] * fe_0 - 5.0 * ta1_z_zzzz_0_1[i] * fe_0 + ta_zzzzz_0_1[i] + ta1_z_zzzzz_0_0[i] * pa_z[i] -
                              ta1_z_zzzzz_0_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
