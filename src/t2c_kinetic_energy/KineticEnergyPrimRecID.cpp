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

#include "KineticEnergyPrimRecID.hpp"

namespace kinrec {  // kinrec namespace

auto
comp_prim_kinetic_energy_id(CSimdArray<double>&       pbuffer,
                            const size_t              idx_kin_id,
                            const size_t              idx_ovl_gd,
                            const size_t              idx_kin_gd,
                            const size_t              idx_kin_hp,
                            const size_t              idx_kin_hd,
                            const size_t              idx_ovl_id,
                            const CSimdArray<double>& factors,
                            const size_t              idx_rpa,
                            const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GD

    auto ts_xxxx_xx = pbuffer.data(idx_ovl_gd);

    auto ts_xxxx_xy = pbuffer.data(idx_ovl_gd + 1);

    auto ts_xxxx_xz = pbuffer.data(idx_ovl_gd + 2);

    auto ts_xxxx_yy = pbuffer.data(idx_ovl_gd + 3);

    auto ts_xxxx_yz = pbuffer.data(idx_ovl_gd + 4);

    auto ts_xxxx_zz = pbuffer.data(idx_ovl_gd + 5);

    auto ts_xxxy_xx = pbuffer.data(idx_ovl_gd + 6);

    auto ts_xxxy_xz = pbuffer.data(idx_ovl_gd + 8);

    auto ts_xxxz_xx = pbuffer.data(idx_ovl_gd + 12);

    auto ts_xxxz_xy = pbuffer.data(idx_ovl_gd + 13);

    auto ts_xxyy_xx = pbuffer.data(idx_ovl_gd + 18);

    auto ts_xxyy_xy = pbuffer.data(idx_ovl_gd + 19);

    auto ts_xxyy_xz = pbuffer.data(idx_ovl_gd + 20);

    auto ts_xxyy_yy = pbuffer.data(idx_ovl_gd + 21);

    auto ts_xxyy_yz = pbuffer.data(idx_ovl_gd + 22);

    auto ts_xxyy_zz = pbuffer.data(idx_ovl_gd + 23);

    auto ts_xxzz_xx = pbuffer.data(idx_ovl_gd + 30);

    auto ts_xxzz_xy = pbuffer.data(idx_ovl_gd + 31);

    auto ts_xxzz_xz = pbuffer.data(idx_ovl_gd + 32);

    auto ts_xxzz_yy = pbuffer.data(idx_ovl_gd + 33);

    auto ts_xxzz_yz = pbuffer.data(idx_ovl_gd + 34);

    auto ts_xxzz_zz = pbuffer.data(idx_ovl_gd + 35);

    auto ts_xyyy_xy = pbuffer.data(idx_ovl_gd + 37);

    auto ts_xyyy_yy = pbuffer.data(idx_ovl_gd + 39);

    auto ts_xyyy_yz = pbuffer.data(idx_ovl_gd + 40);

    auto ts_xyyy_zz = pbuffer.data(idx_ovl_gd + 41);

    auto ts_xzzz_xz = pbuffer.data(idx_ovl_gd + 56);

    auto ts_xzzz_yy = pbuffer.data(idx_ovl_gd + 57);

    auto ts_xzzz_yz = pbuffer.data(idx_ovl_gd + 58);

    auto ts_xzzz_zz = pbuffer.data(idx_ovl_gd + 59);

    auto ts_yyyy_xx = pbuffer.data(idx_ovl_gd + 60);

    auto ts_yyyy_xy = pbuffer.data(idx_ovl_gd + 61);

    auto ts_yyyy_xz = pbuffer.data(idx_ovl_gd + 62);

    auto ts_yyyy_yy = pbuffer.data(idx_ovl_gd + 63);

    auto ts_yyyy_yz = pbuffer.data(idx_ovl_gd + 64);

    auto ts_yyyy_zz = pbuffer.data(idx_ovl_gd + 65);

    auto ts_yyyz_xy = pbuffer.data(idx_ovl_gd + 67);

    auto ts_yyyz_yy = pbuffer.data(idx_ovl_gd + 69);

    auto ts_yyzz_xx = pbuffer.data(idx_ovl_gd + 72);

    auto ts_yyzz_xy = pbuffer.data(idx_ovl_gd + 73);

    auto ts_yyzz_xz = pbuffer.data(idx_ovl_gd + 74);

    auto ts_yyzz_yy = pbuffer.data(idx_ovl_gd + 75);

    auto ts_yyzz_yz = pbuffer.data(idx_ovl_gd + 76);

    auto ts_yyzz_zz = pbuffer.data(idx_ovl_gd + 77);

    auto ts_yzzz_xx = pbuffer.data(idx_ovl_gd + 78);

    auto ts_yzzz_xz = pbuffer.data(idx_ovl_gd + 80);

    auto ts_yzzz_yz = pbuffer.data(idx_ovl_gd + 82);

    auto ts_yzzz_zz = pbuffer.data(idx_ovl_gd + 83);

    auto ts_zzzz_xx = pbuffer.data(idx_ovl_gd + 84);

    auto ts_zzzz_xy = pbuffer.data(idx_ovl_gd + 85);

    auto ts_zzzz_xz = pbuffer.data(idx_ovl_gd + 86);

    auto ts_zzzz_yy = pbuffer.data(idx_ovl_gd + 87);

    auto ts_zzzz_yz = pbuffer.data(idx_ovl_gd + 88);

    auto ts_zzzz_zz = pbuffer.data(idx_ovl_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto tk_xxxx_xx = pbuffer.data(idx_kin_gd);

    auto tk_xxxx_xy = pbuffer.data(idx_kin_gd + 1);

    auto tk_xxxx_xz = pbuffer.data(idx_kin_gd + 2);

    auto tk_xxxx_yy = pbuffer.data(idx_kin_gd + 3);

    auto tk_xxxx_yz = pbuffer.data(idx_kin_gd + 4);

    auto tk_xxxx_zz = pbuffer.data(idx_kin_gd + 5);

    auto tk_xxxy_xx = pbuffer.data(idx_kin_gd + 6);

    auto tk_xxxy_xz = pbuffer.data(idx_kin_gd + 8);

    auto tk_xxxz_xx = pbuffer.data(idx_kin_gd + 12);

    auto tk_xxxz_xy = pbuffer.data(idx_kin_gd + 13);

    auto tk_xxyy_xx = pbuffer.data(idx_kin_gd + 18);

    auto tk_xxyy_xy = pbuffer.data(idx_kin_gd + 19);

    auto tk_xxyy_xz = pbuffer.data(idx_kin_gd + 20);

    auto tk_xxyy_yy = pbuffer.data(idx_kin_gd + 21);

    auto tk_xxyy_yz = pbuffer.data(idx_kin_gd + 22);

    auto tk_xxyy_zz = pbuffer.data(idx_kin_gd + 23);

    auto tk_xxzz_xx = pbuffer.data(idx_kin_gd + 30);

    auto tk_xxzz_xy = pbuffer.data(idx_kin_gd + 31);

    auto tk_xxzz_xz = pbuffer.data(idx_kin_gd + 32);

    auto tk_xxzz_yy = pbuffer.data(idx_kin_gd + 33);

    auto tk_xxzz_yz = pbuffer.data(idx_kin_gd + 34);

    auto tk_xxzz_zz = pbuffer.data(idx_kin_gd + 35);

    auto tk_xyyy_xy = pbuffer.data(idx_kin_gd + 37);

    auto tk_xyyy_yy = pbuffer.data(idx_kin_gd + 39);

    auto tk_xyyy_yz = pbuffer.data(idx_kin_gd + 40);

    auto tk_xyyy_zz = pbuffer.data(idx_kin_gd + 41);

    auto tk_xzzz_xz = pbuffer.data(idx_kin_gd + 56);

    auto tk_xzzz_yy = pbuffer.data(idx_kin_gd + 57);

    auto tk_xzzz_yz = pbuffer.data(idx_kin_gd + 58);

    auto tk_xzzz_zz = pbuffer.data(idx_kin_gd + 59);

    auto tk_yyyy_xx = pbuffer.data(idx_kin_gd + 60);

    auto tk_yyyy_xy = pbuffer.data(idx_kin_gd + 61);

    auto tk_yyyy_xz = pbuffer.data(idx_kin_gd + 62);

    auto tk_yyyy_yy = pbuffer.data(idx_kin_gd + 63);

    auto tk_yyyy_yz = pbuffer.data(idx_kin_gd + 64);

    auto tk_yyyy_zz = pbuffer.data(idx_kin_gd + 65);

    auto tk_yyyz_xy = pbuffer.data(idx_kin_gd + 67);

    auto tk_yyyz_yy = pbuffer.data(idx_kin_gd + 69);

    auto tk_yyzz_xx = pbuffer.data(idx_kin_gd + 72);

    auto tk_yyzz_xy = pbuffer.data(idx_kin_gd + 73);

    auto tk_yyzz_xz = pbuffer.data(idx_kin_gd + 74);

    auto tk_yyzz_yy = pbuffer.data(idx_kin_gd + 75);

    auto tk_yyzz_yz = pbuffer.data(idx_kin_gd + 76);

    auto tk_yyzz_zz = pbuffer.data(idx_kin_gd + 77);

    auto tk_yzzz_xx = pbuffer.data(idx_kin_gd + 78);

    auto tk_yzzz_xz = pbuffer.data(idx_kin_gd + 80);

    auto tk_yzzz_yz = pbuffer.data(idx_kin_gd + 82);

    auto tk_yzzz_zz = pbuffer.data(idx_kin_gd + 83);

    auto tk_zzzz_xx = pbuffer.data(idx_kin_gd + 84);

    auto tk_zzzz_xy = pbuffer.data(idx_kin_gd + 85);

    auto tk_zzzz_xz = pbuffer.data(idx_kin_gd + 86);

    auto tk_zzzz_yy = pbuffer.data(idx_kin_gd + 87);

    auto tk_zzzz_yz = pbuffer.data(idx_kin_gd + 88);

    auto tk_zzzz_zz = pbuffer.data(idx_kin_gd + 89);

    // Set up components of auxiliary buffer : HP

    auto tk_xxxxx_x = pbuffer.data(idx_kin_hp);

    auto tk_xxxxx_y = pbuffer.data(idx_kin_hp + 1);

    auto tk_xxxxx_z = pbuffer.data(idx_kin_hp + 2);

    auto tk_xxxxz_z = pbuffer.data(idx_kin_hp + 8);

    auto tk_xxxyy_x = pbuffer.data(idx_kin_hp + 9);

    auto tk_xxxyy_y = pbuffer.data(idx_kin_hp + 10);

    auto tk_xxxyy_z = pbuffer.data(idx_kin_hp + 11);

    auto tk_xxxzz_x = pbuffer.data(idx_kin_hp + 15);

    auto tk_xxxzz_y = pbuffer.data(idx_kin_hp + 16);

    auto tk_xxxzz_z = pbuffer.data(idx_kin_hp + 17);

    auto tk_xxyyy_x = pbuffer.data(idx_kin_hp + 18);

    auto tk_xxyyy_y = pbuffer.data(idx_kin_hp + 19);

    auto tk_xxyyy_z = pbuffer.data(idx_kin_hp + 20);

    auto tk_xxzzz_x = pbuffer.data(idx_kin_hp + 27);

    auto tk_xxzzz_y = pbuffer.data(idx_kin_hp + 28);

    auto tk_xxzzz_z = pbuffer.data(idx_kin_hp + 29);

    auto tk_xyyyy_y = pbuffer.data(idx_kin_hp + 31);

    auto tk_xzzzz_z = pbuffer.data(idx_kin_hp + 44);

    auto tk_yyyyy_x = pbuffer.data(idx_kin_hp + 45);

    auto tk_yyyyy_y = pbuffer.data(idx_kin_hp + 46);

    auto tk_yyyyy_z = pbuffer.data(idx_kin_hp + 47);

    auto tk_yyyyz_z = pbuffer.data(idx_kin_hp + 50);

    auto tk_yyyzz_x = pbuffer.data(idx_kin_hp + 51);

    auto tk_yyyzz_y = pbuffer.data(idx_kin_hp + 52);

    auto tk_yyyzz_z = pbuffer.data(idx_kin_hp + 53);

    auto tk_yyzzz_x = pbuffer.data(idx_kin_hp + 54);

    auto tk_yyzzz_y = pbuffer.data(idx_kin_hp + 55);

    auto tk_yyzzz_z = pbuffer.data(idx_kin_hp + 56);

    auto tk_yzzzz_y = pbuffer.data(idx_kin_hp + 58);

    auto tk_yzzzz_z = pbuffer.data(idx_kin_hp + 59);

    auto tk_zzzzz_x = pbuffer.data(idx_kin_hp + 60);

    auto tk_zzzzz_y = pbuffer.data(idx_kin_hp + 61);

    auto tk_zzzzz_z = pbuffer.data(idx_kin_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto tk_xxxxx_xx = pbuffer.data(idx_kin_hd);

    auto tk_xxxxx_xy = pbuffer.data(idx_kin_hd + 1);

    auto tk_xxxxx_xz = pbuffer.data(idx_kin_hd + 2);

    auto tk_xxxxx_yy = pbuffer.data(idx_kin_hd + 3);

    auto tk_xxxxx_yz = pbuffer.data(idx_kin_hd + 4);

    auto tk_xxxxx_zz = pbuffer.data(idx_kin_hd + 5);

    auto tk_xxxxy_xx = pbuffer.data(idx_kin_hd + 6);

    auto tk_xxxxy_xy = pbuffer.data(idx_kin_hd + 7);

    auto tk_xxxxy_xz = pbuffer.data(idx_kin_hd + 8);

    auto tk_xxxxy_yy = pbuffer.data(idx_kin_hd + 9);

    auto tk_xxxxz_xx = pbuffer.data(idx_kin_hd + 12);

    auto tk_xxxxz_xy = pbuffer.data(idx_kin_hd + 13);

    auto tk_xxxxz_xz = pbuffer.data(idx_kin_hd + 14);

    auto tk_xxxxz_yz = pbuffer.data(idx_kin_hd + 16);

    auto tk_xxxxz_zz = pbuffer.data(idx_kin_hd + 17);

    auto tk_xxxyy_xx = pbuffer.data(idx_kin_hd + 18);

    auto tk_xxxyy_xy = pbuffer.data(idx_kin_hd + 19);

    auto tk_xxxyy_xz = pbuffer.data(idx_kin_hd + 20);

    auto tk_xxxyy_yy = pbuffer.data(idx_kin_hd + 21);

    auto tk_xxxyy_yz = pbuffer.data(idx_kin_hd + 22);

    auto tk_xxxyy_zz = pbuffer.data(idx_kin_hd + 23);

    auto tk_xxxzz_xx = pbuffer.data(idx_kin_hd + 30);

    auto tk_xxxzz_xy = pbuffer.data(idx_kin_hd + 31);

    auto tk_xxxzz_xz = pbuffer.data(idx_kin_hd + 32);

    auto tk_xxxzz_yy = pbuffer.data(idx_kin_hd + 33);

    auto tk_xxxzz_yz = pbuffer.data(idx_kin_hd + 34);

    auto tk_xxxzz_zz = pbuffer.data(idx_kin_hd + 35);

    auto tk_xxyyy_xx = pbuffer.data(idx_kin_hd + 36);

    auto tk_xxyyy_xy = pbuffer.data(idx_kin_hd + 37);

    auto tk_xxyyy_xz = pbuffer.data(idx_kin_hd + 38);

    auto tk_xxyyy_yy = pbuffer.data(idx_kin_hd + 39);

    auto tk_xxyyy_yz = pbuffer.data(idx_kin_hd + 40);

    auto tk_xxyyy_zz = pbuffer.data(idx_kin_hd + 41);

    auto tk_xxyyz_xy = pbuffer.data(idx_kin_hd + 43);

    auto tk_xxyzz_xx = pbuffer.data(idx_kin_hd + 48);

    auto tk_xxyzz_xz = pbuffer.data(idx_kin_hd + 50);

    auto tk_xxzzz_xx = pbuffer.data(idx_kin_hd + 54);

    auto tk_xxzzz_xy = pbuffer.data(idx_kin_hd + 55);

    auto tk_xxzzz_xz = pbuffer.data(idx_kin_hd + 56);

    auto tk_xxzzz_yy = pbuffer.data(idx_kin_hd + 57);

    auto tk_xxzzz_yz = pbuffer.data(idx_kin_hd + 58);

    auto tk_xxzzz_zz = pbuffer.data(idx_kin_hd + 59);

    auto tk_xyyyy_xx = pbuffer.data(idx_kin_hd + 60);

    auto tk_xyyyy_xy = pbuffer.data(idx_kin_hd + 61);

    auto tk_xyyyy_yy = pbuffer.data(idx_kin_hd + 63);

    auto tk_xyyyy_yz = pbuffer.data(idx_kin_hd + 64);

    auto tk_xyyyy_zz = pbuffer.data(idx_kin_hd + 65);

    auto tk_xyyzz_yy = pbuffer.data(idx_kin_hd + 75);

    auto tk_xyyzz_yz = pbuffer.data(idx_kin_hd + 76);

    auto tk_xyyzz_zz = pbuffer.data(idx_kin_hd + 77);

    auto tk_xzzzz_xx = pbuffer.data(idx_kin_hd + 84);

    auto tk_xzzzz_xz = pbuffer.data(idx_kin_hd + 86);

    auto tk_xzzzz_yy = pbuffer.data(idx_kin_hd + 87);

    auto tk_xzzzz_yz = pbuffer.data(idx_kin_hd + 88);

    auto tk_xzzzz_zz = pbuffer.data(idx_kin_hd + 89);

    auto tk_yyyyy_xx = pbuffer.data(idx_kin_hd + 90);

    auto tk_yyyyy_xy = pbuffer.data(idx_kin_hd + 91);

    auto tk_yyyyy_xz = pbuffer.data(idx_kin_hd + 92);

    auto tk_yyyyy_yy = pbuffer.data(idx_kin_hd + 93);

    auto tk_yyyyy_yz = pbuffer.data(idx_kin_hd + 94);

    auto tk_yyyyy_zz = pbuffer.data(idx_kin_hd + 95);

    auto tk_yyyyz_xy = pbuffer.data(idx_kin_hd + 97);

    auto tk_yyyyz_xz = pbuffer.data(idx_kin_hd + 98);

    auto tk_yyyyz_yy = pbuffer.data(idx_kin_hd + 99);

    auto tk_yyyyz_yz = pbuffer.data(idx_kin_hd + 100);

    auto tk_yyyyz_zz = pbuffer.data(idx_kin_hd + 101);

    auto tk_yyyzz_xx = pbuffer.data(idx_kin_hd + 102);

    auto tk_yyyzz_xy = pbuffer.data(idx_kin_hd + 103);

    auto tk_yyyzz_xz = pbuffer.data(idx_kin_hd + 104);

    auto tk_yyyzz_yy = pbuffer.data(idx_kin_hd + 105);

    auto tk_yyyzz_yz = pbuffer.data(idx_kin_hd + 106);

    auto tk_yyyzz_zz = pbuffer.data(idx_kin_hd + 107);

    auto tk_yyzzz_xx = pbuffer.data(idx_kin_hd + 108);

    auto tk_yyzzz_xy = pbuffer.data(idx_kin_hd + 109);

    auto tk_yyzzz_xz = pbuffer.data(idx_kin_hd + 110);

    auto tk_yyzzz_yy = pbuffer.data(idx_kin_hd + 111);

    auto tk_yyzzz_yz = pbuffer.data(idx_kin_hd + 112);

    auto tk_yyzzz_zz = pbuffer.data(idx_kin_hd + 113);

    auto tk_yzzzz_xx = pbuffer.data(idx_kin_hd + 114);

    auto tk_yzzzz_xy = pbuffer.data(idx_kin_hd + 115);

    auto tk_yzzzz_xz = pbuffer.data(idx_kin_hd + 116);

    auto tk_yzzzz_yy = pbuffer.data(idx_kin_hd + 117);

    auto tk_yzzzz_yz = pbuffer.data(idx_kin_hd + 118);

    auto tk_yzzzz_zz = pbuffer.data(idx_kin_hd + 119);

    auto tk_zzzzz_xx = pbuffer.data(idx_kin_hd + 120);

    auto tk_zzzzz_xy = pbuffer.data(idx_kin_hd + 121);

    auto tk_zzzzz_xz = pbuffer.data(idx_kin_hd + 122);

    auto tk_zzzzz_yy = pbuffer.data(idx_kin_hd + 123);

    auto tk_zzzzz_yz = pbuffer.data(idx_kin_hd + 124);

    auto tk_zzzzz_zz = pbuffer.data(idx_kin_hd + 125);

    // Set up components of auxiliary buffer : ID

    auto ts_xxxxxx_xx = pbuffer.data(idx_ovl_id);

    auto ts_xxxxxx_xy = pbuffer.data(idx_ovl_id + 1);

    auto ts_xxxxxx_xz = pbuffer.data(idx_ovl_id + 2);

    auto ts_xxxxxx_yy = pbuffer.data(idx_ovl_id + 3);

    auto ts_xxxxxx_yz = pbuffer.data(idx_ovl_id + 4);

    auto ts_xxxxxx_zz = pbuffer.data(idx_ovl_id + 5);

    auto ts_xxxxxy_xx = pbuffer.data(idx_ovl_id + 6);

    auto ts_xxxxxy_xy = pbuffer.data(idx_ovl_id + 7);

    auto ts_xxxxxy_xz = pbuffer.data(idx_ovl_id + 8);

    auto ts_xxxxxy_yy = pbuffer.data(idx_ovl_id + 9);

    auto ts_xxxxxy_yz = pbuffer.data(idx_ovl_id + 10);

    auto ts_xxxxxy_zz = pbuffer.data(idx_ovl_id + 11);

    auto ts_xxxxxz_xx = pbuffer.data(idx_ovl_id + 12);

    auto ts_xxxxxz_xy = pbuffer.data(idx_ovl_id + 13);

    auto ts_xxxxxz_xz = pbuffer.data(idx_ovl_id + 14);

    auto ts_xxxxxz_yy = pbuffer.data(idx_ovl_id + 15);

    auto ts_xxxxxz_yz = pbuffer.data(idx_ovl_id + 16);

    auto ts_xxxxxz_zz = pbuffer.data(idx_ovl_id + 17);

    auto ts_xxxxyy_xx = pbuffer.data(idx_ovl_id + 18);

    auto ts_xxxxyy_xy = pbuffer.data(idx_ovl_id + 19);

    auto ts_xxxxyy_xz = pbuffer.data(idx_ovl_id + 20);

    auto ts_xxxxyy_yy = pbuffer.data(idx_ovl_id + 21);

    auto ts_xxxxyy_yz = pbuffer.data(idx_ovl_id + 22);

    auto ts_xxxxyy_zz = pbuffer.data(idx_ovl_id + 23);

    auto ts_xxxxyz_xx = pbuffer.data(idx_ovl_id + 24);

    auto ts_xxxxyz_xy = pbuffer.data(idx_ovl_id + 25);

    auto ts_xxxxyz_xz = pbuffer.data(idx_ovl_id + 26);

    auto ts_xxxxyz_yy = pbuffer.data(idx_ovl_id + 27);

    auto ts_xxxxyz_yz = pbuffer.data(idx_ovl_id + 28);

    auto ts_xxxxyz_zz = pbuffer.data(idx_ovl_id + 29);

    auto ts_xxxxzz_xx = pbuffer.data(idx_ovl_id + 30);

    auto ts_xxxxzz_xy = pbuffer.data(idx_ovl_id + 31);

    auto ts_xxxxzz_xz = pbuffer.data(idx_ovl_id + 32);

    auto ts_xxxxzz_yy = pbuffer.data(idx_ovl_id + 33);

    auto ts_xxxxzz_yz = pbuffer.data(idx_ovl_id + 34);

    auto ts_xxxxzz_zz = pbuffer.data(idx_ovl_id + 35);

    auto ts_xxxyyy_xx = pbuffer.data(idx_ovl_id + 36);

    auto ts_xxxyyy_xy = pbuffer.data(idx_ovl_id + 37);

    auto ts_xxxyyy_xz = pbuffer.data(idx_ovl_id + 38);

    auto ts_xxxyyy_yy = pbuffer.data(idx_ovl_id + 39);

    auto ts_xxxyyy_yz = pbuffer.data(idx_ovl_id + 40);

    auto ts_xxxyyy_zz = pbuffer.data(idx_ovl_id + 41);

    auto ts_xxxyyz_xx = pbuffer.data(idx_ovl_id + 42);

    auto ts_xxxyyz_xy = pbuffer.data(idx_ovl_id + 43);

    auto ts_xxxyyz_xz = pbuffer.data(idx_ovl_id + 44);

    auto ts_xxxyyz_yy = pbuffer.data(idx_ovl_id + 45);

    auto ts_xxxyyz_yz = pbuffer.data(idx_ovl_id + 46);

    auto ts_xxxyyz_zz = pbuffer.data(idx_ovl_id + 47);

    auto ts_xxxyzz_xx = pbuffer.data(idx_ovl_id + 48);

    auto ts_xxxyzz_xy = pbuffer.data(idx_ovl_id + 49);

    auto ts_xxxyzz_xz = pbuffer.data(idx_ovl_id + 50);

    auto ts_xxxyzz_yy = pbuffer.data(idx_ovl_id + 51);

    auto ts_xxxyzz_yz = pbuffer.data(idx_ovl_id + 52);

    auto ts_xxxyzz_zz = pbuffer.data(idx_ovl_id + 53);

    auto ts_xxxzzz_xx = pbuffer.data(idx_ovl_id + 54);

    auto ts_xxxzzz_xy = pbuffer.data(idx_ovl_id + 55);

    auto ts_xxxzzz_xz = pbuffer.data(idx_ovl_id + 56);

    auto ts_xxxzzz_yy = pbuffer.data(idx_ovl_id + 57);

    auto ts_xxxzzz_yz = pbuffer.data(idx_ovl_id + 58);

    auto ts_xxxzzz_zz = pbuffer.data(idx_ovl_id + 59);

    auto ts_xxyyyy_xx = pbuffer.data(idx_ovl_id + 60);

    auto ts_xxyyyy_xy = pbuffer.data(idx_ovl_id + 61);

    auto ts_xxyyyy_xz = pbuffer.data(idx_ovl_id + 62);

    auto ts_xxyyyy_yy = pbuffer.data(idx_ovl_id + 63);

    auto ts_xxyyyy_yz = pbuffer.data(idx_ovl_id + 64);

    auto ts_xxyyyy_zz = pbuffer.data(idx_ovl_id + 65);

    auto ts_xxyyyz_xx = pbuffer.data(idx_ovl_id + 66);

    auto ts_xxyyyz_xy = pbuffer.data(idx_ovl_id + 67);

    auto ts_xxyyyz_xz = pbuffer.data(idx_ovl_id + 68);

    auto ts_xxyyyz_yy = pbuffer.data(idx_ovl_id + 69);

    auto ts_xxyyyz_yz = pbuffer.data(idx_ovl_id + 70);

    auto ts_xxyyyz_zz = pbuffer.data(idx_ovl_id + 71);

    auto ts_xxyyzz_xx = pbuffer.data(idx_ovl_id + 72);

    auto ts_xxyyzz_xy = pbuffer.data(idx_ovl_id + 73);

    auto ts_xxyyzz_xz = pbuffer.data(idx_ovl_id + 74);

    auto ts_xxyyzz_yy = pbuffer.data(idx_ovl_id + 75);

    auto ts_xxyyzz_yz = pbuffer.data(idx_ovl_id + 76);

    auto ts_xxyyzz_zz = pbuffer.data(idx_ovl_id + 77);

    auto ts_xxyzzz_xx = pbuffer.data(idx_ovl_id + 78);

    auto ts_xxyzzz_xy = pbuffer.data(idx_ovl_id + 79);

    auto ts_xxyzzz_xz = pbuffer.data(idx_ovl_id + 80);

    auto ts_xxyzzz_yy = pbuffer.data(idx_ovl_id + 81);

    auto ts_xxyzzz_yz = pbuffer.data(idx_ovl_id + 82);

    auto ts_xxyzzz_zz = pbuffer.data(idx_ovl_id + 83);

    auto ts_xxzzzz_xx = pbuffer.data(idx_ovl_id + 84);

    auto ts_xxzzzz_xy = pbuffer.data(idx_ovl_id + 85);

    auto ts_xxzzzz_xz = pbuffer.data(idx_ovl_id + 86);

    auto ts_xxzzzz_yy = pbuffer.data(idx_ovl_id + 87);

    auto ts_xxzzzz_yz = pbuffer.data(idx_ovl_id + 88);

    auto ts_xxzzzz_zz = pbuffer.data(idx_ovl_id + 89);

    auto ts_xyyyyy_xx = pbuffer.data(idx_ovl_id + 90);

    auto ts_xyyyyy_xy = pbuffer.data(idx_ovl_id + 91);

    auto ts_xyyyyy_xz = pbuffer.data(idx_ovl_id + 92);

    auto ts_xyyyyy_yy = pbuffer.data(idx_ovl_id + 93);

    auto ts_xyyyyy_yz = pbuffer.data(idx_ovl_id + 94);

    auto ts_xyyyyy_zz = pbuffer.data(idx_ovl_id + 95);

    auto ts_xyyyyz_xx = pbuffer.data(idx_ovl_id + 96);

    auto ts_xyyyyz_xy = pbuffer.data(idx_ovl_id + 97);

    auto ts_xyyyyz_xz = pbuffer.data(idx_ovl_id + 98);

    auto ts_xyyyyz_yy = pbuffer.data(idx_ovl_id + 99);

    auto ts_xyyyyz_yz = pbuffer.data(idx_ovl_id + 100);

    auto ts_xyyyyz_zz = pbuffer.data(idx_ovl_id + 101);

    auto ts_xyyyzz_xx = pbuffer.data(idx_ovl_id + 102);

    auto ts_xyyyzz_xy = pbuffer.data(idx_ovl_id + 103);

    auto ts_xyyyzz_xz = pbuffer.data(idx_ovl_id + 104);

    auto ts_xyyyzz_yy = pbuffer.data(idx_ovl_id + 105);

    auto ts_xyyyzz_yz = pbuffer.data(idx_ovl_id + 106);

    auto ts_xyyyzz_zz = pbuffer.data(idx_ovl_id + 107);

    auto ts_xyyzzz_xx = pbuffer.data(idx_ovl_id + 108);

    auto ts_xyyzzz_xy = pbuffer.data(idx_ovl_id + 109);

    auto ts_xyyzzz_xz = pbuffer.data(idx_ovl_id + 110);

    auto ts_xyyzzz_yy = pbuffer.data(idx_ovl_id + 111);

    auto ts_xyyzzz_yz = pbuffer.data(idx_ovl_id + 112);

    auto ts_xyyzzz_zz = pbuffer.data(idx_ovl_id + 113);

    auto ts_xyzzzz_xx = pbuffer.data(idx_ovl_id + 114);

    auto ts_xyzzzz_xy = pbuffer.data(idx_ovl_id + 115);

    auto ts_xyzzzz_xz = pbuffer.data(idx_ovl_id + 116);

    auto ts_xyzzzz_yy = pbuffer.data(idx_ovl_id + 117);

    auto ts_xyzzzz_yz = pbuffer.data(idx_ovl_id + 118);

    auto ts_xyzzzz_zz = pbuffer.data(idx_ovl_id + 119);

    auto ts_xzzzzz_xx = pbuffer.data(idx_ovl_id + 120);

    auto ts_xzzzzz_xy = pbuffer.data(idx_ovl_id + 121);

    auto ts_xzzzzz_xz = pbuffer.data(idx_ovl_id + 122);

    auto ts_xzzzzz_yy = pbuffer.data(idx_ovl_id + 123);

    auto ts_xzzzzz_yz = pbuffer.data(idx_ovl_id + 124);

    auto ts_xzzzzz_zz = pbuffer.data(idx_ovl_id + 125);

    auto ts_yyyyyy_xx = pbuffer.data(idx_ovl_id + 126);

    auto ts_yyyyyy_xy = pbuffer.data(idx_ovl_id + 127);

    auto ts_yyyyyy_xz = pbuffer.data(idx_ovl_id + 128);

    auto ts_yyyyyy_yy = pbuffer.data(idx_ovl_id + 129);

    auto ts_yyyyyy_yz = pbuffer.data(idx_ovl_id + 130);

    auto ts_yyyyyy_zz = pbuffer.data(idx_ovl_id + 131);

    auto ts_yyyyyz_xx = pbuffer.data(idx_ovl_id + 132);

    auto ts_yyyyyz_xy = pbuffer.data(idx_ovl_id + 133);

    auto ts_yyyyyz_xz = pbuffer.data(idx_ovl_id + 134);

    auto ts_yyyyyz_yy = pbuffer.data(idx_ovl_id + 135);

    auto ts_yyyyyz_yz = pbuffer.data(idx_ovl_id + 136);

    auto ts_yyyyyz_zz = pbuffer.data(idx_ovl_id + 137);

    auto ts_yyyyzz_xx = pbuffer.data(idx_ovl_id + 138);

    auto ts_yyyyzz_xy = pbuffer.data(idx_ovl_id + 139);

    auto ts_yyyyzz_xz = pbuffer.data(idx_ovl_id + 140);

    auto ts_yyyyzz_yy = pbuffer.data(idx_ovl_id + 141);

    auto ts_yyyyzz_yz = pbuffer.data(idx_ovl_id + 142);

    auto ts_yyyyzz_zz = pbuffer.data(idx_ovl_id + 143);

    auto ts_yyyzzz_xx = pbuffer.data(idx_ovl_id + 144);

    auto ts_yyyzzz_xy = pbuffer.data(idx_ovl_id + 145);

    auto ts_yyyzzz_xz = pbuffer.data(idx_ovl_id + 146);

    auto ts_yyyzzz_yy = pbuffer.data(idx_ovl_id + 147);

    auto ts_yyyzzz_yz = pbuffer.data(idx_ovl_id + 148);

    auto ts_yyyzzz_zz = pbuffer.data(idx_ovl_id + 149);

    auto ts_yyzzzz_xx = pbuffer.data(idx_ovl_id + 150);

    auto ts_yyzzzz_xy = pbuffer.data(idx_ovl_id + 151);

    auto ts_yyzzzz_xz = pbuffer.data(idx_ovl_id + 152);

    auto ts_yyzzzz_yy = pbuffer.data(idx_ovl_id + 153);

    auto ts_yyzzzz_yz = pbuffer.data(idx_ovl_id + 154);

    auto ts_yyzzzz_zz = pbuffer.data(idx_ovl_id + 155);

    auto ts_yzzzzz_xx = pbuffer.data(idx_ovl_id + 156);

    auto ts_yzzzzz_xy = pbuffer.data(idx_ovl_id + 157);

    auto ts_yzzzzz_xz = pbuffer.data(idx_ovl_id + 158);

    auto ts_yzzzzz_yy = pbuffer.data(idx_ovl_id + 159);

    auto ts_yzzzzz_yz = pbuffer.data(idx_ovl_id + 160);

    auto ts_yzzzzz_zz = pbuffer.data(idx_ovl_id + 161);

    auto ts_zzzzzz_xx = pbuffer.data(idx_ovl_id + 162);

    auto ts_zzzzzz_xy = pbuffer.data(idx_ovl_id + 163);

    auto ts_zzzzzz_xz = pbuffer.data(idx_ovl_id + 164);

    auto ts_zzzzzz_yy = pbuffer.data(idx_ovl_id + 165);

    auto ts_zzzzzz_yz = pbuffer.data(idx_ovl_id + 166);

    auto ts_zzzzzz_zz = pbuffer.data(idx_ovl_id + 167);

    // Set up 0-6 components of targeted buffer : ID

    auto tk_xxxxxx_xx = pbuffer.data(idx_kin_id);

    auto tk_xxxxxx_xy = pbuffer.data(idx_kin_id + 1);

    auto tk_xxxxxx_xz = pbuffer.data(idx_kin_id + 2);

    auto tk_xxxxxx_yy = pbuffer.data(idx_kin_id + 3);

    auto tk_xxxxxx_yz = pbuffer.data(idx_kin_id + 4);

    auto tk_xxxxxx_zz = pbuffer.data(idx_kin_id + 5);

#pragma omp simd aligned(pa_x,             \
                             tk_xxxx_xx,   \
                             tk_xxxx_xy,   \
                             tk_xxxx_xz,   \
                             tk_xxxx_yy,   \
                             tk_xxxx_yz,   \
                             tk_xxxx_zz,   \
                             tk_xxxxx_x,   \
                             tk_xxxxx_xx,  \
                             tk_xxxxx_xy,  \
                             tk_xxxxx_xz,  \
                             tk_xxxxx_y,   \
                             tk_xxxxx_yy,  \
                             tk_xxxxx_yz,  \
                             tk_xxxxx_z,   \
                             tk_xxxxx_zz,  \
                             tk_xxxxxx_xx, \
                             tk_xxxxxx_xy, \
                             tk_xxxxxx_xz, \
                             tk_xxxxxx_yy, \
                             tk_xxxxxx_yz, \
                             tk_xxxxxx_zz, \
                             ts_xxxx_xx,   \
                             ts_xxxx_xy,   \
                             ts_xxxx_xz,   \
                             ts_xxxx_yy,   \
                             ts_xxxx_yz,   \
                             ts_xxxx_zz,   \
                             ts_xxxxxx_xx, \
                             ts_xxxxxx_xy, \
                             ts_xxxxxx_xz, \
                             ts_xxxxxx_yy, \
                             ts_xxxxxx_yz, \
                             ts_xxxxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxxx_xx[i] = -10.0 * ts_xxxx_xx[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xx[i] * fe_0 + 2.0 * tk_xxxxx_x[i] * fe_0 + tk_xxxxx_xx[i] * pa_x[i] +
                          2.0 * ts_xxxxxx_xx[i] * fz_0;

        tk_xxxxxx_xy[i] = -10.0 * ts_xxxx_xy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xy[i] * fe_0 + tk_xxxxx_y[i] * fe_0 + tk_xxxxx_xy[i] * pa_x[i] +
                          2.0 * ts_xxxxxx_xy[i] * fz_0;

        tk_xxxxxx_xz[i] = -10.0 * ts_xxxx_xz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_xz[i] * fe_0 + tk_xxxxx_z[i] * fe_0 + tk_xxxxx_xz[i] * pa_x[i] +
                          2.0 * ts_xxxxxx_xz[i] * fz_0;

        tk_xxxxxx_yy[i] = -10.0 * ts_xxxx_yy[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yy[i] * fe_0 + tk_xxxxx_yy[i] * pa_x[i] + 2.0 * ts_xxxxxx_yy[i] * fz_0;

        tk_xxxxxx_yz[i] = -10.0 * ts_xxxx_yz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_yz[i] * fe_0 + tk_xxxxx_yz[i] * pa_x[i] + 2.0 * ts_xxxxxx_yz[i] * fz_0;

        tk_xxxxxx_zz[i] = -10.0 * ts_xxxx_zz[i] * fbe_0 * fz_0 + 5.0 * tk_xxxx_zz[i] * fe_0 + tk_xxxxx_zz[i] * pa_x[i] + 2.0 * ts_xxxxxx_zz[i] * fz_0;
    }

    // Set up 6-12 components of targeted buffer : ID

    auto tk_xxxxxy_xx = pbuffer.data(idx_kin_id + 6);

    auto tk_xxxxxy_xy = pbuffer.data(idx_kin_id + 7);

    auto tk_xxxxxy_xz = pbuffer.data(idx_kin_id + 8);

    auto tk_xxxxxy_yy = pbuffer.data(idx_kin_id + 9);

    auto tk_xxxxxy_yz = pbuffer.data(idx_kin_id + 10);

    auto tk_xxxxxy_zz = pbuffer.data(idx_kin_id + 11);

#pragma omp simd aligned(pa_y,             \
                             tk_xxxxx_x,   \
                             tk_xxxxx_xx,  \
                             tk_xxxxx_xy,  \
                             tk_xxxxx_xz,  \
                             tk_xxxxx_y,   \
                             tk_xxxxx_yy,  \
                             tk_xxxxx_yz,  \
                             tk_xxxxx_z,   \
                             tk_xxxxx_zz,  \
                             tk_xxxxxy_xx, \
                             tk_xxxxxy_xy, \
                             tk_xxxxxy_xz, \
                             tk_xxxxxy_yy, \
                             tk_xxxxxy_yz, \
                             tk_xxxxxy_zz, \
                             ts_xxxxxy_xx, \
                             ts_xxxxxy_xy, \
                             ts_xxxxxy_xz, \
                             ts_xxxxxy_yy, \
                             ts_xxxxxy_yz, \
                             ts_xxxxxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxy_xx[i] = tk_xxxxx_xx[i] * pa_y[i] + 2.0 * ts_xxxxxy_xx[i] * fz_0;

        tk_xxxxxy_xy[i] = tk_xxxxx_x[i] * fe_0 + tk_xxxxx_xy[i] * pa_y[i] + 2.0 * ts_xxxxxy_xy[i] * fz_0;

        tk_xxxxxy_xz[i] = tk_xxxxx_xz[i] * pa_y[i] + 2.0 * ts_xxxxxy_xz[i] * fz_0;

        tk_xxxxxy_yy[i] = 2.0 * tk_xxxxx_y[i] * fe_0 + tk_xxxxx_yy[i] * pa_y[i] + 2.0 * ts_xxxxxy_yy[i] * fz_0;

        tk_xxxxxy_yz[i] = tk_xxxxx_z[i] * fe_0 + tk_xxxxx_yz[i] * pa_y[i] + 2.0 * ts_xxxxxy_yz[i] * fz_0;

        tk_xxxxxy_zz[i] = tk_xxxxx_zz[i] * pa_y[i] + 2.0 * ts_xxxxxy_zz[i] * fz_0;
    }

    // Set up 12-18 components of targeted buffer : ID

    auto tk_xxxxxz_xx = pbuffer.data(idx_kin_id + 12);

    auto tk_xxxxxz_xy = pbuffer.data(idx_kin_id + 13);

    auto tk_xxxxxz_xz = pbuffer.data(idx_kin_id + 14);

    auto tk_xxxxxz_yy = pbuffer.data(idx_kin_id + 15);

    auto tk_xxxxxz_yz = pbuffer.data(idx_kin_id + 16);

    auto tk_xxxxxz_zz = pbuffer.data(idx_kin_id + 17);

#pragma omp simd aligned(pa_z,             \
                             tk_xxxxx_x,   \
                             tk_xxxxx_xx,  \
                             tk_xxxxx_xy,  \
                             tk_xxxxx_xz,  \
                             tk_xxxxx_y,   \
                             tk_xxxxx_yy,  \
                             tk_xxxxx_yz,  \
                             tk_xxxxx_z,   \
                             tk_xxxxx_zz,  \
                             tk_xxxxxz_xx, \
                             tk_xxxxxz_xy, \
                             tk_xxxxxz_xz, \
                             tk_xxxxxz_yy, \
                             tk_xxxxxz_yz, \
                             tk_xxxxxz_zz, \
                             ts_xxxxxz_xx, \
                             ts_xxxxxz_xy, \
                             ts_xxxxxz_xz, \
                             ts_xxxxxz_yy, \
                             ts_xxxxxz_yz, \
                             ts_xxxxxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxxz_xx[i] = tk_xxxxx_xx[i] * pa_z[i] + 2.0 * ts_xxxxxz_xx[i] * fz_0;

        tk_xxxxxz_xy[i] = tk_xxxxx_xy[i] * pa_z[i] + 2.0 * ts_xxxxxz_xy[i] * fz_0;

        tk_xxxxxz_xz[i] = tk_xxxxx_x[i] * fe_0 + tk_xxxxx_xz[i] * pa_z[i] + 2.0 * ts_xxxxxz_xz[i] * fz_0;

        tk_xxxxxz_yy[i] = tk_xxxxx_yy[i] * pa_z[i] + 2.0 * ts_xxxxxz_yy[i] * fz_0;

        tk_xxxxxz_yz[i] = tk_xxxxx_y[i] * fe_0 + tk_xxxxx_yz[i] * pa_z[i] + 2.0 * ts_xxxxxz_yz[i] * fz_0;

        tk_xxxxxz_zz[i] = 2.0 * tk_xxxxx_z[i] * fe_0 + tk_xxxxx_zz[i] * pa_z[i] + 2.0 * ts_xxxxxz_zz[i] * fz_0;
    }

    // Set up 18-24 components of targeted buffer : ID

    auto tk_xxxxyy_xx = pbuffer.data(idx_kin_id + 18);

    auto tk_xxxxyy_xy = pbuffer.data(idx_kin_id + 19);

    auto tk_xxxxyy_xz = pbuffer.data(idx_kin_id + 20);

    auto tk_xxxxyy_yy = pbuffer.data(idx_kin_id + 21);

    auto tk_xxxxyy_yz = pbuffer.data(idx_kin_id + 22);

    auto tk_xxxxyy_zz = pbuffer.data(idx_kin_id + 23);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xxxx_xx,   \
                             tk_xxxx_xz,   \
                             tk_xxxxy_xx,  \
                             tk_xxxxy_xz,  \
                             tk_xxxxyy_xx, \
                             tk_xxxxyy_xy, \
                             tk_xxxxyy_xz, \
                             tk_xxxxyy_yy, \
                             tk_xxxxyy_yz, \
                             tk_xxxxyy_zz, \
                             tk_xxxyy_xy,  \
                             tk_xxxyy_y,   \
                             tk_xxxyy_yy,  \
                             tk_xxxyy_yz,  \
                             tk_xxxyy_zz,  \
                             tk_xxyy_xy,   \
                             tk_xxyy_yy,   \
                             tk_xxyy_yz,   \
                             tk_xxyy_zz,   \
                             ts_xxxx_xx,   \
                             ts_xxxx_xz,   \
                             ts_xxxxyy_xx, \
                             ts_xxxxyy_xy, \
                             ts_xxxxyy_xz, \
                             ts_xxxxyy_yy, \
                             ts_xxxxyy_yz, \
                             ts_xxxxyy_zz, \
                             ts_xxyy_xy,   \
                             ts_xxyy_yy,   \
                             ts_xxyy_yz,   \
                             ts_xxyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxyy_xx[i] = -2.0 * ts_xxxx_xx[i] * fbe_0 * fz_0 + tk_xxxx_xx[i] * fe_0 + tk_xxxxy_xx[i] * pa_y[i] + 2.0 * ts_xxxxyy_xx[i] * fz_0;

        tk_xxxxyy_xy[i] = -6.0 * ts_xxyy_xy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xy[i] * fe_0 + tk_xxxyy_y[i] * fe_0 + tk_xxxyy_xy[i] * pa_x[i] +
                          2.0 * ts_xxxxyy_xy[i] * fz_0;

        tk_xxxxyy_xz[i] = -2.0 * ts_xxxx_xz[i] * fbe_0 * fz_0 + tk_xxxx_xz[i] * fe_0 + tk_xxxxy_xz[i] * pa_y[i] + 2.0 * ts_xxxxyy_xz[i] * fz_0;

        tk_xxxxyy_yy[i] = -6.0 * ts_xxyy_yy[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yy[i] * fe_0 + tk_xxxyy_yy[i] * pa_x[i] + 2.0 * ts_xxxxyy_yy[i] * fz_0;

        tk_xxxxyy_yz[i] = -6.0 * ts_xxyy_yz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_yz[i] * fe_0 + tk_xxxyy_yz[i] * pa_x[i] + 2.0 * ts_xxxxyy_yz[i] * fz_0;

        tk_xxxxyy_zz[i] = -6.0 * ts_xxyy_zz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_zz[i] * fe_0 + tk_xxxyy_zz[i] * pa_x[i] + 2.0 * ts_xxxxyy_zz[i] * fz_0;
    }

    // Set up 24-30 components of targeted buffer : ID

    auto tk_xxxxyz_xx = pbuffer.data(idx_kin_id + 24);

    auto tk_xxxxyz_xy = pbuffer.data(idx_kin_id + 25);

    auto tk_xxxxyz_xz = pbuffer.data(idx_kin_id + 26);

    auto tk_xxxxyz_yy = pbuffer.data(idx_kin_id + 27);

    auto tk_xxxxyz_yz = pbuffer.data(idx_kin_id + 28);

    auto tk_xxxxyz_zz = pbuffer.data(idx_kin_id + 29);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_xxxxy_xy,  \
                             tk_xxxxy_yy,  \
                             tk_xxxxyz_xx, \
                             tk_xxxxyz_xy, \
                             tk_xxxxyz_xz, \
                             tk_xxxxyz_yy, \
                             tk_xxxxyz_yz, \
                             tk_xxxxyz_zz, \
                             tk_xxxxz_xx,  \
                             tk_xxxxz_xz,  \
                             tk_xxxxz_yz,  \
                             tk_xxxxz_z,   \
                             tk_xxxxz_zz,  \
                             ts_xxxxyz_xx, \
                             ts_xxxxyz_xy, \
                             ts_xxxxyz_xz, \
                             ts_xxxxyz_yy, \
                             ts_xxxxyz_yz, \
                             ts_xxxxyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxxyz_xx[i] = tk_xxxxz_xx[i] * pa_y[i] + 2.0 * ts_xxxxyz_xx[i] * fz_0;

        tk_xxxxyz_xy[i] = tk_xxxxy_xy[i] * pa_z[i] + 2.0 * ts_xxxxyz_xy[i] * fz_0;

        tk_xxxxyz_xz[i] = tk_xxxxz_xz[i] * pa_y[i] + 2.0 * ts_xxxxyz_xz[i] * fz_0;

        tk_xxxxyz_yy[i] = tk_xxxxy_yy[i] * pa_z[i] + 2.0 * ts_xxxxyz_yy[i] * fz_0;

        tk_xxxxyz_yz[i] = tk_xxxxz_z[i] * fe_0 + tk_xxxxz_yz[i] * pa_y[i] + 2.0 * ts_xxxxyz_yz[i] * fz_0;

        tk_xxxxyz_zz[i] = tk_xxxxz_zz[i] * pa_y[i] + 2.0 * ts_xxxxyz_zz[i] * fz_0;
    }

    // Set up 30-36 components of targeted buffer : ID

    auto tk_xxxxzz_xx = pbuffer.data(idx_kin_id + 30);

    auto tk_xxxxzz_xy = pbuffer.data(idx_kin_id + 31);

    auto tk_xxxxzz_xz = pbuffer.data(idx_kin_id + 32);

    auto tk_xxxxzz_yy = pbuffer.data(idx_kin_id + 33);

    auto tk_xxxxzz_yz = pbuffer.data(idx_kin_id + 34);

    auto tk_xxxxzz_zz = pbuffer.data(idx_kin_id + 35);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xxxx_xx,   \
                             tk_xxxx_xy,   \
                             tk_xxxxz_xx,  \
                             tk_xxxxz_xy,  \
                             tk_xxxxzz_xx, \
                             tk_xxxxzz_xy, \
                             tk_xxxxzz_xz, \
                             tk_xxxxzz_yy, \
                             tk_xxxxzz_yz, \
                             tk_xxxxzz_zz, \
                             tk_xxxzz_xz,  \
                             tk_xxxzz_yy,  \
                             tk_xxxzz_yz,  \
                             tk_xxxzz_z,   \
                             tk_xxxzz_zz,  \
                             tk_xxzz_xz,   \
                             tk_xxzz_yy,   \
                             tk_xxzz_yz,   \
                             tk_xxzz_zz,   \
                             ts_xxxx_xx,   \
                             ts_xxxx_xy,   \
                             ts_xxxxzz_xx, \
                             ts_xxxxzz_xy, \
                             ts_xxxxzz_xz, \
                             ts_xxxxzz_yy, \
                             ts_xxxxzz_yz, \
                             ts_xxxxzz_zz, \
                             ts_xxzz_xz,   \
                             ts_xxzz_yy,   \
                             ts_xxzz_yz,   \
                             ts_xxzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxxzz_xx[i] = -2.0 * ts_xxxx_xx[i] * fbe_0 * fz_0 + tk_xxxx_xx[i] * fe_0 + tk_xxxxz_xx[i] * pa_z[i] + 2.0 * ts_xxxxzz_xx[i] * fz_0;

        tk_xxxxzz_xy[i] = -2.0 * ts_xxxx_xy[i] * fbe_0 * fz_0 + tk_xxxx_xy[i] * fe_0 + tk_xxxxz_xy[i] * pa_z[i] + 2.0 * ts_xxxxzz_xy[i] * fz_0;

        tk_xxxxzz_xz[i] = -6.0 * ts_xxzz_xz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xz[i] * fe_0 + tk_xxxzz_z[i] * fe_0 + tk_xxxzz_xz[i] * pa_x[i] +
                          2.0 * ts_xxxxzz_xz[i] * fz_0;

        tk_xxxxzz_yy[i] = -6.0 * ts_xxzz_yy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yy[i] * fe_0 + tk_xxxzz_yy[i] * pa_x[i] + 2.0 * ts_xxxxzz_yy[i] * fz_0;

        tk_xxxxzz_yz[i] = -6.0 * ts_xxzz_yz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_yz[i] * fe_0 + tk_xxxzz_yz[i] * pa_x[i] + 2.0 * ts_xxxxzz_yz[i] * fz_0;

        tk_xxxxzz_zz[i] = -6.0 * ts_xxzz_zz[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_zz[i] * fe_0 + tk_xxxzz_zz[i] * pa_x[i] + 2.0 * ts_xxxxzz_zz[i] * fz_0;
    }

    // Set up 36-42 components of targeted buffer : ID

    auto tk_xxxyyy_xx = pbuffer.data(idx_kin_id + 36);

    auto tk_xxxyyy_xy = pbuffer.data(idx_kin_id + 37);

    auto tk_xxxyyy_xz = pbuffer.data(idx_kin_id + 38);

    auto tk_xxxyyy_yy = pbuffer.data(idx_kin_id + 39);

    auto tk_xxxyyy_yz = pbuffer.data(idx_kin_id + 40);

    auto tk_xxxyyy_zz = pbuffer.data(idx_kin_id + 41);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xxxy_xx,   \
                             tk_xxxy_xz,   \
                             tk_xxxyy_xx,  \
                             tk_xxxyy_xz,  \
                             tk_xxxyyy_xx, \
                             tk_xxxyyy_xy, \
                             tk_xxxyyy_xz, \
                             tk_xxxyyy_yy, \
                             tk_xxxyyy_yz, \
                             tk_xxxyyy_zz, \
                             tk_xxyyy_xy,  \
                             tk_xxyyy_y,   \
                             tk_xxyyy_yy,  \
                             tk_xxyyy_yz,  \
                             tk_xxyyy_zz,  \
                             tk_xyyy_xy,   \
                             tk_xyyy_yy,   \
                             tk_xyyy_yz,   \
                             tk_xyyy_zz,   \
                             ts_xxxy_xx,   \
                             ts_xxxy_xz,   \
                             ts_xxxyyy_xx, \
                             ts_xxxyyy_xy, \
                             ts_xxxyyy_xz, \
                             ts_xxxyyy_yy, \
                             ts_xxxyyy_yz, \
                             ts_xxxyyy_zz, \
                             ts_xyyy_xy,   \
                             ts_xyyy_yy,   \
                             ts_xyyy_yz,   \
                             ts_xyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxyyy_xx[i] = -4.0 * ts_xxxy_xx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xx[i] * fe_0 + tk_xxxyy_xx[i] * pa_y[i] + 2.0 * ts_xxxyyy_xx[i] * fz_0;

        tk_xxxyyy_xy[i] = -4.0 * ts_xyyy_xy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_xy[i] * fe_0 + tk_xxyyy_y[i] * fe_0 + tk_xxyyy_xy[i] * pa_x[i] +
                          2.0 * ts_xxxyyy_xy[i] * fz_0;

        tk_xxxyyy_xz[i] = -4.0 * ts_xxxy_xz[i] * fbe_0 * fz_0 + 2.0 * tk_xxxy_xz[i] * fe_0 + tk_xxxyy_xz[i] * pa_y[i] + 2.0 * ts_xxxyyy_xz[i] * fz_0;

        tk_xxxyyy_yy[i] = -4.0 * ts_xyyy_yy[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yy[i] * fe_0 + tk_xxyyy_yy[i] * pa_x[i] + 2.0 * ts_xxxyyy_yy[i] * fz_0;

        tk_xxxyyy_yz[i] = -4.0 * ts_xyyy_yz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_yz[i] * fe_0 + tk_xxyyy_yz[i] * pa_x[i] + 2.0 * ts_xxxyyy_yz[i] * fz_0;

        tk_xxxyyy_zz[i] = -4.0 * ts_xyyy_zz[i] * fbe_0 * fz_0 + 2.0 * tk_xyyy_zz[i] * fe_0 + tk_xxyyy_zz[i] * pa_x[i] + 2.0 * ts_xxxyyy_zz[i] * fz_0;
    }

    // Set up 42-48 components of targeted buffer : ID

    auto tk_xxxyyz_xx = pbuffer.data(idx_kin_id + 42);

    auto tk_xxxyyz_xy = pbuffer.data(idx_kin_id + 43);

    auto tk_xxxyyz_xz = pbuffer.data(idx_kin_id + 44);

    auto tk_xxxyyz_yy = pbuffer.data(idx_kin_id + 45);

    auto tk_xxxyyz_yz = pbuffer.data(idx_kin_id + 46);

    auto tk_xxxyyz_zz = pbuffer.data(idx_kin_id + 47);

#pragma omp simd aligned(pa_z,             \
                             tk_xxxyy_x,   \
                             tk_xxxyy_xx,  \
                             tk_xxxyy_xy,  \
                             tk_xxxyy_xz,  \
                             tk_xxxyy_y,   \
                             tk_xxxyy_yy,  \
                             tk_xxxyy_yz,  \
                             tk_xxxyy_z,   \
                             tk_xxxyy_zz,  \
                             tk_xxxyyz_xx, \
                             tk_xxxyyz_xy, \
                             tk_xxxyyz_xz, \
                             tk_xxxyyz_yy, \
                             tk_xxxyyz_yz, \
                             tk_xxxyyz_zz, \
                             ts_xxxyyz_xx, \
                             ts_xxxyyz_xy, \
                             ts_xxxyyz_xz, \
                             ts_xxxyyz_yy, \
                             ts_xxxyyz_yz, \
                             ts_xxxyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyyz_xx[i] = tk_xxxyy_xx[i] * pa_z[i] + 2.0 * ts_xxxyyz_xx[i] * fz_0;

        tk_xxxyyz_xy[i] = tk_xxxyy_xy[i] * pa_z[i] + 2.0 * ts_xxxyyz_xy[i] * fz_0;

        tk_xxxyyz_xz[i] = tk_xxxyy_x[i] * fe_0 + tk_xxxyy_xz[i] * pa_z[i] + 2.0 * ts_xxxyyz_xz[i] * fz_0;

        tk_xxxyyz_yy[i] = tk_xxxyy_yy[i] * pa_z[i] + 2.0 * ts_xxxyyz_yy[i] * fz_0;

        tk_xxxyyz_yz[i] = tk_xxxyy_y[i] * fe_0 + tk_xxxyy_yz[i] * pa_z[i] + 2.0 * ts_xxxyyz_yz[i] * fz_0;

        tk_xxxyyz_zz[i] = 2.0 * tk_xxxyy_z[i] * fe_0 + tk_xxxyy_zz[i] * pa_z[i] + 2.0 * ts_xxxyyz_zz[i] * fz_0;
    }

    // Set up 48-54 components of targeted buffer : ID

    auto tk_xxxyzz_xx = pbuffer.data(idx_kin_id + 48);

    auto tk_xxxyzz_xy = pbuffer.data(idx_kin_id + 49);

    auto tk_xxxyzz_xz = pbuffer.data(idx_kin_id + 50);

    auto tk_xxxyzz_yy = pbuffer.data(idx_kin_id + 51);

    auto tk_xxxyzz_yz = pbuffer.data(idx_kin_id + 52);

    auto tk_xxxyzz_zz = pbuffer.data(idx_kin_id + 53);

#pragma omp simd aligned(pa_y,             \
                             tk_xxxyzz_xx, \
                             tk_xxxyzz_xy, \
                             tk_xxxyzz_xz, \
                             tk_xxxyzz_yy, \
                             tk_xxxyzz_yz, \
                             tk_xxxyzz_zz, \
                             tk_xxxzz_x,   \
                             tk_xxxzz_xx,  \
                             tk_xxxzz_xy,  \
                             tk_xxxzz_xz,  \
                             tk_xxxzz_y,   \
                             tk_xxxzz_yy,  \
                             tk_xxxzz_yz,  \
                             tk_xxxzz_z,   \
                             tk_xxxzz_zz,  \
                             ts_xxxyzz_xx, \
                             ts_xxxyzz_xy, \
                             ts_xxxyzz_xz, \
                             ts_xxxyzz_yy, \
                             ts_xxxyzz_yz, \
                             ts_xxxyzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxxyzz_xx[i] = tk_xxxzz_xx[i] * pa_y[i] + 2.0 * ts_xxxyzz_xx[i] * fz_0;

        tk_xxxyzz_xy[i] = tk_xxxzz_x[i] * fe_0 + tk_xxxzz_xy[i] * pa_y[i] + 2.0 * ts_xxxyzz_xy[i] * fz_0;

        tk_xxxyzz_xz[i] = tk_xxxzz_xz[i] * pa_y[i] + 2.0 * ts_xxxyzz_xz[i] * fz_0;

        tk_xxxyzz_yy[i] = 2.0 * tk_xxxzz_y[i] * fe_0 + tk_xxxzz_yy[i] * pa_y[i] + 2.0 * ts_xxxyzz_yy[i] * fz_0;

        tk_xxxyzz_yz[i] = tk_xxxzz_z[i] * fe_0 + tk_xxxzz_yz[i] * pa_y[i] + 2.0 * ts_xxxyzz_yz[i] * fz_0;

        tk_xxxyzz_zz[i] = tk_xxxzz_zz[i] * pa_y[i] + 2.0 * ts_xxxyzz_zz[i] * fz_0;
    }

    // Set up 54-60 components of targeted buffer : ID

    auto tk_xxxzzz_xx = pbuffer.data(idx_kin_id + 54);

    auto tk_xxxzzz_xy = pbuffer.data(idx_kin_id + 55);

    auto tk_xxxzzz_xz = pbuffer.data(idx_kin_id + 56);

    auto tk_xxxzzz_yy = pbuffer.data(idx_kin_id + 57);

    auto tk_xxxzzz_yz = pbuffer.data(idx_kin_id + 58);

    auto tk_xxxzzz_zz = pbuffer.data(idx_kin_id + 59);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xxxz_xx,   \
                             tk_xxxz_xy,   \
                             tk_xxxzz_xx,  \
                             tk_xxxzz_xy,  \
                             tk_xxxzzz_xx, \
                             tk_xxxzzz_xy, \
                             tk_xxxzzz_xz, \
                             tk_xxxzzz_yy, \
                             tk_xxxzzz_yz, \
                             tk_xxxzzz_zz, \
                             tk_xxzzz_xz,  \
                             tk_xxzzz_yy,  \
                             tk_xxzzz_yz,  \
                             tk_xxzzz_z,   \
                             tk_xxzzz_zz,  \
                             tk_xzzz_xz,   \
                             tk_xzzz_yy,   \
                             tk_xzzz_yz,   \
                             tk_xzzz_zz,   \
                             ts_xxxz_xx,   \
                             ts_xxxz_xy,   \
                             ts_xxxzzz_xx, \
                             ts_xxxzzz_xy, \
                             ts_xxxzzz_xz, \
                             ts_xxxzzz_yy, \
                             ts_xxxzzz_yz, \
                             ts_xxxzzz_zz, \
                             ts_xzzz_xz,   \
                             ts_xzzz_yy,   \
                             ts_xzzz_yz,   \
                             ts_xzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxxzzz_xx[i] = -4.0 * ts_xxxz_xx[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xx[i] * fe_0 + tk_xxxzz_xx[i] * pa_z[i] + 2.0 * ts_xxxzzz_xx[i] * fz_0;

        tk_xxxzzz_xy[i] = -4.0 * ts_xxxz_xy[i] * fbe_0 * fz_0 + 2.0 * tk_xxxz_xy[i] * fe_0 + tk_xxxzz_xy[i] * pa_z[i] + 2.0 * ts_xxxzzz_xy[i] * fz_0;

        tk_xxxzzz_xz[i] = -4.0 * ts_xzzz_xz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_xz[i] * fe_0 + tk_xxzzz_z[i] * fe_0 + tk_xxzzz_xz[i] * pa_x[i] +
                          2.0 * ts_xxxzzz_xz[i] * fz_0;

        tk_xxxzzz_yy[i] = -4.0 * ts_xzzz_yy[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yy[i] * fe_0 + tk_xxzzz_yy[i] * pa_x[i] + 2.0 * ts_xxxzzz_yy[i] * fz_0;

        tk_xxxzzz_yz[i] = -4.0 * ts_xzzz_yz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_yz[i] * fe_0 + tk_xxzzz_yz[i] * pa_x[i] + 2.0 * ts_xxxzzz_yz[i] * fz_0;

        tk_xxxzzz_zz[i] = -4.0 * ts_xzzz_zz[i] * fbe_0 * fz_0 + 2.0 * tk_xzzz_zz[i] * fe_0 + tk_xxzzz_zz[i] * pa_x[i] + 2.0 * ts_xxxzzz_zz[i] * fz_0;
    }

    // Set up 60-66 components of targeted buffer : ID

    auto tk_xxyyyy_xx = pbuffer.data(idx_kin_id + 60);

    auto tk_xxyyyy_xy = pbuffer.data(idx_kin_id + 61);

    auto tk_xxyyyy_xz = pbuffer.data(idx_kin_id + 62);

    auto tk_xxyyyy_yy = pbuffer.data(idx_kin_id + 63);

    auto tk_xxyyyy_yz = pbuffer.data(idx_kin_id + 64);

    auto tk_xxyyyy_zz = pbuffer.data(idx_kin_id + 65);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xxyy_xx,   \
                             tk_xxyy_xz,   \
                             tk_xxyyy_xx,  \
                             tk_xxyyy_xz,  \
                             tk_xxyyyy_xx, \
                             tk_xxyyyy_xy, \
                             tk_xxyyyy_xz, \
                             tk_xxyyyy_yy, \
                             tk_xxyyyy_yz, \
                             tk_xxyyyy_zz, \
                             tk_xyyyy_xy,  \
                             tk_xyyyy_y,   \
                             tk_xyyyy_yy,  \
                             tk_xyyyy_yz,  \
                             tk_xyyyy_zz,  \
                             tk_yyyy_xy,   \
                             tk_yyyy_yy,   \
                             tk_yyyy_yz,   \
                             tk_yyyy_zz,   \
                             ts_xxyy_xx,   \
                             ts_xxyy_xz,   \
                             ts_xxyyyy_xx, \
                             ts_xxyyyy_xy, \
                             ts_xxyyyy_xz, \
                             ts_xxyyyy_yy, \
                             ts_xxyyyy_yz, \
                             ts_xxyyyy_zz, \
                             ts_yyyy_xy,   \
                             ts_yyyy_yy,   \
                             ts_yyyy_yz,   \
                             ts_yyyy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyyy_xx[i] = -6.0 * ts_xxyy_xx[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xx[i] * fe_0 + tk_xxyyy_xx[i] * pa_y[i] + 2.0 * ts_xxyyyy_xx[i] * fz_0;

        tk_xxyyyy_xy[i] = -2.0 * ts_yyyy_xy[i] * fbe_0 * fz_0 + tk_yyyy_xy[i] * fe_0 + tk_xyyyy_y[i] * fe_0 + tk_xyyyy_xy[i] * pa_x[i] +
                          2.0 * ts_xxyyyy_xy[i] * fz_0;

        tk_xxyyyy_xz[i] = -6.0 * ts_xxyy_xz[i] * fbe_0 * fz_0 + 3.0 * tk_xxyy_xz[i] * fe_0 + tk_xxyyy_xz[i] * pa_y[i] + 2.0 * ts_xxyyyy_xz[i] * fz_0;

        tk_xxyyyy_yy[i] = -2.0 * ts_yyyy_yy[i] * fbe_0 * fz_0 + tk_yyyy_yy[i] * fe_0 + tk_xyyyy_yy[i] * pa_x[i] + 2.0 * ts_xxyyyy_yy[i] * fz_0;

        tk_xxyyyy_yz[i] = -2.0 * ts_yyyy_yz[i] * fbe_0 * fz_0 + tk_yyyy_yz[i] * fe_0 + tk_xyyyy_yz[i] * pa_x[i] + 2.0 * ts_xxyyyy_yz[i] * fz_0;

        tk_xxyyyy_zz[i] = -2.0 * ts_yyyy_zz[i] * fbe_0 * fz_0 + tk_yyyy_zz[i] * fe_0 + tk_xyyyy_zz[i] * pa_x[i] + 2.0 * ts_xxyyyy_zz[i] * fz_0;
    }

    // Set up 66-72 components of targeted buffer : ID

    auto tk_xxyyyz_xx = pbuffer.data(idx_kin_id + 66);

    auto tk_xxyyyz_xy = pbuffer.data(idx_kin_id + 67);

    auto tk_xxyyyz_xz = pbuffer.data(idx_kin_id + 68);

    auto tk_xxyyyz_yy = pbuffer.data(idx_kin_id + 69);

    auto tk_xxyyyz_yz = pbuffer.data(idx_kin_id + 70);

    auto tk_xxyyyz_zz = pbuffer.data(idx_kin_id + 71);

#pragma omp simd aligned(pa_z,             \
                             tk_xxyyy_x,   \
                             tk_xxyyy_xx,  \
                             tk_xxyyy_xy,  \
                             tk_xxyyy_xz,  \
                             tk_xxyyy_y,   \
                             tk_xxyyy_yy,  \
                             tk_xxyyy_yz,  \
                             tk_xxyyy_z,   \
                             tk_xxyyy_zz,  \
                             tk_xxyyyz_xx, \
                             tk_xxyyyz_xy, \
                             tk_xxyyyz_xz, \
                             tk_xxyyyz_yy, \
                             tk_xxyyyz_yz, \
                             tk_xxyyyz_zz, \
                             ts_xxyyyz_xx, \
                             ts_xxyyyz_xy, \
                             ts_xxyyyz_xz, \
                             ts_xxyyyz_yy, \
                             ts_xxyyyz_yz, \
                             ts_xxyyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyyyz_xx[i] = tk_xxyyy_xx[i] * pa_z[i] + 2.0 * ts_xxyyyz_xx[i] * fz_0;

        tk_xxyyyz_xy[i] = tk_xxyyy_xy[i] * pa_z[i] + 2.0 * ts_xxyyyz_xy[i] * fz_0;

        tk_xxyyyz_xz[i] = tk_xxyyy_x[i] * fe_0 + tk_xxyyy_xz[i] * pa_z[i] + 2.0 * ts_xxyyyz_xz[i] * fz_0;

        tk_xxyyyz_yy[i] = tk_xxyyy_yy[i] * pa_z[i] + 2.0 * ts_xxyyyz_yy[i] * fz_0;

        tk_xxyyyz_yz[i] = tk_xxyyy_y[i] * fe_0 + tk_xxyyy_yz[i] * pa_z[i] + 2.0 * ts_xxyyyz_yz[i] * fz_0;

        tk_xxyyyz_zz[i] = 2.0 * tk_xxyyy_z[i] * fe_0 + tk_xxyyy_zz[i] * pa_z[i] + 2.0 * ts_xxyyyz_zz[i] * fz_0;
    }

    // Set up 72-78 components of targeted buffer : ID

    auto tk_xxyyzz_xx = pbuffer.data(idx_kin_id + 72);

    auto tk_xxyyzz_xy = pbuffer.data(idx_kin_id + 73);

    auto tk_xxyyzz_xz = pbuffer.data(idx_kin_id + 74);

    auto tk_xxyyzz_yy = pbuffer.data(idx_kin_id + 75);

    auto tk_xxyyzz_yz = pbuffer.data(idx_kin_id + 76);

    auto tk_xxyyzz_zz = pbuffer.data(idx_kin_id + 77);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             tk_xxyy_xy,   \
                             tk_xxyyz_xy,  \
                             tk_xxyyzz_xx, \
                             tk_xxyyzz_xy, \
                             tk_xxyyzz_xz, \
                             tk_xxyyzz_yy, \
                             tk_xxyyzz_yz, \
                             tk_xxyyzz_zz, \
                             tk_xxyzz_xx,  \
                             tk_xxyzz_xz,  \
                             tk_xxzz_xx,   \
                             tk_xxzz_xz,   \
                             tk_xyyzz_yy,  \
                             tk_xyyzz_yz,  \
                             tk_xyyzz_zz,  \
                             tk_yyzz_yy,   \
                             tk_yyzz_yz,   \
                             tk_yyzz_zz,   \
                             ts_xxyy_xy,   \
                             ts_xxyyzz_xx, \
                             ts_xxyyzz_xy, \
                             ts_xxyyzz_xz, \
                             ts_xxyyzz_yy, \
                             ts_xxyyzz_yz, \
                             ts_xxyyzz_zz, \
                             ts_xxzz_xx,   \
                             ts_xxzz_xz,   \
                             ts_yyzz_yy,   \
                             ts_yyzz_yz,   \
                             ts_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxyyzz_xx[i] = -2.0 * ts_xxzz_xx[i] * fbe_0 * fz_0 + tk_xxzz_xx[i] * fe_0 + tk_xxyzz_xx[i] * pa_y[i] + 2.0 * ts_xxyyzz_xx[i] * fz_0;

        tk_xxyyzz_xy[i] = -2.0 * ts_xxyy_xy[i] * fbe_0 * fz_0 + tk_xxyy_xy[i] * fe_0 + tk_xxyyz_xy[i] * pa_z[i] + 2.0 * ts_xxyyzz_xy[i] * fz_0;

        tk_xxyyzz_xz[i] = -2.0 * ts_xxzz_xz[i] * fbe_0 * fz_0 + tk_xxzz_xz[i] * fe_0 + tk_xxyzz_xz[i] * pa_y[i] + 2.0 * ts_xxyyzz_xz[i] * fz_0;

        tk_xxyyzz_yy[i] = -2.0 * ts_yyzz_yy[i] * fbe_0 * fz_0 + tk_yyzz_yy[i] * fe_0 + tk_xyyzz_yy[i] * pa_x[i] + 2.0 * ts_xxyyzz_yy[i] * fz_0;

        tk_xxyyzz_yz[i] = -2.0 * ts_yyzz_yz[i] * fbe_0 * fz_0 + tk_yyzz_yz[i] * fe_0 + tk_xyyzz_yz[i] * pa_x[i] + 2.0 * ts_xxyyzz_yz[i] * fz_0;

        tk_xxyyzz_zz[i] = -2.0 * ts_yyzz_zz[i] * fbe_0 * fz_0 + tk_yyzz_zz[i] * fe_0 + tk_xyyzz_zz[i] * pa_x[i] + 2.0 * ts_xxyyzz_zz[i] * fz_0;
    }

    // Set up 78-84 components of targeted buffer : ID

    auto tk_xxyzzz_xx = pbuffer.data(idx_kin_id + 78);

    auto tk_xxyzzz_xy = pbuffer.data(idx_kin_id + 79);

    auto tk_xxyzzz_xz = pbuffer.data(idx_kin_id + 80);

    auto tk_xxyzzz_yy = pbuffer.data(idx_kin_id + 81);

    auto tk_xxyzzz_yz = pbuffer.data(idx_kin_id + 82);

    auto tk_xxyzzz_zz = pbuffer.data(idx_kin_id + 83);

#pragma omp simd aligned(pa_y,             \
                             tk_xxyzzz_xx, \
                             tk_xxyzzz_xy, \
                             tk_xxyzzz_xz, \
                             tk_xxyzzz_yy, \
                             tk_xxyzzz_yz, \
                             tk_xxyzzz_zz, \
                             tk_xxzzz_x,   \
                             tk_xxzzz_xx,  \
                             tk_xxzzz_xy,  \
                             tk_xxzzz_xz,  \
                             tk_xxzzz_y,   \
                             tk_xxzzz_yy,  \
                             tk_xxzzz_yz,  \
                             tk_xxzzz_z,   \
                             tk_xxzzz_zz,  \
                             ts_xxyzzz_xx, \
                             ts_xxyzzz_xy, \
                             ts_xxyzzz_xz, \
                             ts_xxyzzz_yy, \
                             ts_xxyzzz_yz, \
                             ts_xxyzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xxyzzz_xx[i] = tk_xxzzz_xx[i] * pa_y[i] + 2.0 * ts_xxyzzz_xx[i] * fz_0;

        tk_xxyzzz_xy[i] = tk_xxzzz_x[i] * fe_0 + tk_xxzzz_xy[i] * pa_y[i] + 2.0 * ts_xxyzzz_xy[i] * fz_0;

        tk_xxyzzz_xz[i] = tk_xxzzz_xz[i] * pa_y[i] + 2.0 * ts_xxyzzz_xz[i] * fz_0;

        tk_xxyzzz_yy[i] = 2.0 * tk_xxzzz_y[i] * fe_0 + tk_xxzzz_yy[i] * pa_y[i] + 2.0 * ts_xxyzzz_yy[i] * fz_0;

        tk_xxyzzz_yz[i] = tk_xxzzz_z[i] * fe_0 + tk_xxzzz_yz[i] * pa_y[i] + 2.0 * ts_xxyzzz_yz[i] * fz_0;

        tk_xxyzzz_zz[i] = tk_xxzzz_zz[i] * pa_y[i] + 2.0 * ts_xxyzzz_zz[i] * fz_0;
    }

    // Set up 84-90 components of targeted buffer : ID

    auto tk_xxzzzz_xx = pbuffer.data(idx_kin_id + 84);

    auto tk_xxzzzz_xy = pbuffer.data(idx_kin_id + 85);

    auto tk_xxzzzz_xz = pbuffer.data(idx_kin_id + 86);

    auto tk_xxzzzz_yy = pbuffer.data(idx_kin_id + 87);

    auto tk_xxzzzz_yz = pbuffer.data(idx_kin_id + 88);

    auto tk_xxzzzz_zz = pbuffer.data(idx_kin_id + 89);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xxzz_xx,   \
                             tk_xxzz_xy,   \
                             tk_xxzzz_xx,  \
                             tk_xxzzz_xy,  \
                             tk_xxzzzz_xx, \
                             tk_xxzzzz_xy, \
                             tk_xxzzzz_xz, \
                             tk_xxzzzz_yy, \
                             tk_xxzzzz_yz, \
                             tk_xxzzzz_zz, \
                             tk_xzzzz_xz,  \
                             tk_xzzzz_yy,  \
                             tk_xzzzz_yz,  \
                             tk_xzzzz_z,   \
                             tk_xzzzz_zz,  \
                             tk_zzzz_xz,   \
                             tk_zzzz_yy,   \
                             tk_zzzz_yz,   \
                             tk_zzzz_zz,   \
                             ts_xxzz_xx,   \
                             ts_xxzz_xy,   \
                             ts_xxzzzz_xx, \
                             ts_xxzzzz_xy, \
                             ts_xxzzzz_xz, \
                             ts_xxzzzz_yy, \
                             ts_xxzzzz_yz, \
                             ts_xxzzzz_zz, \
                             ts_zzzz_xz,   \
                             ts_zzzz_yy,   \
                             ts_zzzz_yz,   \
                             ts_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_xxzzzz_xx[i] = -6.0 * ts_xxzz_xx[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xx[i] * fe_0 + tk_xxzzz_xx[i] * pa_z[i] + 2.0 * ts_xxzzzz_xx[i] * fz_0;

        tk_xxzzzz_xy[i] = -6.0 * ts_xxzz_xy[i] * fbe_0 * fz_0 + 3.0 * tk_xxzz_xy[i] * fe_0 + tk_xxzzz_xy[i] * pa_z[i] + 2.0 * ts_xxzzzz_xy[i] * fz_0;

        tk_xxzzzz_xz[i] = -2.0 * ts_zzzz_xz[i] * fbe_0 * fz_0 + tk_zzzz_xz[i] * fe_0 + tk_xzzzz_z[i] * fe_0 + tk_xzzzz_xz[i] * pa_x[i] +
                          2.0 * ts_xxzzzz_xz[i] * fz_0;

        tk_xxzzzz_yy[i] = -2.0 * ts_zzzz_yy[i] * fbe_0 * fz_0 + tk_zzzz_yy[i] * fe_0 + tk_xzzzz_yy[i] * pa_x[i] + 2.0 * ts_xxzzzz_yy[i] * fz_0;

        tk_xxzzzz_yz[i] = -2.0 * ts_zzzz_yz[i] * fbe_0 * fz_0 + tk_zzzz_yz[i] * fe_0 + tk_xzzzz_yz[i] * pa_x[i] + 2.0 * ts_xxzzzz_yz[i] * fz_0;

        tk_xxzzzz_zz[i] = -2.0 * ts_zzzz_zz[i] * fbe_0 * fz_0 + tk_zzzz_zz[i] * fe_0 + tk_xzzzz_zz[i] * pa_x[i] + 2.0 * ts_xxzzzz_zz[i] * fz_0;
    }

    // Set up 90-96 components of targeted buffer : ID

    auto tk_xyyyyy_xx = pbuffer.data(idx_kin_id + 90);

    auto tk_xyyyyy_xy = pbuffer.data(idx_kin_id + 91);

    auto tk_xyyyyy_xz = pbuffer.data(idx_kin_id + 92);

    auto tk_xyyyyy_yy = pbuffer.data(idx_kin_id + 93);

    auto tk_xyyyyy_yz = pbuffer.data(idx_kin_id + 94);

    auto tk_xyyyyy_zz = pbuffer.data(idx_kin_id + 95);

#pragma omp simd aligned(pa_x,             \
                             tk_xyyyyy_xx, \
                             tk_xyyyyy_xy, \
                             tk_xyyyyy_xz, \
                             tk_xyyyyy_yy, \
                             tk_xyyyyy_yz, \
                             tk_xyyyyy_zz, \
                             tk_yyyyy_x,   \
                             tk_yyyyy_xx,  \
                             tk_yyyyy_xy,  \
                             tk_yyyyy_xz,  \
                             tk_yyyyy_y,   \
                             tk_yyyyy_yy,  \
                             tk_yyyyy_yz,  \
                             tk_yyyyy_z,   \
                             tk_yyyyy_zz,  \
                             ts_xyyyyy_xx, \
                             ts_xyyyyy_xy, \
                             ts_xyyyyy_xz, \
                             ts_xyyyyy_yy, \
                             ts_xyyyyy_yz, \
                             ts_xyyyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyy_xx[i] = 2.0 * tk_yyyyy_x[i] * fe_0 + tk_yyyyy_xx[i] * pa_x[i] + 2.0 * ts_xyyyyy_xx[i] * fz_0;

        tk_xyyyyy_xy[i] = tk_yyyyy_y[i] * fe_0 + tk_yyyyy_xy[i] * pa_x[i] + 2.0 * ts_xyyyyy_xy[i] * fz_0;

        tk_xyyyyy_xz[i] = tk_yyyyy_z[i] * fe_0 + tk_yyyyy_xz[i] * pa_x[i] + 2.0 * ts_xyyyyy_xz[i] * fz_0;

        tk_xyyyyy_yy[i] = tk_yyyyy_yy[i] * pa_x[i] + 2.0 * ts_xyyyyy_yy[i] * fz_0;

        tk_xyyyyy_yz[i] = tk_yyyyy_yz[i] * pa_x[i] + 2.0 * ts_xyyyyy_yz[i] * fz_0;

        tk_xyyyyy_zz[i] = tk_yyyyy_zz[i] * pa_x[i] + 2.0 * ts_xyyyyy_zz[i] * fz_0;
    }

    // Set up 96-102 components of targeted buffer : ID

    auto tk_xyyyyz_xx = pbuffer.data(idx_kin_id + 96);

    auto tk_xyyyyz_xy = pbuffer.data(idx_kin_id + 97);

    auto tk_xyyyyz_xz = pbuffer.data(idx_kin_id + 98);

    auto tk_xyyyyz_yy = pbuffer.data(idx_kin_id + 99);

    auto tk_xyyyyz_yz = pbuffer.data(idx_kin_id + 100);

    auto tk_xyyyyz_zz = pbuffer.data(idx_kin_id + 101);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tk_xyyyy_xx,  \
                             tk_xyyyy_xy,  \
                             tk_xyyyyz_xx, \
                             tk_xyyyyz_xy, \
                             tk_xyyyyz_xz, \
                             tk_xyyyyz_yy, \
                             tk_xyyyyz_yz, \
                             tk_xyyyyz_zz, \
                             tk_yyyyz_xz,  \
                             tk_yyyyz_yy,  \
                             tk_yyyyz_yz,  \
                             tk_yyyyz_z,   \
                             tk_yyyyz_zz,  \
                             ts_xyyyyz_xx, \
                             ts_xyyyyz_xy, \
                             ts_xyyyyz_xz, \
                             ts_xyyyyz_yy, \
                             ts_xyyyyz_yz, \
                             ts_xyyyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyyz_xx[i] = tk_xyyyy_xx[i] * pa_z[i] + 2.0 * ts_xyyyyz_xx[i] * fz_0;

        tk_xyyyyz_xy[i] = tk_xyyyy_xy[i] * pa_z[i] + 2.0 * ts_xyyyyz_xy[i] * fz_0;

        tk_xyyyyz_xz[i] = tk_yyyyz_z[i] * fe_0 + tk_yyyyz_xz[i] * pa_x[i] + 2.0 * ts_xyyyyz_xz[i] * fz_0;

        tk_xyyyyz_yy[i] = tk_yyyyz_yy[i] * pa_x[i] + 2.0 * ts_xyyyyz_yy[i] * fz_0;

        tk_xyyyyz_yz[i] = tk_yyyyz_yz[i] * pa_x[i] + 2.0 * ts_xyyyyz_yz[i] * fz_0;

        tk_xyyyyz_zz[i] = tk_yyyyz_zz[i] * pa_x[i] + 2.0 * ts_xyyyyz_zz[i] * fz_0;
    }

    // Set up 102-108 components of targeted buffer : ID

    auto tk_xyyyzz_xx = pbuffer.data(idx_kin_id + 102);

    auto tk_xyyyzz_xy = pbuffer.data(idx_kin_id + 103);

    auto tk_xyyyzz_xz = pbuffer.data(idx_kin_id + 104);

    auto tk_xyyyzz_yy = pbuffer.data(idx_kin_id + 105);

    auto tk_xyyyzz_yz = pbuffer.data(idx_kin_id + 106);

    auto tk_xyyyzz_zz = pbuffer.data(idx_kin_id + 107);

#pragma omp simd aligned(pa_x,             \
                             tk_xyyyzz_xx, \
                             tk_xyyyzz_xy, \
                             tk_xyyyzz_xz, \
                             tk_xyyyzz_yy, \
                             tk_xyyyzz_yz, \
                             tk_xyyyzz_zz, \
                             tk_yyyzz_x,   \
                             tk_yyyzz_xx,  \
                             tk_yyyzz_xy,  \
                             tk_yyyzz_xz,  \
                             tk_yyyzz_y,   \
                             tk_yyyzz_yy,  \
                             tk_yyyzz_yz,  \
                             tk_yyyzz_z,   \
                             tk_yyyzz_zz,  \
                             ts_xyyyzz_xx, \
                             ts_xyyyzz_xy, \
                             ts_xyyyzz_xz, \
                             ts_xyyyzz_yy, \
                             ts_xyyyzz_yz, \
                             ts_xyyyzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyyzz_xx[i] = 2.0 * tk_yyyzz_x[i] * fe_0 + tk_yyyzz_xx[i] * pa_x[i] + 2.0 * ts_xyyyzz_xx[i] * fz_0;

        tk_xyyyzz_xy[i] = tk_yyyzz_y[i] * fe_0 + tk_yyyzz_xy[i] * pa_x[i] + 2.0 * ts_xyyyzz_xy[i] * fz_0;

        tk_xyyyzz_xz[i] = tk_yyyzz_z[i] * fe_0 + tk_yyyzz_xz[i] * pa_x[i] + 2.0 * ts_xyyyzz_xz[i] * fz_0;

        tk_xyyyzz_yy[i] = tk_yyyzz_yy[i] * pa_x[i] + 2.0 * ts_xyyyzz_yy[i] * fz_0;

        tk_xyyyzz_yz[i] = tk_yyyzz_yz[i] * pa_x[i] + 2.0 * ts_xyyyzz_yz[i] * fz_0;

        tk_xyyyzz_zz[i] = tk_yyyzz_zz[i] * pa_x[i] + 2.0 * ts_xyyyzz_zz[i] * fz_0;
    }

    // Set up 108-114 components of targeted buffer : ID

    auto tk_xyyzzz_xx = pbuffer.data(idx_kin_id + 108);

    auto tk_xyyzzz_xy = pbuffer.data(idx_kin_id + 109);

    auto tk_xyyzzz_xz = pbuffer.data(idx_kin_id + 110);

    auto tk_xyyzzz_yy = pbuffer.data(idx_kin_id + 111);

    auto tk_xyyzzz_yz = pbuffer.data(idx_kin_id + 112);

    auto tk_xyyzzz_zz = pbuffer.data(idx_kin_id + 113);

#pragma omp simd aligned(pa_x,             \
                             tk_xyyzzz_xx, \
                             tk_xyyzzz_xy, \
                             tk_xyyzzz_xz, \
                             tk_xyyzzz_yy, \
                             tk_xyyzzz_yz, \
                             tk_xyyzzz_zz, \
                             tk_yyzzz_x,   \
                             tk_yyzzz_xx,  \
                             tk_yyzzz_xy,  \
                             tk_yyzzz_xz,  \
                             tk_yyzzz_y,   \
                             tk_yyzzz_yy,  \
                             tk_yyzzz_yz,  \
                             tk_yyzzz_z,   \
                             tk_yyzzz_zz,  \
                             ts_xyyzzz_xx, \
                             ts_xyyzzz_xy, \
                             ts_xyyzzz_xz, \
                             ts_xyyzzz_yy, \
                             ts_xyyzzz_yz, \
                             ts_xyyzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyyzzz_xx[i] = 2.0 * tk_yyzzz_x[i] * fe_0 + tk_yyzzz_xx[i] * pa_x[i] + 2.0 * ts_xyyzzz_xx[i] * fz_0;

        tk_xyyzzz_xy[i] = tk_yyzzz_y[i] * fe_0 + tk_yyzzz_xy[i] * pa_x[i] + 2.0 * ts_xyyzzz_xy[i] * fz_0;

        tk_xyyzzz_xz[i] = tk_yyzzz_z[i] * fe_0 + tk_yyzzz_xz[i] * pa_x[i] + 2.0 * ts_xyyzzz_xz[i] * fz_0;

        tk_xyyzzz_yy[i] = tk_yyzzz_yy[i] * pa_x[i] + 2.0 * ts_xyyzzz_yy[i] * fz_0;

        tk_xyyzzz_yz[i] = tk_yyzzz_yz[i] * pa_x[i] + 2.0 * ts_xyyzzz_yz[i] * fz_0;

        tk_xyyzzz_zz[i] = tk_yyzzz_zz[i] * pa_x[i] + 2.0 * ts_xyyzzz_zz[i] * fz_0;
    }

    // Set up 114-120 components of targeted buffer : ID

    auto tk_xyzzzz_xx = pbuffer.data(idx_kin_id + 114);

    auto tk_xyzzzz_xy = pbuffer.data(idx_kin_id + 115);

    auto tk_xyzzzz_xz = pbuffer.data(idx_kin_id + 116);

    auto tk_xyzzzz_yy = pbuffer.data(idx_kin_id + 117);

    auto tk_xyzzzz_yz = pbuffer.data(idx_kin_id + 118);

    auto tk_xyzzzz_zz = pbuffer.data(idx_kin_id + 119);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tk_xyzzzz_xx, \
                             tk_xyzzzz_xy, \
                             tk_xyzzzz_xz, \
                             tk_xyzzzz_yy, \
                             tk_xyzzzz_yz, \
                             tk_xyzzzz_zz, \
                             tk_xzzzz_xx,  \
                             tk_xzzzz_xz,  \
                             tk_yzzzz_xy,  \
                             tk_yzzzz_y,   \
                             tk_yzzzz_yy,  \
                             tk_yzzzz_yz,  \
                             tk_yzzzz_zz,  \
                             ts_xyzzzz_xx, \
                             ts_xyzzzz_xy, \
                             ts_xyzzzz_xz, \
                             ts_xyzzzz_yy, \
                             ts_xyzzzz_yz, \
                             ts_xyzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xyzzzz_xx[i] = tk_xzzzz_xx[i] * pa_y[i] + 2.0 * ts_xyzzzz_xx[i] * fz_0;

        tk_xyzzzz_xy[i] = tk_yzzzz_y[i] * fe_0 + tk_yzzzz_xy[i] * pa_x[i] + 2.0 * ts_xyzzzz_xy[i] * fz_0;

        tk_xyzzzz_xz[i] = tk_xzzzz_xz[i] * pa_y[i] + 2.0 * ts_xyzzzz_xz[i] * fz_0;

        tk_xyzzzz_yy[i] = tk_yzzzz_yy[i] * pa_x[i] + 2.0 * ts_xyzzzz_yy[i] * fz_0;

        tk_xyzzzz_yz[i] = tk_yzzzz_yz[i] * pa_x[i] + 2.0 * ts_xyzzzz_yz[i] * fz_0;

        tk_xyzzzz_zz[i] = tk_yzzzz_zz[i] * pa_x[i] + 2.0 * ts_xyzzzz_zz[i] * fz_0;
    }

    // Set up 120-126 components of targeted buffer : ID

    auto tk_xzzzzz_xx = pbuffer.data(idx_kin_id + 120);

    auto tk_xzzzzz_xy = pbuffer.data(idx_kin_id + 121);

    auto tk_xzzzzz_xz = pbuffer.data(idx_kin_id + 122);

    auto tk_xzzzzz_yy = pbuffer.data(idx_kin_id + 123);

    auto tk_xzzzzz_yz = pbuffer.data(idx_kin_id + 124);

    auto tk_xzzzzz_zz = pbuffer.data(idx_kin_id + 125);

#pragma omp simd aligned(pa_x,             \
                             tk_xzzzzz_xx, \
                             tk_xzzzzz_xy, \
                             tk_xzzzzz_xz, \
                             tk_xzzzzz_yy, \
                             tk_xzzzzz_yz, \
                             tk_xzzzzz_zz, \
                             tk_zzzzz_x,   \
                             tk_zzzzz_xx,  \
                             tk_zzzzz_xy,  \
                             tk_zzzzz_xz,  \
                             tk_zzzzz_y,   \
                             tk_zzzzz_yy,  \
                             tk_zzzzz_yz,  \
                             tk_zzzzz_z,   \
                             tk_zzzzz_zz,  \
                             ts_xzzzzz_xx, \
                             ts_xzzzzz_xy, \
                             ts_xzzzzz_xz, \
                             ts_xzzzzz_yy, \
                             ts_xzzzzz_yz, \
                             ts_xzzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_xzzzzz_xx[i] = 2.0 * tk_zzzzz_x[i] * fe_0 + tk_zzzzz_xx[i] * pa_x[i] + 2.0 * ts_xzzzzz_xx[i] * fz_0;

        tk_xzzzzz_xy[i] = tk_zzzzz_y[i] * fe_0 + tk_zzzzz_xy[i] * pa_x[i] + 2.0 * ts_xzzzzz_xy[i] * fz_0;

        tk_xzzzzz_xz[i] = tk_zzzzz_z[i] * fe_0 + tk_zzzzz_xz[i] * pa_x[i] + 2.0 * ts_xzzzzz_xz[i] * fz_0;

        tk_xzzzzz_yy[i] = tk_zzzzz_yy[i] * pa_x[i] + 2.0 * ts_xzzzzz_yy[i] * fz_0;

        tk_xzzzzz_yz[i] = tk_zzzzz_yz[i] * pa_x[i] + 2.0 * ts_xzzzzz_yz[i] * fz_0;

        tk_xzzzzz_zz[i] = tk_zzzzz_zz[i] * pa_x[i] + 2.0 * ts_xzzzzz_zz[i] * fz_0;
    }

    // Set up 126-132 components of targeted buffer : ID

    auto tk_yyyyyy_xx = pbuffer.data(idx_kin_id + 126);

    auto tk_yyyyyy_xy = pbuffer.data(idx_kin_id + 127);

    auto tk_yyyyyy_xz = pbuffer.data(idx_kin_id + 128);

    auto tk_yyyyyy_yy = pbuffer.data(idx_kin_id + 129);

    auto tk_yyyyyy_yz = pbuffer.data(idx_kin_id + 130);

    auto tk_yyyyyy_zz = pbuffer.data(idx_kin_id + 131);

#pragma omp simd aligned(pa_y,             \
                             tk_yyyy_xx,   \
                             tk_yyyy_xy,   \
                             tk_yyyy_xz,   \
                             tk_yyyy_yy,   \
                             tk_yyyy_yz,   \
                             tk_yyyy_zz,   \
                             tk_yyyyy_x,   \
                             tk_yyyyy_xx,  \
                             tk_yyyyy_xy,  \
                             tk_yyyyy_xz,  \
                             tk_yyyyy_y,   \
                             tk_yyyyy_yy,  \
                             tk_yyyyy_yz,  \
                             tk_yyyyy_z,   \
                             tk_yyyyy_zz,  \
                             tk_yyyyyy_xx, \
                             tk_yyyyyy_xy, \
                             tk_yyyyyy_xz, \
                             tk_yyyyyy_yy, \
                             tk_yyyyyy_yz, \
                             tk_yyyyyy_zz, \
                             ts_yyyy_xx,   \
                             ts_yyyy_xy,   \
                             ts_yyyy_xz,   \
                             ts_yyyy_yy,   \
                             ts_yyyy_yz,   \
                             ts_yyyy_zz,   \
                             ts_yyyyyy_xx, \
                             ts_yyyyyy_xy, \
                             ts_yyyyyy_xz, \
                             ts_yyyyyy_yy, \
                             ts_yyyyyy_yz, \
                             ts_yyyyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyyy_xx[i] = -10.0 * ts_yyyy_xx[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xx[i] * fe_0 + tk_yyyyy_xx[i] * pa_y[i] + 2.0 * ts_yyyyyy_xx[i] * fz_0;

        tk_yyyyyy_xy[i] = -10.0 * ts_yyyy_xy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xy[i] * fe_0 + tk_yyyyy_x[i] * fe_0 + tk_yyyyy_xy[i] * pa_y[i] +
                          2.0 * ts_yyyyyy_xy[i] * fz_0;

        tk_yyyyyy_xz[i] = -10.0 * ts_yyyy_xz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_xz[i] * fe_0 + tk_yyyyy_xz[i] * pa_y[i] + 2.0 * ts_yyyyyy_xz[i] * fz_0;

        tk_yyyyyy_yy[i] = -10.0 * ts_yyyy_yy[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yy[i] * fe_0 + 2.0 * tk_yyyyy_y[i] * fe_0 + tk_yyyyy_yy[i] * pa_y[i] +
                          2.0 * ts_yyyyyy_yy[i] * fz_0;

        tk_yyyyyy_yz[i] = -10.0 * ts_yyyy_yz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_yz[i] * fe_0 + tk_yyyyy_z[i] * fe_0 + tk_yyyyy_yz[i] * pa_y[i] +
                          2.0 * ts_yyyyyy_yz[i] * fz_0;

        tk_yyyyyy_zz[i] = -10.0 * ts_yyyy_zz[i] * fbe_0 * fz_0 + 5.0 * tk_yyyy_zz[i] * fe_0 + tk_yyyyy_zz[i] * pa_y[i] + 2.0 * ts_yyyyyy_zz[i] * fz_0;
    }

    // Set up 132-138 components of targeted buffer : ID

    auto tk_yyyyyz_xx = pbuffer.data(idx_kin_id + 132);

    auto tk_yyyyyz_xy = pbuffer.data(idx_kin_id + 133);

    auto tk_yyyyyz_xz = pbuffer.data(idx_kin_id + 134);

    auto tk_yyyyyz_yy = pbuffer.data(idx_kin_id + 135);

    auto tk_yyyyyz_yz = pbuffer.data(idx_kin_id + 136);

    auto tk_yyyyyz_zz = pbuffer.data(idx_kin_id + 137);

#pragma omp simd aligned(pa_z,             \
                             tk_yyyyy_x,   \
                             tk_yyyyy_xx,  \
                             tk_yyyyy_xy,  \
                             tk_yyyyy_xz,  \
                             tk_yyyyy_y,   \
                             tk_yyyyy_yy,  \
                             tk_yyyyy_yz,  \
                             tk_yyyyy_z,   \
                             tk_yyyyy_zz,  \
                             tk_yyyyyz_xx, \
                             tk_yyyyyz_xy, \
                             tk_yyyyyz_xz, \
                             tk_yyyyyz_yy, \
                             tk_yyyyyz_yz, \
                             tk_yyyyyz_zz, \
                             ts_yyyyyz_xx, \
                             ts_yyyyyz_xy, \
                             ts_yyyyyz_xz, \
                             ts_yyyyyz_yy, \
                             ts_yyyyyz_yz, \
                             ts_yyyyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yyyyyz_xx[i] = tk_yyyyy_xx[i] * pa_z[i] + 2.0 * ts_yyyyyz_xx[i] * fz_0;

        tk_yyyyyz_xy[i] = tk_yyyyy_xy[i] * pa_z[i] + 2.0 * ts_yyyyyz_xy[i] * fz_0;

        tk_yyyyyz_xz[i] = tk_yyyyy_x[i] * fe_0 + tk_yyyyy_xz[i] * pa_z[i] + 2.0 * ts_yyyyyz_xz[i] * fz_0;

        tk_yyyyyz_yy[i] = tk_yyyyy_yy[i] * pa_z[i] + 2.0 * ts_yyyyyz_yy[i] * fz_0;

        tk_yyyyyz_yz[i] = tk_yyyyy_y[i] * fe_0 + tk_yyyyy_yz[i] * pa_z[i] + 2.0 * ts_yyyyyz_yz[i] * fz_0;

        tk_yyyyyz_zz[i] = 2.0 * tk_yyyyy_z[i] * fe_0 + tk_yyyyy_zz[i] * pa_z[i] + 2.0 * ts_yyyyyz_zz[i] * fz_0;
    }

    // Set up 138-144 components of targeted buffer : ID

    auto tk_yyyyzz_xx = pbuffer.data(idx_kin_id + 138);

    auto tk_yyyyzz_xy = pbuffer.data(idx_kin_id + 139);

    auto tk_yyyyzz_xz = pbuffer.data(idx_kin_id + 140);

    auto tk_yyyyzz_yy = pbuffer.data(idx_kin_id + 141);

    auto tk_yyyyzz_yz = pbuffer.data(idx_kin_id + 142);

    auto tk_yyyyzz_zz = pbuffer.data(idx_kin_id + 143);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_yyyy_xy,   \
                             tk_yyyy_yy,   \
                             tk_yyyyz_xy,  \
                             tk_yyyyz_yy,  \
                             tk_yyyyzz_xx, \
                             tk_yyyyzz_xy, \
                             tk_yyyyzz_xz, \
                             tk_yyyyzz_yy, \
                             tk_yyyyzz_yz, \
                             tk_yyyyzz_zz, \
                             tk_yyyzz_xx,  \
                             tk_yyyzz_xz,  \
                             tk_yyyzz_yz,  \
                             tk_yyyzz_z,   \
                             tk_yyyzz_zz,  \
                             tk_yyzz_xx,   \
                             tk_yyzz_xz,   \
                             tk_yyzz_yz,   \
                             tk_yyzz_zz,   \
                             ts_yyyy_xy,   \
                             ts_yyyy_yy,   \
                             ts_yyyyzz_xx, \
                             ts_yyyyzz_xy, \
                             ts_yyyyzz_xz, \
                             ts_yyyyzz_yy, \
                             ts_yyyyzz_yz, \
                             ts_yyyyzz_zz, \
                             ts_yyzz_xx,   \
                             ts_yyzz_xz,   \
                             ts_yyzz_yz,   \
                             ts_yyzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyyzz_xx[i] = -6.0 * ts_yyzz_xx[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xx[i] * fe_0 + tk_yyyzz_xx[i] * pa_y[i] + 2.0 * ts_yyyyzz_xx[i] * fz_0;

        tk_yyyyzz_xy[i] = -2.0 * ts_yyyy_xy[i] * fbe_0 * fz_0 + tk_yyyy_xy[i] * fe_0 + tk_yyyyz_xy[i] * pa_z[i] + 2.0 * ts_yyyyzz_xy[i] * fz_0;

        tk_yyyyzz_xz[i] = -6.0 * ts_yyzz_xz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xz[i] * fe_0 + tk_yyyzz_xz[i] * pa_y[i] + 2.0 * ts_yyyyzz_xz[i] * fz_0;

        tk_yyyyzz_yy[i] = -2.0 * ts_yyyy_yy[i] * fbe_0 * fz_0 + tk_yyyy_yy[i] * fe_0 + tk_yyyyz_yy[i] * pa_z[i] + 2.0 * ts_yyyyzz_yy[i] * fz_0;

        tk_yyyyzz_yz[i] = -6.0 * ts_yyzz_yz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yz[i] * fe_0 + tk_yyyzz_z[i] * fe_0 + tk_yyyzz_yz[i] * pa_y[i] +
                          2.0 * ts_yyyyzz_yz[i] * fz_0;

        tk_yyyyzz_zz[i] = -6.0 * ts_yyzz_zz[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_zz[i] * fe_0 + tk_yyyzz_zz[i] * pa_y[i] + 2.0 * ts_yyyyzz_zz[i] * fz_0;
    }

    // Set up 144-150 components of targeted buffer : ID

    auto tk_yyyzzz_xx = pbuffer.data(idx_kin_id + 144);

    auto tk_yyyzzz_xy = pbuffer.data(idx_kin_id + 145);

    auto tk_yyyzzz_xz = pbuffer.data(idx_kin_id + 146);

    auto tk_yyyzzz_yy = pbuffer.data(idx_kin_id + 147);

    auto tk_yyyzzz_yz = pbuffer.data(idx_kin_id + 148);

    auto tk_yyyzzz_zz = pbuffer.data(idx_kin_id + 149);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_yyyz_xy,   \
                             tk_yyyz_yy,   \
                             tk_yyyzz_xy,  \
                             tk_yyyzz_yy,  \
                             tk_yyyzzz_xx, \
                             tk_yyyzzz_xy, \
                             tk_yyyzzz_xz, \
                             tk_yyyzzz_yy, \
                             tk_yyyzzz_yz, \
                             tk_yyyzzz_zz, \
                             tk_yyzzz_xx,  \
                             tk_yyzzz_xz,  \
                             tk_yyzzz_yz,  \
                             tk_yyzzz_z,   \
                             tk_yyzzz_zz,  \
                             tk_yzzz_xx,   \
                             tk_yzzz_xz,   \
                             tk_yzzz_yz,   \
                             tk_yzzz_zz,   \
                             ts_yyyz_xy,   \
                             ts_yyyz_yy,   \
                             ts_yyyzzz_xx, \
                             ts_yyyzzz_xy, \
                             ts_yyyzzz_xz, \
                             ts_yyyzzz_yy, \
                             ts_yyyzzz_yz, \
                             ts_yyyzzz_zz, \
                             ts_yzzz_xx,   \
                             ts_yzzz_xz,   \
                             ts_yzzz_yz,   \
                             ts_yzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyyzzz_xx[i] = -4.0 * ts_yzzz_xx[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xx[i] * fe_0 + tk_yyzzz_xx[i] * pa_y[i] + 2.0 * ts_yyyzzz_xx[i] * fz_0;

        tk_yyyzzz_xy[i] = -4.0 * ts_yyyz_xy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_xy[i] * fe_0 + tk_yyyzz_xy[i] * pa_z[i] + 2.0 * ts_yyyzzz_xy[i] * fz_0;

        tk_yyyzzz_xz[i] = -4.0 * ts_yzzz_xz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_xz[i] * fe_0 + tk_yyzzz_xz[i] * pa_y[i] + 2.0 * ts_yyyzzz_xz[i] * fz_0;

        tk_yyyzzz_yy[i] = -4.0 * ts_yyyz_yy[i] * fbe_0 * fz_0 + 2.0 * tk_yyyz_yy[i] * fe_0 + tk_yyyzz_yy[i] * pa_z[i] + 2.0 * ts_yyyzzz_yy[i] * fz_0;

        tk_yyyzzz_yz[i] = -4.0 * ts_yzzz_yz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_yz[i] * fe_0 + tk_yyzzz_z[i] * fe_0 + tk_yyzzz_yz[i] * pa_y[i] +
                          2.0 * ts_yyyzzz_yz[i] * fz_0;

        tk_yyyzzz_zz[i] = -4.0 * ts_yzzz_zz[i] * fbe_0 * fz_0 + 2.0 * tk_yzzz_zz[i] * fe_0 + tk_yyzzz_zz[i] * pa_y[i] + 2.0 * ts_yyyzzz_zz[i] * fz_0;
    }

    // Set up 150-156 components of targeted buffer : ID

    auto tk_yyzzzz_xx = pbuffer.data(idx_kin_id + 150);

    auto tk_yyzzzz_xy = pbuffer.data(idx_kin_id + 151);

    auto tk_yyzzzz_xz = pbuffer.data(idx_kin_id + 152);

    auto tk_yyzzzz_yy = pbuffer.data(idx_kin_id + 153);

    auto tk_yyzzzz_yz = pbuffer.data(idx_kin_id + 154);

    auto tk_yyzzzz_zz = pbuffer.data(idx_kin_id + 155);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tk_yyzz_xy,   \
                             tk_yyzz_yy,   \
                             tk_yyzzz_xy,  \
                             tk_yyzzz_yy,  \
                             tk_yyzzzz_xx, \
                             tk_yyzzzz_xy, \
                             tk_yyzzzz_xz, \
                             tk_yyzzzz_yy, \
                             tk_yyzzzz_yz, \
                             tk_yyzzzz_zz, \
                             tk_yzzzz_xx,  \
                             tk_yzzzz_xz,  \
                             tk_yzzzz_yz,  \
                             tk_yzzzz_z,   \
                             tk_yzzzz_zz,  \
                             tk_zzzz_xx,   \
                             tk_zzzz_xz,   \
                             tk_zzzz_yz,   \
                             tk_zzzz_zz,   \
                             ts_yyzz_xy,   \
                             ts_yyzz_yy,   \
                             ts_yyzzzz_xx, \
                             ts_yyzzzz_xy, \
                             ts_yyzzzz_xz, \
                             ts_yyzzzz_yy, \
                             ts_yyzzzz_yz, \
                             ts_yyzzzz_zz, \
                             ts_zzzz_xx,   \
                             ts_zzzz_xz,   \
                             ts_zzzz_yz,   \
                             ts_zzzz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_yyzzzz_xx[i] = -2.0 * ts_zzzz_xx[i] * fbe_0 * fz_0 + tk_zzzz_xx[i] * fe_0 + tk_yzzzz_xx[i] * pa_y[i] + 2.0 * ts_yyzzzz_xx[i] * fz_0;

        tk_yyzzzz_xy[i] = -6.0 * ts_yyzz_xy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_xy[i] * fe_0 + tk_yyzzz_xy[i] * pa_z[i] + 2.0 * ts_yyzzzz_xy[i] * fz_0;

        tk_yyzzzz_xz[i] = -2.0 * ts_zzzz_xz[i] * fbe_0 * fz_0 + tk_zzzz_xz[i] * fe_0 + tk_yzzzz_xz[i] * pa_y[i] + 2.0 * ts_yyzzzz_xz[i] * fz_0;

        tk_yyzzzz_yy[i] = -6.0 * ts_yyzz_yy[i] * fbe_0 * fz_0 + 3.0 * tk_yyzz_yy[i] * fe_0 + tk_yyzzz_yy[i] * pa_z[i] + 2.0 * ts_yyzzzz_yy[i] * fz_0;

        tk_yyzzzz_yz[i] = -2.0 * ts_zzzz_yz[i] * fbe_0 * fz_0 + tk_zzzz_yz[i] * fe_0 + tk_yzzzz_z[i] * fe_0 + tk_yzzzz_yz[i] * pa_y[i] +
                          2.0 * ts_yyzzzz_yz[i] * fz_0;

        tk_yyzzzz_zz[i] = -2.0 * ts_zzzz_zz[i] * fbe_0 * fz_0 + tk_zzzz_zz[i] * fe_0 + tk_yzzzz_zz[i] * pa_y[i] + 2.0 * ts_yyzzzz_zz[i] * fz_0;
    }

    // Set up 156-162 components of targeted buffer : ID

    auto tk_yzzzzz_xx = pbuffer.data(idx_kin_id + 156);

    auto tk_yzzzzz_xy = pbuffer.data(idx_kin_id + 157);

    auto tk_yzzzzz_xz = pbuffer.data(idx_kin_id + 158);

    auto tk_yzzzzz_yy = pbuffer.data(idx_kin_id + 159);

    auto tk_yzzzzz_yz = pbuffer.data(idx_kin_id + 160);

    auto tk_yzzzzz_zz = pbuffer.data(idx_kin_id + 161);

#pragma omp simd aligned(pa_y,             \
                             tk_yzzzzz_xx, \
                             tk_yzzzzz_xy, \
                             tk_yzzzzz_xz, \
                             tk_yzzzzz_yy, \
                             tk_yzzzzz_yz, \
                             tk_yzzzzz_zz, \
                             tk_zzzzz_x,   \
                             tk_zzzzz_xx,  \
                             tk_zzzzz_xy,  \
                             tk_zzzzz_xz,  \
                             tk_zzzzz_y,   \
                             tk_zzzzz_yy,  \
                             tk_zzzzz_yz,  \
                             tk_zzzzz_z,   \
                             tk_zzzzz_zz,  \
                             ts_yzzzzz_xx, \
                             ts_yzzzzz_xy, \
                             ts_yzzzzz_xz, \
                             ts_yzzzzz_yy, \
                             ts_yzzzzz_yz, \
                             ts_yzzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        tk_yzzzzz_xx[i] = tk_zzzzz_xx[i] * pa_y[i] + 2.0 * ts_yzzzzz_xx[i] * fz_0;

        tk_yzzzzz_xy[i] = tk_zzzzz_x[i] * fe_0 + tk_zzzzz_xy[i] * pa_y[i] + 2.0 * ts_yzzzzz_xy[i] * fz_0;

        tk_yzzzzz_xz[i] = tk_zzzzz_xz[i] * pa_y[i] + 2.0 * ts_yzzzzz_xz[i] * fz_0;

        tk_yzzzzz_yy[i] = 2.0 * tk_zzzzz_y[i] * fe_0 + tk_zzzzz_yy[i] * pa_y[i] + 2.0 * ts_yzzzzz_yy[i] * fz_0;

        tk_yzzzzz_yz[i] = tk_zzzzz_z[i] * fe_0 + tk_zzzzz_yz[i] * pa_y[i] + 2.0 * ts_yzzzzz_yz[i] * fz_0;

        tk_yzzzzz_zz[i] = tk_zzzzz_zz[i] * pa_y[i] + 2.0 * ts_yzzzzz_zz[i] * fz_0;
    }

    // Set up 162-168 components of targeted buffer : ID

    auto tk_zzzzzz_xx = pbuffer.data(idx_kin_id + 162);

    auto tk_zzzzzz_xy = pbuffer.data(idx_kin_id + 163);

    auto tk_zzzzzz_xz = pbuffer.data(idx_kin_id + 164);

    auto tk_zzzzzz_yy = pbuffer.data(idx_kin_id + 165);

    auto tk_zzzzzz_yz = pbuffer.data(idx_kin_id + 166);

    auto tk_zzzzzz_zz = pbuffer.data(idx_kin_id + 167);

#pragma omp simd aligned(pa_z,             \
                             tk_zzzz_xx,   \
                             tk_zzzz_xy,   \
                             tk_zzzz_xz,   \
                             tk_zzzz_yy,   \
                             tk_zzzz_yz,   \
                             tk_zzzz_zz,   \
                             tk_zzzzz_x,   \
                             tk_zzzzz_xx,  \
                             tk_zzzzz_xy,  \
                             tk_zzzzz_xz,  \
                             tk_zzzzz_y,   \
                             tk_zzzzz_yy,  \
                             tk_zzzzz_yz,  \
                             tk_zzzzz_z,   \
                             tk_zzzzz_zz,  \
                             tk_zzzzzz_xx, \
                             tk_zzzzzz_xy, \
                             tk_zzzzzz_xz, \
                             tk_zzzzzz_yy, \
                             tk_zzzzzz_yz, \
                             tk_zzzzzz_zz, \
                             ts_zzzz_xx,   \
                             ts_zzzz_xy,   \
                             ts_zzzz_xz,   \
                             ts_zzzz_yy,   \
                             ts_zzzz_yz,   \
                             ts_zzzz_zz,   \
                             ts_zzzzzz_xx, \
                             ts_zzzzzz_xy, \
                             ts_zzzzzz_xz, \
                             ts_zzzzzz_yy, \
                             ts_zzzzzz_yz, \
                             ts_zzzzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;

        const double fbe_0 = 0.5 / a_exp;

        tk_zzzzzz_xx[i] = -10.0 * ts_zzzz_xx[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xx[i] * fe_0 + tk_zzzzz_xx[i] * pa_z[i] + 2.0 * ts_zzzzzz_xx[i] * fz_0;

        tk_zzzzzz_xy[i] = -10.0 * ts_zzzz_xy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xy[i] * fe_0 + tk_zzzzz_xy[i] * pa_z[i] + 2.0 * ts_zzzzzz_xy[i] * fz_0;

        tk_zzzzzz_xz[i] = -10.0 * ts_zzzz_xz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_xz[i] * fe_0 + tk_zzzzz_x[i] * fe_0 + tk_zzzzz_xz[i] * pa_z[i] +
                          2.0 * ts_zzzzzz_xz[i] * fz_0;

        tk_zzzzzz_yy[i] = -10.0 * ts_zzzz_yy[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yy[i] * fe_0 + tk_zzzzz_yy[i] * pa_z[i] + 2.0 * ts_zzzzzz_yy[i] * fz_0;

        tk_zzzzzz_yz[i] = -10.0 * ts_zzzz_yz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_yz[i] * fe_0 + tk_zzzzz_y[i] * fe_0 + tk_zzzzz_yz[i] * pa_z[i] +
                          2.0 * ts_zzzzzz_yz[i] * fz_0;

        tk_zzzzzz_zz[i] = -10.0 * ts_zzzz_zz[i] * fbe_0 * fz_0 + 5.0 * tk_zzzz_zz[i] * fe_0 + 2.0 * tk_zzzzz_z[i] * fe_0 + tk_zzzzz_zz[i] * pa_z[i] +
                          2.0 * ts_zzzzzz_zz[i] * fz_0;
    }
}

}  // namespace kinrec
