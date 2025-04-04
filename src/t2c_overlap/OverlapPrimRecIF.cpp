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

#include "OverlapPrimRecIF.hpp"

namespace ovlrec { // ovlrec namespace

auto
comp_prim_overlap_if(CSimdArray<double>& pbuffer, 
                     const size_t idx_ovl_if,
                     const size_t idx_ovl_gf,
                     const size_t idx_ovl_hd,
                     const size_t idx_ovl_hf,
                     const CSimdArray<double>& factors,
                     const size_t idx_rpa,
                     const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

    auto ts_xxxy_xxx = pbuffer.data(idx_ovl_gf + 10);

    auto ts_xxxy_xxz = pbuffer.data(idx_ovl_gf + 12);

    auto ts_xxxy_xzz = pbuffer.data(idx_ovl_gf + 15);

    auto ts_xxxy_yyy = pbuffer.data(idx_ovl_gf + 16);

    auto ts_xxxy_yyz = pbuffer.data(idx_ovl_gf + 17);

    auto ts_xxxy_yzz = pbuffer.data(idx_ovl_gf + 18);

    auto ts_xxxz_xxx = pbuffer.data(idx_ovl_gf + 20);

    auto ts_xxxz_xxy = pbuffer.data(idx_ovl_gf + 21);

    auto ts_xxxz_xxz = pbuffer.data(idx_ovl_gf + 22);

    auto ts_xxxz_xyy = pbuffer.data(idx_ovl_gf + 23);

    auto ts_xxxz_xzz = pbuffer.data(idx_ovl_gf + 25);

    auto ts_xxxz_yyz = pbuffer.data(idx_ovl_gf + 27);

    auto ts_xxxz_yzz = pbuffer.data(idx_ovl_gf + 28);

    auto ts_xxxz_zzz = pbuffer.data(idx_ovl_gf + 29);

    auto ts_xxyy_xxx = pbuffer.data(idx_ovl_gf + 30);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xxz = pbuffer.data(idx_ovl_gf + 32);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_xyz = pbuffer.data(idx_ovl_gf + 34);

    auto ts_xxyy_xzz = pbuffer.data(idx_ovl_gf + 35);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxyy_zzz = pbuffer.data(idx_ovl_gf + 39);

    auto ts_xxyz_xxz = pbuffer.data(idx_ovl_gf + 42);

    auto ts_xxyz_xzz = pbuffer.data(idx_ovl_gf + 45);

    auto ts_xxyz_yyz = pbuffer.data(idx_ovl_gf + 47);

    auto ts_xxyz_yzz = pbuffer.data(idx_ovl_gf + 48);

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxy = pbuffer.data(idx_ovl_gf + 51);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xyy = pbuffer.data(idx_ovl_gf + 53);

    auto ts_xxzz_xyz = pbuffer.data(idx_ovl_gf + 54);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyy = pbuffer.data(idx_ovl_gf + 56);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

    auto ts_xyyy_xxy = pbuffer.data(idx_ovl_gf + 61);

    auto ts_xyyy_xyy = pbuffer.data(idx_ovl_gf + 63);

    auto ts_xyyy_xyz = pbuffer.data(idx_ovl_gf + 64);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xyyy_zzz = pbuffer.data(idx_ovl_gf + 69);

    auto ts_xyyz_yyz = pbuffer.data(idx_ovl_gf + 77);

    auto ts_xyyz_yzz = pbuffer.data(idx_ovl_gf + 78);

    auto ts_xyyz_zzz = pbuffer.data(idx_ovl_gf + 79);

    auto ts_xyzz_yyy = pbuffer.data(idx_ovl_gf + 86);

    auto ts_xyzz_yyz = pbuffer.data(idx_ovl_gf + 87);

    auto ts_xyzz_yzz = pbuffer.data(idx_ovl_gf + 88);

    auto ts_xzzz_xxz = pbuffer.data(idx_ovl_gf + 92);

    auto ts_xzzz_xyz = pbuffer.data(idx_ovl_gf + 94);

    auto ts_xzzz_xzz = pbuffer.data(idx_ovl_gf + 95);

    auto ts_xzzz_yyy = pbuffer.data(idx_ovl_gf + 96);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

    auto ts_yyyz_xxy = pbuffer.data(idx_ovl_gf + 111);

    auto ts_yyyz_xxz = pbuffer.data(idx_ovl_gf + 112);

    auto ts_yyyz_xyy = pbuffer.data(idx_ovl_gf + 113);

    auto ts_yyyz_xzz = pbuffer.data(idx_ovl_gf + 115);

    auto ts_yyyz_yyy = pbuffer.data(idx_ovl_gf + 116);

    auto ts_yyyz_yyz = pbuffer.data(idx_ovl_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_ovl_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_ovl_gf + 119);

    auto ts_yyzz_xxx = pbuffer.data(idx_ovl_gf + 120);

    auto ts_yyzz_xxy = pbuffer.data(idx_ovl_gf + 121);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyy = pbuffer.data(idx_ovl_gf + 123);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

    auto ts_yzzz_xxx = pbuffer.data(idx_ovl_gf + 130);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xyz = pbuffer.data(idx_ovl_gf + 134);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_ovl_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

    // Set up components of auxiliary buffer : HD

    auto ts_xxxxx_xx = pbuffer.data(idx_ovl_hd);

    auto ts_xxxxx_xy = pbuffer.data(idx_ovl_hd + 1);

    auto ts_xxxxx_xz = pbuffer.data(idx_ovl_hd + 2);

    auto ts_xxxxx_yy = pbuffer.data(idx_ovl_hd + 3);

    auto ts_xxxxx_yz = pbuffer.data(idx_ovl_hd + 4);

    auto ts_xxxxx_zz = pbuffer.data(idx_ovl_hd + 5);

    auto ts_xxxxz_xz = pbuffer.data(idx_ovl_hd + 14);

    auto ts_xxxyy_xy = pbuffer.data(idx_ovl_hd + 19);

    auto ts_xxxyy_yy = pbuffer.data(idx_ovl_hd + 21);

    auto ts_xxxyy_yz = pbuffer.data(idx_ovl_hd + 22);

    auto ts_xxxzz_xx = pbuffer.data(idx_ovl_hd + 30);

    auto ts_xxxzz_xy = pbuffer.data(idx_ovl_hd + 31);

    auto ts_xxxzz_xz = pbuffer.data(idx_ovl_hd + 32);

    auto ts_xxxzz_yz = pbuffer.data(idx_ovl_hd + 34);

    auto ts_xxxzz_zz = pbuffer.data(idx_ovl_hd + 35);

    auto ts_xxyyy_xy = pbuffer.data(idx_ovl_hd + 37);

    auto ts_xxyyy_yy = pbuffer.data(idx_ovl_hd + 39);

    auto ts_xxyyy_yz = pbuffer.data(idx_ovl_hd + 40);

    auto ts_xxzzz_xx = pbuffer.data(idx_ovl_hd + 54);

    auto ts_xxzzz_xy = pbuffer.data(idx_ovl_hd + 55);

    auto ts_xxzzz_xz = pbuffer.data(idx_ovl_hd + 56);

    auto ts_xxzzz_yz = pbuffer.data(idx_ovl_hd + 58);

    auto ts_xxzzz_zz = pbuffer.data(idx_ovl_hd + 59);

    auto ts_xyyyy_xy = pbuffer.data(idx_ovl_hd + 61);

    auto ts_xyyyy_yy = pbuffer.data(idx_ovl_hd + 63);

    auto ts_xyyyy_yz = pbuffer.data(idx_ovl_hd + 64);

    auto ts_xyyzz_yz = pbuffer.data(idx_ovl_hd + 76);

    auto ts_xzzzz_xz = pbuffer.data(idx_ovl_hd + 86);

    auto ts_xzzzz_yz = pbuffer.data(idx_ovl_hd + 88);

    auto ts_xzzzz_zz = pbuffer.data(idx_ovl_hd + 89);

    auto ts_yyyyy_xx = pbuffer.data(idx_ovl_hd + 90);

    auto ts_yyyyy_xy = pbuffer.data(idx_ovl_hd + 91);

    auto ts_yyyyy_xz = pbuffer.data(idx_ovl_hd + 92);

    auto ts_yyyyy_yy = pbuffer.data(idx_ovl_hd + 93);

    auto ts_yyyyy_yz = pbuffer.data(idx_ovl_hd + 94);

    auto ts_yyyyy_zz = pbuffer.data(idx_ovl_hd + 95);

    auto ts_yyyyz_xz = pbuffer.data(idx_ovl_hd + 98);

    auto ts_yyyyz_yz = pbuffer.data(idx_ovl_hd + 100);

    auto ts_yyyyz_zz = pbuffer.data(idx_ovl_hd + 101);

    auto ts_yyyzz_xx = pbuffer.data(idx_ovl_hd + 102);

    auto ts_yyyzz_xy = pbuffer.data(idx_ovl_hd + 103);

    auto ts_yyyzz_xz = pbuffer.data(idx_ovl_hd + 104);

    auto ts_yyyzz_yy = pbuffer.data(idx_ovl_hd + 105);

    auto ts_yyyzz_yz = pbuffer.data(idx_ovl_hd + 106);

    auto ts_yyyzz_zz = pbuffer.data(idx_ovl_hd + 107);

    auto ts_yyzzz_xx = pbuffer.data(idx_ovl_hd + 108);

    auto ts_yyzzz_xy = pbuffer.data(idx_ovl_hd + 109);

    auto ts_yyzzz_xz = pbuffer.data(idx_ovl_hd + 110);

    auto ts_yyzzz_yy = pbuffer.data(idx_ovl_hd + 111);

    auto ts_yyzzz_yz = pbuffer.data(idx_ovl_hd + 112);

    auto ts_yyzzz_zz = pbuffer.data(idx_ovl_hd + 113);

    auto ts_yzzzz_xy = pbuffer.data(idx_ovl_hd + 115);

    auto ts_yzzzz_xz = pbuffer.data(idx_ovl_hd + 116);

    auto ts_yzzzz_yy = pbuffer.data(idx_ovl_hd + 117);

    auto ts_yzzzz_yz = pbuffer.data(idx_ovl_hd + 118);

    auto ts_yzzzz_zz = pbuffer.data(idx_ovl_hd + 119);

    auto ts_zzzzz_xx = pbuffer.data(idx_ovl_hd + 120);

    auto ts_zzzzz_xy = pbuffer.data(idx_ovl_hd + 121);

    auto ts_zzzzz_xz = pbuffer.data(idx_ovl_hd + 122);

    auto ts_zzzzz_yy = pbuffer.data(idx_ovl_hd + 123);

    auto ts_zzzzz_yz = pbuffer.data(idx_ovl_hd + 124);

    auto ts_zzzzz_zz = pbuffer.data(idx_ovl_hd + 125);

    // Set up components of auxiliary buffer : HF

    auto ts_xxxxx_xxx = pbuffer.data(idx_ovl_hf);

    auto ts_xxxxx_xxy = pbuffer.data(idx_ovl_hf + 1);

    auto ts_xxxxx_xxz = pbuffer.data(idx_ovl_hf + 2);

    auto ts_xxxxx_xyy = pbuffer.data(idx_ovl_hf + 3);

    auto ts_xxxxx_xyz = pbuffer.data(idx_ovl_hf + 4);

    auto ts_xxxxx_xzz = pbuffer.data(idx_ovl_hf + 5);

    auto ts_xxxxx_yyy = pbuffer.data(idx_ovl_hf + 6);

    auto ts_xxxxx_yyz = pbuffer.data(idx_ovl_hf + 7);

    auto ts_xxxxx_yzz = pbuffer.data(idx_ovl_hf + 8);

    auto ts_xxxxx_zzz = pbuffer.data(idx_ovl_hf + 9);

    auto ts_xxxxy_xxx = pbuffer.data(idx_ovl_hf + 10);

    auto ts_xxxxy_xxy = pbuffer.data(idx_ovl_hf + 11);

    auto ts_xxxxy_xxz = pbuffer.data(idx_ovl_hf + 12);

    auto ts_xxxxy_xyy = pbuffer.data(idx_ovl_hf + 13);

    auto ts_xxxxy_xzz = pbuffer.data(idx_ovl_hf + 15);

    auto ts_xxxxy_yyy = pbuffer.data(idx_ovl_hf + 16);

    auto ts_xxxxy_yyz = pbuffer.data(idx_ovl_hf + 17);

    auto ts_xxxxy_yzz = pbuffer.data(idx_ovl_hf + 18);

    auto ts_xxxxz_xxx = pbuffer.data(idx_ovl_hf + 20);

    auto ts_xxxxz_xxy = pbuffer.data(idx_ovl_hf + 21);

    auto ts_xxxxz_xxz = pbuffer.data(idx_ovl_hf + 22);

    auto ts_xxxxz_xyy = pbuffer.data(idx_ovl_hf + 23);

    auto ts_xxxxz_xyz = pbuffer.data(idx_ovl_hf + 24);

    auto ts_xxxxz_xzz = pbuffer.data(idx_ovl_hf + 25);

    auto ts_xxxxz_yyz = pbuffer.data(idx_ovl_hf + 27);

    auto ts_xxxxz_yzz = pbuffer.data(idx_ovl_hf + 28);

    auto ts_xxxxz_zzz = pbuffer.data(idx_ovl_hf + 29);

    auto ts_xxxyy_xxx = pbuffer.data(idx_ovl_hf + 30);

    auto ts_xxxyy_xxy = pbuffer.data(idx_ovl_hf + 31);

    auto ts_xxxyy_xxz = pbuffer.data(idx_ovl_hf + 32);

    auto ts_xxxyy_xyy = pbuffer.data(idx_ovl_hf + 33);

    auto ts_xxxyy_xyz = pbuffer.data(idx_ovl_hf + 34);

    auto ts_xxxyy_xzz = pbuffer.data(idx_ovl_hf + 35);

    auto ts_xxxyy_yyy = pbuffer.data(idx_ovl_hf + 36);

    auto ts_xxxyy_yyz = pbuffer.data(idx_ovl_hf + 37);

    auto ts_xxxyy_yzz = pbuffer.data(idx_ovl_hf + 38);

    auto ts_xxxyy_zzz = pbuffer.data(idx_ovl_hf + 39);

    auto ts_xxxyz_xxz = pbuffer.data(idx_ovl_hf + 42);

    auto ts_xxxyz_xzz = pbuffer.data(idx_ovl_hf + 45);

    auto ts_xxxyz_yyz = pbuffer.data(idx_ovl_hf + 47);

    auto ts_xxxyz_yzz = pbuffer.data(idx_ovl_hf + 48);

    auto ts_xxxzz_xxx = pbuffer.data(idx_ovl_hf + 50);

    auto ts_xxxzz_xxy = pbuffer.data(idx_ovl_hf + 51);

    auto ts_xxxzz_xxz = pbuffer.data(idx_ovl_hf + 52);

    auto ts_xxxzz_xyy = pbuffer.data(idx_ovl_hf + 53);

    auto ts_xxxzz_xyz = pbuffer.data(idx_ovl_hf + 54);

    auto ts_xxxzz_xzz = pbuffer.data(idx_ovl_hf + 55);

    auto ts_xxxzz_yyy = pbuffer.data(idx_ovl_hf + 56);

    auto ts_xxxzz_yyz = pbuffer.data(idx_ovl_hf + 57);

    auto ts_xxxzz_yzz = pbuffer.data(idx_ovl_hf + 58);

    auto ts_xxxzz_zzz = pbuffer.data(idx_ovl_hf + 59);

    auto ts_xxyyy_xxx = pbuffer.data(idx_ovl_hf + 60);

    auto ts_xxyyy_xxy = pbuffer.data(idx_ovl_hf + 61);

    auto ts_xxyyy_xxz = pbuffer.data(idx_ovl_hf + 62);

    auto ts_xxyyy_xyy = pbuffer.data(idx_ovl_hf + 63);

    auto ts_xxyyy_xyz = pbuffer.data(idx_ovl_hf + 64);

    auto ts_xxyyy_xzz = pbuffer.data(idx_ovl_hf + 65);

    auto ts_xxyyy_yyy = pbuffer.data(idx_ovl_hf + 66);

    auto ts_xxyyy_yyz = pbuffer.data(idx_ovl_hf + 67);

    auto ts_xxyyy_yzz = pbuffer.data(idx_ovl_hf + 68);

    auto ts_xxyyy_zzz = pbuffer.data(idx_ovl_hf + 69);

    auto ts_xxyyz_xxy = pbuffer.data(idx_ovl_hf + 71);

    auto ts_xxyyz_xxz = pbuffer.data(idx_ovl_hf + 72);

    auto ts_xxyyz_xyy = pbuffer.data(idx_ovl_hf + 73);

    auto ts_xxyyz_xzz = pbuffer.data(idx_ovl_hf + 75);

    auto ts_xxyyz_yyz = pbuffer.data(idx_ovl_hf + 77);

    auto ts_xxyyz_yzz = pbuffer.data(idx_ovl_hf + 78);

    auto ts_xxyyz_zzz = pbuffer.data(idx_ovl_hf + 79);

    auto ts_xxyzz_xxx = pbuffer.data(idx_ovl_hf + 80);

    auto ts_xxyzz_xxz = pbuffer.data(idx_ovl_hf + 82);

    auto ts_xxyzz_xzz = pbuffer.data(idx_ovl_hf + 85);

    auto ts_xxyzz_yyy = pbuffer.data(idx_ovl_hf + 86);

    auto ts_xxyzz_yyz = pbuffer.data(idx_ovl_hf + 87);

    auto ts_xxyzz_yzz = pbuffer.data(idx_ovl_hf + 88);

    auto ts_xxzzz_xxx = pbuffer.data(idx_ovl_hf + 90);

    auto ts_xxzzz_xxy = pbuffer.data(idx_ovl_hf + 91);

    auto ts_xxzzz_xxz = pbuffer.data(idx_ovl_hf + 92);

    auto ts_xxzzz_xyy = pbuffer.data(idx_ovl_hf + 93);

    auto ts_xxzzz_xyz = pbuffer.data(idx_ovl_hf + 94);

    auto ts_xxzzz_xzz = pbuffer.data(idx_ovl_hf + 95);

    auto ts_xxzzz_yyy = pbuffer.data(idx_ovl_hf + 96);

    auto ts_xxzzz_yyz = pbuffer.data(idx_ovl_hf + 97);

    auto ts_xxzzz_yzz = pbuffer.data(idx_ovl_hf + 98);

    auto ts_xxzzz_zzz = pbuffer.data(idx_ovl_hf + 99);

    auto ts_xyyyy_xxx = pbuffer.data(idx_ovl_hf + 100);

    auto ts_xyyyy_xxy = pbuffer.data(idx_ovl_hf + 101);

    auto ts_xyyyy_xyy = pbuffer.data(idx_ovl_hf + 103);

    auto ts_xyyyy_xyz = pbuffer.data(idx_ovl_hf + 104);

    auto ts_xyyyy_yyy = pbuffer.data(idx_ovl_hf + 106);

    auto ts_xyyyy_yyz = pbuffer.data(idx_ovl_hf + 107);

    auto ts_xyyyy_yzz = pbuffer.data(idx_ovl_hf + 108);

    auto ts_xyyyy_zzz = pbuffer.data(idx_ovl_hf + 109);

    auto ts_xyyyz_yyz = pbuffer.data(idx_ovl_hf + 117);

    auto ts_xyyyz_yzz = pbuffer.data(idx_ovl_hf + 118);

    auto ts_xyyyz_zzz = pbuffer.data(idx_ovl_hf + 119);

    auto ts_xyyzz_xyz = pbuffer.data(idx_ovl_hf + 124);

    auto ts_xyyzz_yyy = pbuffer.data(idx_ovl_hf + 126);

    auto ts_xyyzz_yyz = pbuffer.data(idx_ovl_hf + 127);

    auto ts_xyyzz_yzz = pbuffer.data(idx_ovl_hf + 128);

    auto ts_xyyzz_zzz = pbuffer.data(idx_ovl_hf + 129);

    auto ts_xyzzz_yyy = pbuffer.data(idx_ovl_hf + 136);

    auto ts_xyzzz_yyz = pbuffer.data(idx_ovl_hf + 137);

    auto ts_xyzzz_yzz = pbuffer.data(idx_ovl_hf + 138);

    auto ts_xzzzz_xxx = pbuffer.data(idx_ovl_hf + 140);

    auto ts_xzzzz_xxz = pbuffer.data(idx_ovl_hf + 142);

    auto ts_xzzzz_xyz = pbuffer.data(idx_ovl_hf + 144);

    auto ts_xzzzz_xzz = pbuffer.data(idx_ovl_hf + 145);

    auto ts_xzzzz_yyy = pbuffer.data(idx_ovl_hf + 146);

    auto ts_xzzzz_yyz = pbuffer.data(idx_ovl_hf + 147);

    auto ts_xzzzz_yzz = pbuffer.data(idx_ovl_hf + 148);

    auto ts_xzzzz_zzz = pbuffer.data(idx_ovl_hf + 149);

    auto ts_yyyyy_xxx = pbuffer.data(idx_ovl_hf + 150);

    auto ts_yyyyy_xxy = pbuffer.data(idx_ovl_hf + 151);

    auto ts_yyyyy_xxz = pbuffer.data(idx_ovl_hf + 152);

    auto ts_yyyyy_xyy = pbuffer.data(idx_ovl_hf + 153);

    auto ts_yyyyy_xyz = pbuffer.data(idx_ovl_hf + 154);

    auto ts_yyyyy_xzz = pbuffer.data(idx_ovl_hf + 155);

    auto ts_yyyyy_yyy = pbuffer.data(idx_ovl_hf + 156);

    auto ts_yyyyy_yyz = pbuffer.data(idx_ovl_hf + 157);

    auto ts_yyyyy_yzz = pbuffer.data(idx_ovl_hf + 158);

    auto ts_yyyyy_zzz = pbuffer.data(idx_ovl_hf + 159);

    auto ts_yyyyz_xxy = pbuffer.data(idx_ovl_hf + 161);

    auto ts_yyyyz_xxz = pbuffer.data(idx_ovl_hf + 162);

    auto ts_yyyyz_xyy = pbuffer.data(idx_ovl_hf + 163);

    auto ts_yyyyz_xyz = pbuffer.data(idx_ovl_hf + 164);

    auto ts_yyyyz_xzz = pbuffer.data(idx_ovl_hf + 165);

    auto ts_yyyyz_yyy = pbuffer.data(idx_ovl_hf + 166);

    auto ts_yyyyz_yyz = pbuffer.data(idx_ovl_hf + 167);

    auto ts_yyyyz_yzz = pbuffer.data(idx_ovl_hf + 168);

    auto ts_yyyyz_zzz = pbuffer.data(idx_ovl_hf + 169);

    auto ts_yyyzz_xxx = pbuffer.data(idx_ovl_hf + 170);

    auto ts_yyyzz_xxy = pbuffer.data(idx_ovl_hf + 171);

    auto ts_yyyzz_xxz = pbuffer.data(idx_ovl_hf + 172);

    auto ts_yyyzz_xyy = pbuffer.data(idx_ovl_hf + 173);

    auto ts_yyyzz_xyz = pbuffer.data(idx_ovl_hf + 174);

    auto ts_yyyzz_xzz = pbuffer.data(idx_ovl_hf + 175);

    auto ts_yyyzz_yyy = pbuffer.data(idx_ovl_hf + 176);

    auto ts_yyyzz_yyz = pbuffer.data(idx_ovl_hf + 177);

    auto ts_yyyzz_yzz = pbuffer.data(idx_ovl_hf + 178);

    auto ts_yyyzz_zzz = pbuffer.data(idx_ovl_hf + 179);

    auto ts_yyzzz_xxx = pbuffer.data(idx_ovl_hf + 180);

    auto ts_yyzzz_xxy = pbuffer.data(idx_ovl_hf + 181);

    auto ts_yyzzz_xxz = pbuffer.data(idx_ovl_hf + 182);

    auto ts_yyzzz_xyy = pbuffer.data(idx_ovl_hf + 183);

    auto ts_yyzzz_xyz = pbuffer.data(idx_ovl_hf + 184);

    auto ts_yyzzz_xzz = pbuffer.data(idx_ovl_hf + 185);

    auto ts_yyzzz_yyy = pbuffer.data(idx_ovl_hf + 186);

    auto ts_yyzzz_yyz = pbuffer.data(idx_ovl_hf + 187);

    auto ts_yyzzz_yzz = pbuffer.data(idx_ovl_hf + 188);

    auto ts_yyzzz_zzz = pbuffer.data(idx_ovl_hf + 189);

    auto ts_yzzzz_xxx = pbuffer.data(idx_ovl_hf + 190);

    auto ts_yzzzz_xxy = pbuffer.data(idx_ovl_hf + 191);

    auto ts_yzzzz_xxz = pbuffer.data(idx_ovl_hf + 192);

    auto ts_yzzzz_xyy = pbuffer.data(idx_ovl_hf + 193);

    auto ts_yzzzz_xyz = pbuffer.data(idx_ovl_hf + 194);

    auto ts_yzzzz_xzz = pbuffer.data(idx_ovl_hf + 195);

    auto ts_yzzzz_yyy = pbuffer.data(idx_ovl_hf + 196);

    auto ts_yzzzz_yyz = pbuffer.data(idx_ovl_hf + 197);

    auto ts_yzzzz_yzz = pbuffer.data(idx_ovl_hf + 198);

    auto ts_yzzzz_zzz = pbuffer.data(idx_ovl_hf + 199);

    auto ts_zzzzz_xxx = pbuffer.data(idx_ovl_hf + 200);

    auto ts_zzzzz_xxy = pbuffer.data(idx_ovl_hf + 201);

    auto ts_zzzzz_xxz = pbuffer.data(idx_ovl_hf + 202);

    auto ts_zzzzz_xyy = pbuffer.data(idx_ovl_hf + 203);

    auto ts_zzzzz_xyz = pbuffer.data(idx_ovl_hf + 204);

    auto ts_zzzzz_xzz = pbuffer.data(idx_ovl_hf + 205);

    auto ts_zzzzz_yyy = pbuffer.data(idx_ovl_hf + 206);

    auto ts_zzzzz_yyz = pbuffer.data(idx_ovl_hf + 207);

    auto ts_zzzzz_yzz = pbuffer.data(idx_ovl_hf + 208);

    auto ts_zzzzz_zzz = pbuffer.data(idx_ovl_hf + 209);

    // Set up 0-10 components of targeted buffer : IF

    auto ts_xxxxxx_xxx = pbuffer.data(idx_ovl_if);

    auto ts_xxxxxx_xxy = pbuffer.data(idx_ovl_if + 1);

    auto ts_xxxxxx_xxz = pbuffer.data(idx_ovl_if + 2);

    auto ts_xxxxxx_xyy = pbuffer.data(idx_ovl_if + 3);

    auto ts_xxxxxx_xyz = pbuffer.data(idx_ovl_if + 4);

    auto ts_xxxxxx_xzz = pbuffer.data(idx_ovl_if + 5);

    auto ts_xxxxxx_yyy = pbuffer.data(idx_ovl_if + 6);

    auto ts_xxxxxx_yyz = pbuffer.data(idx_ovl_if + 7);

    auto ts_xxxxxx_yzz = pbuffer.data(idx_ovl_if + 8);

    auto ts_xxxxxx_zzz = pbuffer.data(idx_ovl_if + 9);

    #pragma omp simd aligned(pa_x, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xxz, ts_xxxx_xyy, ts_xxxx_xyz, ts_xxxx_xzz, ts_xxxx_yyy, ts_xxxx_yyz, ts_xxxx_yzz, ts_xxxx_zzz, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_yy, ts_xxxxx_yyy, ts_xxxxx_yyz, ts_xxxxx_yz, ts_xxxxx_yzz, ts_xxxxx_zz, ts_xxxxx_zzz, ts_xxxxxx_xxx, ts_xxxxxx_xxy, ts_xxxxxx_xxz, ts_xxxxxx_xyy, ts_xxxxxx_xyz, ts_xxxxxx_xzz, ts_xxxxxx_yyy, ts_xxxxxx_yyz, ts_xxxxxx_yzz, ts_xxxxxx_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxx_xxx[i] = 5.0 * ts_xxxx_xxx[i] * fe_0 + 3.0 * ts_xxxxx_xx[i] * fe_0 + ts_xxxxx_xxx[i] * pa_x[i];

        ts_xxxxxx_xxy[i] = 5.0 * ts_xxxx_xxy[i] * fe_0 + 2.0 * ts_xxxxx_xy[i] * fe_0 + ts_xxxxx_xxy[i] * pa_x[i];

        ts_xxxxxx_xxz[i] = 5.0 * ts_xxxx_xxz[i] * fe_0 + 2.0 * ts_xxxxx_xz[i] * fe_0 + ts_xxxxx_xxz[i] * pa_x[i];

        ts_xxxxxx_xyy[i] = 5.0 * ts_xxxx_xyy[i] * fe_0 + ts_xxxxx_yy[i] * fe_0 + ts_xxxxx_xyy[i] * pa_x[i];

        ts_xxxxxx_xyz[i] = 5.0 * ts_xxxx_xyz[i] * fe_0 + ts_xxxxx_yz[i] * fe_0 + ts_xxxxx_xyz[i] * pa_x[i];

        ts_xxxxxx_xzz[i] = 5.0 * ts_xxxx_xzz[i] * fe_0 + ts_xxxxx_zz[i] * fe_0 + ts_xxxxx_xzz[i] * pa_x[i];

        ts_xxxxxx_yyy[i] = 5.0 * ts_xxxx_yyy[i] * fe_0 + ts_xxxxx_yyy[i] * pa_x[i];

        ts_xxxxxx_yyz[i] = 5.0 * ts_xxxx_yyz[i] * fe_0 + ts_xxxxx_yyz[i] * pa_x[i];

        ts_xxxxxx_yzz[i] = 5.0 * ts_xxxx_yzz[i] * fe_0 + ts_xxxxx_yzz[i] * pa_x[i];

        ts_xxxxxx_zzz[i] = 5.0 * ts_xxxx_zzz[i] * fe_0 + ts_xxxxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : IF

    auto ts_xxxxxy_xxx = pbuffer.data(idx_ovl_if + 10);

    auto ts_xxxxxy_xxy = pbuffer.data(idx_ovl_if + 11);

    auto ts_xxxxxy_xxz = pbuffer.data(idx_ovl_if + 12);

    auto ts_xxxxxy_xyy = pbuffer.data(idx_ovl_if + 13);

    auto ts_xxxxxy_xyz = pbuffer.data(idx_ovl_if + 14);

    auto ts_xxxxxy_xzz = pbuffer.data(idx_ovl_if + 15);

    auto ts_xxxxxy_yyy = pbuffer.data(idx_ovl_if + 16);

    auto ts_xxxxxy_yyz = pbuffer.data(idx_ovl_if + 17);

    auto ts_xxxxxy_yzz = pbuffer.data(idx_ovl_if + 18);

    auto ts_xxxxxy_zzz = pbuffer.data(idx_ovl_if + 19);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_zzz, ts_xxxxxy_xxx, ts_xxxxxy_xxy, ts_xxxxxy_xxz, ts_xxxxxy_xyy, ts_xxxxxy_xyz, ts_xxxxxy_xzz, ts_xxxxxy_yyy, ts_xxxxxy_yyz, ts_xxxxxy_yzz, ts_xxxxxy_zzz, ts_xxxxy_yyy, ts_xxxxy_yyz, ts_xxxxy_yzz, ts_xxxy_yyy, ts_xxxy_yyz, ts_xxxy_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxy_xxx[i] = ts_xxxxx_xxx[i] * pa_y[i];

        ts_xxxxxy_xxy[i] = ts_xxxxx_xx[i] * fe_0 + ts_xxxxx_xxy[i] * pa_y[i];

        ts_xxxxxy_xxz[i] = ts_xxxxx_xxz[i] * pa_y[i];

        ts_xxxxxy_xyy[i] = 2.0 * ts_xxxxx_xy[i] * fe_0 + ts_xxxxx_xyy[i] * pa_y[i];

        ts_xxxxxy_xyz[i] = ts_xxxxx_xz[i] * fe_0 + ts_xxxxx_xyz[i] * pa_y[i];

        ts_xxxxxy_xzz[i] = ts_xxxxx_xzz[i] * pa_y[i];

        ts_xxxxxy_yyy[i] = 4.0 * ts_xxxy_yyy[i] * fe_0 + ts_xxxxy_yyy[i] * pa_x[i];

        ts_xxxxxy_yyz[i] = 4.0 * ts_xxxy_yyz[i] * fe_0 + ts_xxxxy_yyz[i] * pa_x[i];

        ts_xxxxxy_yzz[i] = 4.0 * ts_xxxy_yzz[i] * fe_0 + ts_xxxxy_yzz[i] * pa_x[i];

        ts_xxxxxy_zzz[i] = ts_xxxxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : IF

    auto ts_xxxxxz_xxx = pbuffer.data(idx_ovl_if + 20);

    auto ts_xxxxxz_xxy = pbuffer.data(idx_ovl_if + 21);

    auto ts_xxxxxz_xxz = pbuffer.data(idx_ovl_if + 22);

    auto ts_xxxxxz_xyy = pbuffer.data(idx_ovl_if + 23);

    auto ts_xxxxxz_xyz = pbuffer.data(idx_ovl_if + 24);

    auto ts_xxxxxz_xzz = pbuffer.data(idx_ovl_if + 25);

    auto ts_xxxxxz_yyy = pbuffer.data(idx_ovl_if + 26);

    auto ts_xxxxxz_yyz = pbuffer.data(idx_ovl_if + 27);

    auto ts_xxxxxz_yzz = pbuffer.data(idx_ovl_if + 28);

    auto ts_xxxxxz_zzz = pbuffer.data(idx_ovl_if + 29);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxxx_xx, ts_xxxxx_xxx, ts_xxxxx_xxy, ts_xxxxx_xxz, ts_xxxxx_xy, ts_xxxxx_xyy, ts_xxxxx_xyz, ts_xxxxx_xz, ts_xxxxx_xzz, ts_xxxxx_yyy, ts_xxxxxz_xxx, ts_xxxxxz_xxy, ts_xxxxxz_xxz, ts_xxxxxz_xyy, ts_xxxxxz_xyz, ts_xxxxxz_xzz, ts_xxxxxz_yyy, ts_xxxxxz_yyz, ts_xxxxxz_yzz, ts_xxxxxz_zzz, ts_xxxxz_yyz, ts_xxxxz_yzz, ts_xxxxz_zzz, ts_xxxz_yyz, ts_xxxz_yzz, ts_xxxz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxxz_xxx[i] = ts_xxxxx_xxx[i] * pa_z[i];

        ts_xxxxxz_xxy[i] = ts_xxxxx_xxy[i] * pa_z[i];

        ts_xxxxxz_xxz[i] = ts_xxxxx_xx[i] * fe_0 + ts_xxxxx_xxz[i] * pa_z[i];

        ts_xxxxxz_xyy[i] = ts_xxxxx_xyy[i] * pa_z[i];

        ts_xxxxxz_xyz[i] = ts_xxxxx_xy[i] * fe_0 + ts_xxxxx_xyz[i] * pa_z[i];

        ts_xxxxxz_xzz[i] = 2.0 * ts_xxxxx_xz[i] * fe_0 + ts_xxxxx_xzz[i] * pa_z[i];

        ts_xxxxxz_yyy[i] = ts_xxxxx_yyy[i] * pa_z[i];

        ts_xxxxxz_yyz[i] = 4.0 * ts_xxxz_yyz[i] * fe_0 + ts_xxxxz_yyz[i] * pa_x[i];

        ts_xxxxxz_yzz[i] = 4.0 * ts_xxxz_yzz[i] * fe_0 + ts_xxxxz_yzz[i] * pa_x[i];

        ts_xxxxxz_zzz[i] = 4.0 * ts_xxxz_zzz[i] * fe_0 + ts_xxxxz_zzz[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : IF

    auto ts_xxxxyy_xxx = pbuffer.data(idx_ovl_if + 30);

    auto ts_xxxxyy_xxy = pbuffer.data(idx_ovl_if + 31);

    auto ts_xxxxyy_xxz = pbuffer.data(idx_ovl_if + 32);

    auto ts_xxxxyy_xyy = pbuffer.data(idx_ovl_if + 33);

    auto ts_xxxxyy_xyz = pbuffer.data(idx_ovl_if + 34);

    auto ts_xxxxyy_xzz = pbuffer.data(idx_ovl_if + 35);

    auto ts_xxxxyy_yyy = pbuffer.data(idx_ovl_if + 36);

    auto ts_xxxxyy_yyz = pbuffer.data(idx_ovl_if + 37);

    auto ts_xxxxyy_yzz = pbuffer.data(idx_ovl_if + 38);

    auto ts_xxxxyy_zzz = pbuffer.data(idx_ovl_if + 39);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxx_xxx, ts_xxxx_xxz, ts_xxxx_xzz, ts_xxxxy_xxx, ts_xxxxy_xxz, ts_xxxxy_xzz, ts_xxxxyy_xxx, ts_xxxxyy_xxy, ts_xxxxyy_xxz, ts_xxxxyy_xyy, ts_xxxxyy_xyz, ts_xxxxyy_xzz, ts_xxxxyy_yyy, ts_xxxxyy_yyz, ts_xxxxyy_yzz, ts_xxxxyy_zzz, ts_xxxyy_xxy, ts_xxxyy_xy, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_yy, ts_xxxyy_yyy, ts_xxxyy_yyz, ts_xxxyy_yz, ts_xxxyy_yzz, ts_xxxyy_zzz, ts_xxyy_xxy, ts_xxyy_xyy, ts_xxyy_xyz, ts_xxyy_yyy, ts_xxyy_yyz, ts_xxyy_yzz, ts_xxyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyy_xxx[i] = ts_xxxx_xxx[i] * fe_0 + ts_xxxxy_xxx[i] * pa_y[i];

        ts_xxxxyy_xxy[i] = 3.0 * ts_xxyy_xxy[i] * fe_0 + 2.0 * ts_xxxyy_xy[i] * fe_0 + ts_xxxyy_xxy[i] * pa_x[i];

        ts_xxxxyy_xxz[i] = ts_xxxx_xxz[i] * fe_0 + ts_xxxxy_xxz[i] * pa_y[i];

        ts_xxxxyy_xyy[i] = 3.0 * ts_xxyy_xyy[i] * fe_0 + ts_xxxyy_yy[i] * fe_0 + ts_xxxyy_xyy[i] * pa_x[i];

        ts_xxxxyy_xyz[i] = 3.0 * ts_xxyy_xyz[i] * fe_0 + ts_xxxyy_yz[i] * fe_0 + ts_xxxyy_xyz[i] * pa_x[i];

        ts_xxxxyy_xzz[i] = ts_xxxx_xzz[i] * fe_0 + ts_xxxxy_xzz[i] * pa_y[i];

        ts_xxxxyy_yyy[i] = 3.0 * ts_xxyy_yyy[i] * fe_0 + ts_xxxyy_yyy[i] * pa_x[i];

        ts_xxxxyy_yyz[i] = 3.0 * ts_xxyy_yyz[i] * fe_0 + ts_xxxyy_yyz[i] * pa_x[i];

        ts_xxxxyy_yzz[i] = 3.0 * ts_xxyy_yzz[i] * fe_0 + ts_xxxyy_yzz[i] * pa_x[i];

        ts_xxxxyy_zzz[i] = 3.0 * ts_xxyy_zzz[i] * fe_0 + ts_xxxyy_zzz[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : IF

    auto ts_xxxxyz_xxx = pbuffer.data(idx_ovl_if + 40);

    auto ts_xxxxyz_xxy = pbuffer.data(idx_ovl_if + 41);

    auto ts_xxxxyz_xxz = pbuffer.data(idx_ovl_if + 42);

    auto ts_xxxxyz_xyy = pbuffer.data(idx_ovl_if + 43);

    auto ts_xxxxyz_xyz = pbuffer.data(idx_ovl_if + 44);

    auto ts_xxxxyz_xzz = pbuffer.data(idx_ovl_if + 45);

    auto ts_xxxxyz_yyy = pbuffer.data(idx_ovl_if + 46);

    auto ts_xxxxyz_yyz = pbuffer.data(idx_ovl_if + 47);

    auto ts_xxxxyz_yzz = pbuffer.data(idx_ovl_if + 48);

    auto ts_xxxxyz_zzz = pbuffer.data(idx_ovl_if + 49);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxxy_xxy, ts_xxxxy_xyy, ts_xxxxy_yyy, ts_xxxxyz_xxx, ts_xxxxyz_xxy, ts_xxxxyz_xxz, ts_xxxxyz_xyy, ts_xxxxyz_xyz, ts_xxxxyz_xzz, ts_xxxxyz_yyy, ts_xxxxyz_yyz, ts_xxxxyz_yzz, ts_xxxxyz_zzz, ts_xxxxz_xxx, ts_xxxxz_xxz, ts_xxxxz_xyz, ts_xxxxz_xz, ts_xxxxz_xzz, ts_xxxxz_zzz, ts_xxxyz_yyz, ts_xxxyz_yzz, ts_xxyz_yyz, ts_xxyz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxyz_xxx[i] = ts_xxxxz_xxx[i] * pa_y[i];

        ts_xxxxyz_xxy[i] = ts_xxxxy_xxy[i] * pa_z[i];

        ts_xxxxyz_xxz[i] = ts_xxxxz_xxz[i] * pa_y[i];

        ts_xxxxyz_xyy[i] = ts_xxxxy_xyy[i] * pa_z[i];

        ts_xxxxyz_xyz[i] = ts_xxxxz_xz[i] * fe_0 + ts_xxxxz_xyz[i] * pa_y[i];

        ts_xxxxyz_xzz[i] = ts_xxxxz_xzz[i] * pa_y[i];

        ts_xxxxyz_yyy[i] = ts_xxxxy_yyy[i] * pa_z[i];

        ts_xxxxyz_yyz[i] = 3.0 * ts_xxyz_yyz[i] * fe_0 + ts_xxxyz_yyz[i] * pa_x[i];

        ts_xxxxyz_yzz[i] = 3.0 * ts_xxyz_yzz[i] * fe_0 + ts_xxxyz_yzz[i] * pa_x[i];

        ts_xxxxyz_zzz[i] = ts_xxxxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : IF

    auto ts_xxxxzz_xxx = pbuffer.data(idx_ovl_if + 50);

    auto ts_xxxxzz_xxy = pbuffer.data(idx_ovl_if + 51);

    auto ts_xxxxzz_xxz = pbuffer.data(idx_ovl_if + 52);

    auto ts_xxxxzz_xyy = pbuffer.data(idx_ovl_if + 53);

    auto ts_xxxxzz_xyz = pbuffer.data(idx_ovl_if + 54);

    auto ts_xxxxzz_xzz = pbuffer.data(idx_ovl_if + 55);

    auto ts_xxxxzz_yyy = pbuffer.data(idx_ovl_if + 56);

    auto ts_xxxxzz_yyz = pbuffer.data(idx_ovl_if + 57);

    auto ts_xxxxzz_yzz = pbuffer.data(idx_ovl_if + 58);

    auto ts_xxxxzz_zzz = pbuffer.data(idx_ovl_if + 59);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxx_xxx, ts_xxxx_xxy, ts_xxxx_xyy, ts_xxxxz_xxx, ts_xxxxz_xxy, ts_xxxxz_xyy, ts_xxxxzz_xxx, ts_xxxxzz_xxy, ts_xxxxzz_xxz, ts_xxxxzz_xyy, ts_xxxxzz_xyz, ts_xxxxzz_xzz, ts_xxxxzz_yyy, ts_xxxxzz_yyz, ts_xxxxzz_yzz, ts_xxxxzz_zzz, ts_xxxzz_xxz, ts_xxxzz_xyz, ts_xxxzz_xz, ts_xxxzz_xzz, ts_xxxzz_yyy, ts_xxxzz_yyz, ts_xxxzz_yz, ts_xxxzz_yzz, ts_xxxzz_zz, ts_xxxzz_zzz, ts_xxzz_xxz, ts_xxzz_xyz, ts_xxzz_xzz, ts_xxzz_yyy, ts_xxzz_yyz, ts_xxzz_yzz, ts_xxzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxxzz_xxx[i] = ts_xxxx_xxx[i] * fe_0 + ts_xxxxz_xxx[i] * pa_z[i];

        ts_xxxxzz_xxy[i] = ts_xxxx_xxy[i] * fe_0 + ts_xxxxz_xxy[i] * pa_z[i];

        ts_xxxxzz_xxz[i] = 3.0 * ts_xxzz_xxz[i] * fe_0 + 2.0 * ts_xxxzz_xz[i] * fe_0 + ts_xxxzz_xxz[i] * pa_x[i];

        ts_xxxxzz_xyy[i] = ts_xxxx_xyy[i] * fe_0 + ts_xxxxz_xyy[i] * pa_z[i];

        ts_xxxxzz_xyz[i] = 3.0 * ts_xxzz_xyz[i] * fe_0 + ts_xxxzz_yz[i] * fe_0 + ts_xxxzz_xyz[i] * pa_x[i];

        ts_xxxxzz_xzz[i] = 3.0 * ts_xxzz_xzz[i] * fe_0 + ts_xxxzz_zz[i] * fe_0 + ts_xxxzz_xzz[i] * pa_x[i];

        ts_xxxxzz_yyy[i] = 3.0 * ts_xxzz_yyy[i] * fe_0 + ts_xxxzz_yyy[i] * pa_x[i];

        ts_xxxxzz_yyz[i] = 3.0 * ts_xxzz_yyz[i] * fe_0 + ts_xxxzz_yyz[i] * pa_x[i];

        ts_xxxxzz_yzz[i] = 3.0 * ts_xxzz_yzz[i] * fe_0 + ts_xxxzz_yzz[i] * pa_x[i];

        ts_xxxxzz_zzz[i] = 3.0 * ts_xxzz_zzz[i] * fe_0 + ts_xxxzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : IF

    auto ts_xxxyyy_xxx = pbuffer.data(idx_ovl_if + 60);

    auto ts_xxxyyy_xxy = pbuffer.data(idx_ovl_if + 61);

    auto ts_xxxyyy_xxz = pbuffer.data(idx_ovl_if + 62);

    auto ts_xxxyyy_xyy = pbuffer.data(idx_ovl_if + 63);

    auto ts_xxxyyy_xyz = pbuffer.data(idx_ovl_if + 64);

    auto ts_xxxyyy_xzz = pbuffer.data(idx_ovl_if + 65);

    auto ts_xxxyyy_yyy = pbuffer.data(idx_ovl_if + 66);

    auto ts_xxxyyy_yyz = pbuffer.data(idx_ovl_if + 67);

    auto ts_xxxyyy_yzz = pbuffer.data(idx_ovl_if + 68);

    auto ts_xxxyyy_zzz = pbuffer.data(idx_ovl_if + 69);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxy_xxx, ts_xxxy_xxz, ts_xxxy_xzz, ts_xxxyy_xxx, ts_xxxyy_xxz, ts_xxxyy_xzz, ts_xxxyyy_xxx, ts_xxxyyy_xxy, ts_xxxyyy_xxz, ts_xxxyyy_xyy, ts_xxxyyy_xyz, ts_xxxyyy_xzz, ts_xxxyyy_yyy, ts_xxxyyy_yyz, ts_xxxyyy_yzz, ts_xxxyyy_zzz, ts_xxyyy_xxy, ts_xxyyy_xy, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_yy, ts_xxyyy_yyy, ts_xxyyy_yyz, ts_xxyyy_yz, ts_xxyyy_yzz, ts_xxyyy_zzz, ts_xyyy_xxy, ts_xyyy_xyy, ts_xyyy_xyz, ts_xyyy_yyy, ts_xyyy_yyz, ts_xyyy_yzz, ts_xyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyy_xxx[i] = 2.0 * ts_xxxy_xxx[i] * fe_0 + ts_xxxyy_xxx[i] * pa_y[i];

        ts_xxxyyy_xxy[i] = 2.0 * ts_xyyy_xxy[i] * fe_0 + 2.0 * ts_xxyyy_xy[i] * fe_0 + ts_xxyyy_xxy[i] * pa_x[i];

        ts_xxxyyy_xxz[i] = 2.0 * ts_xxxy_xxz[i] * fe_0 + ts_xxxyy_xxz[i] * pa_y[i];

        ts_xxxyyy_xyy[i] = 2.0 * ts_xyyy_xyy[i] * fe_0 + ts_xxyyy_yy[i] * fe_0 + ts_xxyyy_xyy[i] * pa_x[i];

        ts_xxxyyy_xyz[i] = 2.0 * ts_xyyy_xyz[i] * fe_0 + ts_xxyyy_yz[i] * fe_0 + ts_xxyyy_xyz[i] * pa_x[i];

        ts_xxxyyy_xzz[i] = 2.0 * ts_xxxy_xzz[i] * fe_0 + ts_xxxyy_xzz[i] * pa_y[i];

        ts_xxxyyy_yyy[i] = 2.0 * ts_xyyy_yyy[i] * fe_0 + ts_xxyyy_yyy[i] * pa_x[i];

        ts_xxxyyy_yyz[i] = 2.0 * ts_xyyy_yyz[i] * fe_0 + ts_xxyyy_yyz[i] * pa_x[i];

        ts_xxxyyy_yzz[i] = 2.0 * ts_xyyy_yzz[i] * fe_0 + ts_xxyyy_yzz[i] * pa_x[i];

        ts_xxxyyy_zzz[i] = 2.0 * ts_xyyy_zzz[i] * fe_0 + ts_xxyyy_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : IF

    auto ts_xxxyyz_xxx = pbuffer.data(idx_ovl_if + 70);

    auto ts_xxxyyz_xxy = pbuffer.data(idx_ovl_if + 71);

    auto ts_xxxyyz_xxz = pbuffer.data(idx_ovl_if + 72);

    auto ts_xxxyyz_xyy = pbuffer.data(idx_ovl_if + 73);

    auto ts_xxxyyz_xyz = pbuffer.data(idx_ovl_if + 74);

    auto ts_xxxyyz_xzz = pbuffer.data(idx_ovl_if + 75);

    auto ts_xxxyyz_yyy = pbuffer.data(idx_ovl_if + 76);

    auto ts_xxxyyz_yyz = pbuffer.data(idx_ovl_if + 77);

    auto ts_xxxyyz_yzz = pbuffer.data(idx_ovl_if + 78);

    auto ts_xxxyyz_zzz = pbuffer.data(idx_ovl_if + 79);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxxyy_xxx, ts_xxxyy_xxy, ts_xxxyy_xy, ts_xxxyy_xyy, ts_xxxyy_xyz, ts_xxxyy_yyy, ts_xxxyyz_xxx, ts_xxxyyz_xxy, ts_xxxyyz_xxz, ts_xxxyyz_xyy, ts_xxxyyz_xyz, ts_xxxyyz_xzz, ts_xxxyyz_yyy, ts_xxxyyz_yyz, ts_xxxyyz_yzz, ts_xxxyyz_zzz, ts_xxxyz_xxz, ts_xxxyz_xzz, ts_xxxz_xxz, ts_xxxz_xzz, ts_xxyyz_yyz, ts_xxyyz_yzz, ts_xxyyz_zzz, ts_xyyz_yyz, ts_xyyz_yzz, ts_xyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyyz_xxx[i] = ts_xxxyy_xxx[i] * pa_z[i];

        ts_xxxyyz_xxy[i] = ts_xxxyy_xxy[i] * pa_z[i];

        ts_xxxyyz_xxz[i] = ts_xxxz_xxz[i] * fe_0 + ts_xxxyz_xxz[i] * pa_y[i];

        ts_xxxyyz_xyy[i] = ts_xxxyy_xyy[i] * pa_z[i];

        ts_xxxyyz_xyz[i] = ts_xxxyy_xy[i] * fe_0 + ts_xxxyy_xyz[i] * pa_z[i];

        ts_xxxyyz_xzz[i] = ts_xxxz_xzz[i] * fe_0 + ts_xxxyz_xzz[i] * pa_y[i];

        ts_xxxyyz_yyy[i] = ts_xxxyy_yyy[i] * pa_z[i];

        ts_xxxyyz_yyz[i] = 2.0 * ts_xyyz_yyz[i] * fe_0 + ts_xxyyz_yyz[i] * pa_x[i];

        ts_xxxyyz_yzz[i] = 2.0 * ts_xyyz_yzz[i] * fe_0 + ts_xxyyz_yzz[i] * pa_x[i];

        ts_xxxyyz_zzz[i] = 2.0 * ts_xyyz_zzz[i] * fe_0 + ts_xxyyz_zzz[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : IF

    auto ts_xxxyzz_xxx = pbuffer.data(idx_ovl_if + 80);

    auto ts_xxxyzz_xxy = pbuffer.data(idx_ovl_if + 81);

    auto ts_xxxyzz_xxz = pbuffer.data(idx_ovl_if + 82);

    auto ts_xxxyzz_xyy = pbuffer.data(idx_ovl_if + 83);

    auto ts_xxxyzz_xyz = pbuffer.data(idx_ovl_if + 84);

    auto ts_xxxyzz_xzz = pbuffer.data(idx_ovl_if + 85);

    auto ts_xxxyzz_yyy = pbuffer.data(idx_ovl_if + 86);

    auto ts_xxxyzz_yyz = pbuffer.data(idx_ovl_if + 87);

    auto ts_xxxyzz_yzz = pbuffer.data(idx_ovl_if + 88);

    auto ts_xxxyzz_zzz = pbuffer.data(idx_ovl_if + 89);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxxyzz_xxx, ts_xxxyzz_xxy, ts_xxxyzz_xxz, ts_xxxyzz_xyy, ts_xxxyzz_xyz, ts_xxxyzz_xzz, ts_xxxyzz_yyy, ts_xxxyzz_yyz, ts_xxxyzz_yzz, ts_xxxyzz_zzz, ts_xxxzz_xx, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xxz, ts_xxxzz_xy, ts_xxxzz_xyy, ts_xxxzz_xyz, ts_xxxzz_xz, ts_xxxzz_xzz, ts_xxxzz_zzz, ts_xxyzz_yyy, ts_xxyzz_yyz, ts_xxyzz_yzz, ts_xyzz_yyy, ts_xyzz_yyz, ts_xyzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxyzz_xxx[i] = ts_xxxzz_xxx[i] * pa_y[i];

        ts_xxxyzz_xxy[i] = ts_xxxzz_xx[i] * fe_0 + ts_xxxzz_xxy[i] * pa_y[i];

        ts_xxxyzz_xxz[i] = ts_xxxzz_xxz[i] * pa_y[i];

        ts_xxxyzz_xyy[i] = 2.0 * ts_xxxzz_xy[i] * fe_0 + ts_xxxzz_xyy[i] * pa_y[i];

        ts_xxxyzz_xyz[i] = ts_xxxzz_xz[i] * fe_0 + ts_xxxzz_xyz[i] * pa_y[i];

        ts_xxxyzz_xzz[i] = ts_xxxzz_xzz[i] * pa_y[i];

        ts_xxxyzz_yyy[i] = 2.0 * ts_xyzz_yyy[i] * fe_0 + ts_xxyzz_yyy[i] * pa_x[i];

        ts_xxxyzz_yyz[i] = 2.0 * ts_xyzz_yyz[i] * fe_0 + ts_xxyzz_yyz[i] * pa_x[i];

        ts_xxxyzz_yzz[i] = 2.0 * ts_xyzz_yzz[i] * fe_0 + ts_xxyzz_yzz[i] * pa_x[i];

        ts_xxxyzz_zzz[i] = ts_xxxzz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : IF

    auto ts_xxxzzz_xxx = pbuffer.data(idx_ovl_if + 90);

    auto ts_xxxzzz_xxy = pbuffer.data(idx_ovl_if + 91);

    auto ts_xxxzzz_xxz = pbuffer.data(idx_ovl_if + 92);

    auto ts_xxxzzz_xyy = pbuffer.data(idx_ovl_if + 93);

    auto ts_xxxzzz_xyz = pbuffer.data(idx_ovl_if + 94);

    auto ts_xxxzzz_xzz = pbuffer.data(idx_ovl_if + 95);

    auto ts_xxxzzz_yyy = pbuffer.data(idx_ovl_if + 96);

    auto ts_xxxzzz_yyz = pbuffer.data(idx_ovl_if + 97);

    auto ts_xxxzzz_yzz = pbuffer.data(idx_ovl_if + 98);

    auto ts_xxxzzz_zzz = pbuffer.data(idx_ovl_if + 99);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxxz_xxx, ts_xxxz_xxy, ts_xxxz_xyy, ts_xxxzz_xxx, ts_xxxzz_xxy, ts_xxxzz_xyy, ts_xxxzzz_xxx, ts_xxxzzz_xxy, ts_xxxzzz_xxz, ts_xxxzzz_xyy, ts_xxxzzz_xyz, ts_xxxzzz_xzz, ts_xxxzzz_yyy, ts_xxxzzz_yyz, ts_xxxzzz_yzz, ts_xxxzzz_zzz, ts_xxzzz_xxz, ts_xxzzz_xyz, ts_xxzzz_xz, ts_xxzzz_xzz, ts_xxzzz_yyy, ts_xxzzz_yyz, ts_xxzzz_yz, ts_xxzzz_yzz, ts_xxzzz_zz, ts_xxzzz_zzz, ts_xzzz_xxz, ts_xzzz_xyz, ts_xzzz_xzz, ts_xzzz_yyy, ts_xzzz_yyz, ts_xzzz_yzz, ts_xzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxxzzz_xxx[i] = 2.0 * ts_xxxz_xxx[i] * fe_0 + ts_xxxzz_xxx[i] * pa_z[i];

        ts_xxxzzz_xxy[i] = 2.0 * ts_xxxz_xxy[i] * fe_0 + ts_xxxzz_xxy[i] * pa_z[i];

        ts_xxxzzz_xxz[i] = 2.0 * ts_xzzz_xxz[i] * fe_0 + 2.0 * ts_xxzzz_xz[i] * fe_0 + ts_xxzzz_xxz[i] * pa_x[i];

        ts_xxxzzz_xyy[i] = 2.0 * ts_xxxz_xyy[i] * fe_0 + ts_xxxzz_xyy[i] * pa_z[i];

        ts_xxxzzz_xyz[i] = 2.0 * ts_xzzz_xyz[i] * fe_0 + ts_xxzzz_yz[i] * fe_0 + ts_xxzzz_xyz[i] * pa_x[i];

        ts_xxxzzz_xzz[i] = 2.0 * ts_xzzz_xzz[i] * fe_0 + ts_xxzzz_zz[i] * fe_0 + ts_xxzzz_xzz[i] * pa_x[i];

        ts_xxxzzz_yyy[i] = 2.0 * ts_xzzz_yyy[i] * fe_0 + ts_xxzzz_yyy[i] * pa_x[i];

        ts_xxxzzz_yyz[i] = 2.0 * ts_xzzz_yyz[i] * fe_0 + ts_xxzzz_yyz[i] * pa_x[i];

        ts_xxxzzz_yzz[i] = 2.0 * ts_xzzz_yzz[i] * fe_0 + ts_xxzzz_yzz[i] * pa_x[i];

        ts_xxxzzz_zzz[i] = 2.0 * ts_xzzz_zzz[i] * fe_0 + ts_xxzzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : IF

    auto ts_xxyyyy_xxx = pbuffer.data(idx_ovl_if + 100);

    auto ts_xxyyyy_xxy = pbuffer.data(idx_ovl_if + 101);

    auto ts_xxyyyy_xxz = pbuffer.data(idx_ovl_if + 102);

    auto ts_xxyyyy_xyy = pbuffer.data(idx_ovl_if + 103);

    auto ts_xxyyyy_xyz = pbuffer.data(idx_ovl_if + 104);

    auto ts_xxyyyy_xzz = pbuffer.data(idx_ovl_if + 105);

    auto ts_xxyyyy_yyy = pbuffer.data(idx_ovl_if + 106);

    auto ts_xxyyyy_yyz = pbuffer.data(idx_ovl_if + 107);

    auto ts_xxyyyy_yzz = pbuffer.data(idx_ovl_if + 108);

    auto ts_xxyyyy_zzz = pbuffer.data(idx_ovl_if + 109);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyy_xxx, ts_xxyy_xxz, ts_xxyy_xzz, ts_xxyyy_xxx, ts_xxyyy_xxz, ts_xxyyy_xzz, ts_xxyyyy_xxx, ts_xxyyyy_xxy, ts_xxyyyy_xxz, ts_xxyyyy_xyy, ts_xxyyyy_xyz, ts_xxyyyy_xzz, ts_xxyyyy_yyy, ts_xxyyyy_yyz, ts_xxyyyy_yzz, ts_xxyyyy_zzz, ts_xyyyy_xxy, ts_xyyyy_xy, ts_xyyyy_xyy, ts_xyyyy_xyz, ts_xyyyy_yy, ts_xyyyy_yyy, ts_xyyyy_yyz, ts_xyyyy_yz, ts_xyyyy_yzz, ts_xyyyy_zzz, ts_yyyy_xxy, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yzz, ts_yyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyy_xxx[i] = 3.0 * ts_xxyy_xxx[i] * fe_0 + ts_xxyyy_xxx[i] * pa_y[i];

        ts_xxyyyy_xxy[i] = ts_yyyy_xxy[i] * fe_0 + 2.0 * ts_xyyyy_xy[i] * fe_0 + ts_xyyyy_xxy[i] * pa_x[i];

        ts_xxyyyy_xxz[i] = 3.0 * ts_xxyy_xxz[i] * fe_0 + ts_xxyyy_xxz[i] * pa_y[i];

        ts_xxyyyy_xyy[i] = ts_yyyy_xyy[i] * fe_0 + ts_xyyyy_yy[i] * fe_0 + ts_xyyyy_xyy[i] * pa_x[i];

        ts_xxyyyy_xyz[i] = ts_yyyy_xyz[i] * fe_0 + ts_xyyyy_yz[i] * fe_0 + ts_xyyyy_xyz[i] * pa_x[i];

        ts_xxyyyy_xzz[i] = 3.0 * ts_xxyy_xzz[i] * fe_0 + ts_xxyyy_xzz[i] * pa_y[i];

        ts_xxyyyy_yyy[i] = ts_yyyy_yyy[i] * fe_0 + ts_xyyyy_yyy[i] * pa_x[i];

        ts_xxyyyy_yyz[i] = ts_yyyy_yyz[i] * fe_0 + ts_xyyyy_yyz[i] * pa_x[i];

        ts_xxyyyy_yzz[i] = ts_yyyy_yzz[i] * fe_0 + ts_xyyyy_yzz[i] * pa_x[i];

        ts_xxyyyy_zzz[i] = ts_yyyy_zzz[i] * fe_0 + ts_xyyyy_zzz[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : IF

    auto ts_xxyyyz_xxx = pbuffer.data(idx_ovl_if + 110);

    auto ts_xxyyyz_xxy = pbuffer.data(idx_ovl_if + 111);

    auto ts_xxyyyz_xxz = pbuffer.data(idx_ovl_if + 112);

    auto ts_xxyyyz_xyy = pbuffer.data(idx_ovl_if + 113);

    auto ts_xxyyyz_xyz = pbuffer.data(idx_ovl_if + 114);

    auto ts_xxyyyz_xzz = pbuffer.data(idx_ovl_if + 115);

    auto ts_xxyyyz_yyy = pbuffer.data(idx_ovl_if + 116);

    auto ts_xxyyyz_yyz = pbuffer.data(idx_ovl_if + 117);

    auto ts_xxyyyz_yzz = pbuffer.data(idx_ovl_if + 118);

    auto ts_xxyyyz_zzz = pbuffer.data(idx_ovl_if + 119);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyyy_xxx, ts_xxyyy_xxy, ts_xxyyy_xy, ts_xxyyy_xyy, ts_xxyyy_xyz, ts_xxyyy_yyy, ts_xxyyyz_xxx, ts_xxyyyz_xxy, ts_xxyyyz_xxz, ts_xxyyyz_xyy, ts_xxyyyz_xyz, ts_xxyyyz_xzz, ts_xxyyyz_yyy, ts_xxyyyz_yyz, ts_xxyyyz_yzz, ts_xxyyyz_zzz, ts_xxyyz_xxz, ts_xxyyz_xzz, ts_xxyz_xxz, ts_xxyz_xzz, ts_xyyyz_yyz, ts_xyyyz_yzz, ts_xyyyz_zzz, ts_yyyz_yyz, ts_yyyz_yzz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyyz_xxx[i] = ts_xxyyy_xxx[i] * pa_z[i];

        ts_xxyyyz_xxy[i] = ts_xxyyy_xxy[i] * pa_z[i];

        ts_xxyyyz_xxz[i] = 2.0 * ts_xxyz_xxz[i] * fe_0 + ts_xxyyz_xxz[i] * pa_y[i];

        ts_xxyyyz_xyy[i] = ts_xxyyy_xyy[i] * pa_z[i];

        ts_xxyyyz_xyz[i] = ts_xxyyy_xy[i] * fe_0 + ts_xxyyy_xyz[i] * pa_z[i];

        ts_xxyyyz_xzz[i] = 2.0 * ts_xxyz_xzz[i] * fe_0 + ts_xxyyz_xzz[i] * pa_y[i];

        ts_xxyyyz_yyy[i] = ts_xxyyy_yyy[i] * pa_z[i];

        ts_xxyyyz_yyz[i] = ts_yyyz_yyz[i] * fe_0 + ts_xyyyz_yyz[i] * pa_x[i];

        ts_xxyyyz_yzz[i] = ts_yyyz_yzz[i] * fe_0 + ts_xyyyz_yzz[i] * pa_x[i];

        ts_xxyyyz_zzz[i] = ts_yyyz_zzz[i] * fe_0 + ts_xyyyz_zzz[i] * pa_x[i];
    }

    // Set up 120-130 components of targeted buffer : IF

    auto ts_xxyyzz_xxx = pbuffer.data(idx_ovl_if + 120);

    auto ts_xxyyzz_xxy = pbuffer.data(idx_ovl_if + 121);

    auto ts_xxyyzz_xxz = pbuffer.data(idx_ovl_if + 122);

    auto ts_xxyyzz_xyy = pbuffer.data(idx_ovl_if + 123);

    auto ts_xxyyzz_xyz = pbuffer.data(idx_ovl_if + 124);

    auto ts_xxyyzz_xzz = pbuffer.data(idx_ovl_if + 125);

    auto ts_xxyyzz_yyy = pbuffer.data(idx_ovl_if + 126);

    auto ts_xxyyzz_yyz = pbuffer.data(idx_ovl_if + 127);

    auto ts_xxyyzz_yzz = pbuffer.data(idx_ovl_if + 128);

    auto ts_xxyyzz_zzz = pbuffer.data(idx_ovl_if + 129);

    #pragma omp simd aligned(pa_x, pa_y, pa_z, ts_xxyy_xxy, ts_xxyy_xyy, ts_xxyyz_xxy, ts_xxyyz_xyy, ts_xxyyzz_xxx, ts_xxyyzz_xxy, ts_xxyyzz_xxz, ts_xxyyzz_xyy, ts_xxyyzz_xyz, ts_xxyyzz_xzz, ts_xxyyzz_yyy, ts_xxyyzz_yyz, ts_xxyyzz_yzz, ts_xxyyzz_zzz, ts_xxyzz_xxx, ts_xxyzz_xxz, ts_xxyzz_xzz, ts_xxzz_xxx, ts_xxzz_xxz, ts_xxzz_xzz, ts_xyyzz_xyz, ts_xyyzz_yyy, ts_xyyzz_yyz, ts_xyyzz_yz, ts_xyyzz_yzz, ts_xyyzz_zzz, ts_yyzz_xyz, ts_yyzz_yyy, ts_yyzz_yyz, ts_yyzz_yzz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyyzz_xxx[i] = ts_xxzz_xxx[i] * fe_0 + ts_xxyzz_xxx[i] * pa_y[i];

        ts_xxyyzz_xxy[i] = ts_xxyy_xxy[i] * fe_0 + ts_xxyyz_xxy[i] * pa_z[i];

        ts_xxyyzz_xxz[i] = ts_xxzz_xxz[i] * fe_0 + ts_xxyzz_xxz[i] * pa_y[i];

        ts_xxyyzz_xyy[i] = ts_xxyy_xyy[i] * fe_0 + ts_xxyyz_xyy[i] * pa_z[i];

        ts_xxyyzz_xyz[i] = ts_yyzz_xyz[i] * fe_0 + ts_xyyzz_yz[i] * fe_0 + ts_xyyzz_xyz[i] * pa_x[i];

        ts_xxyyzz_xzz[i] = ts_xxzz_xzz[i] * fe_0 + ts_xxyzz_xzz[i] * pa_y[i];

        ts_xxyyzz_yyy[i] = ts_yyzz_yyy[i] * fe_0 + ts_xyyzz_yyy[i] * pa_x[i];

        ts_xxyyzz_yyz[i] = ts_yyzz_yyz[i] * fe_0 + ts_xyyzz_yyz[i] * pa_x[i];

        ts_xxyyzz_yzz[i] = ts_yyzz_yzz[i] * fe_0 + ts_xyyzz_yzz[i] * pa_x[i];

        ts_xxyyzz_zzz[i] = ts_yyzz_zzz[i] * fe_0 + ts_xyyzz_zzz[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : IF

    auto ts_xxyzzz_xxx = pbuffer.data(idx_ovl_if + 130);

    auto ts_xxyzzz_xxy = pbuffer.data(idx_ovl_if + 131);

    auto ts_xxyzzz_xxz = pbuffer.data(idx_ovl_if + 132);

    auto ts_xxyzzz_xyy = pbuffer.data(idx_ovl_if + 133);

    auto ts_xxyzzz_xyz = pbuffer.data(idx_ovl_if + 134);

    auto ts_xxyzzz_xzz = pbuffer.data(idx_ovl_if + 135);

    auto ts_xxyzzz_yyy = pbuffer.data(idx_ovl_if + 136);

    auto ts_xxyzzz_yyz = pbuffer.data(idx_ovl_if + 137);

    auto ts_xxyzzz_yzz = pbuffer.data(idx_ovl_if + 138);

    auto ts_xxyzzz_zzz = pbuffer.data(idx_ovl_if + 139);

    #pragma omp simd aligned(pa_x, pa_y, ts_xxyzzz_xxx, ts_xxyzzz_xxy, ts_xxyzzz_xxz, ts_xxyzzz_xyy, ts_xxyzzz_xyz, ts_xxyzzz_xzz, ts_xxyzzz_yyy, ts_xxyzzz_yyz, ts_xxyzzz_yzz, ts_xxyzzz_zzz, ts_xxzzz_xx, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xxz, ts_xxzzz_xy, ts_xxzzz_xyy, ts_xxzzz_xyz, ts_xxzzz_xz, ts_xxzzz_xzz, ts_xxzzz_zzz, ts_xyzzz_yyy, ts_xyzzz_yyz, ts_xyzzz_yzz, ts_yzzz_yyy, ts_yzzz_yyz, ts_yzzz_yzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxyzzz_xxx[i] = ts_xxzzz_xxx[i] * pa_y[i];

        ts_xxyzzz_xxy[i] = ts_xxzzz_xx[i] * fe_0 + ts_xxzzz_xxy[i] * pa_y[i];

        ts_xxyzzz_xxz[i] = ts_xxzzz_xxz[i] * pa_y[i];

        ts_xxyzzz_xyy[i] = 2.0 * ts_xxzzz_xy[i] * fe_0 + ts_xxzzz_xyy[i] * pa_y[i];

        ts_xxyzzz_xyz[i] = ts_xxzzz_xz[i] * fe_0 + ts_xxzzz_xyz[i] * pa_y[i];

        ts_xxyzzz_xzz[i] = ts_xxzzz_xzz[i] * pa_y[i];

        ts_xxyzzz_yyy[i] = ts_yzzz_yyy[i] * fe_0 + ts_xyzzz_yyy[i] * pa_x[i];

        ts_xxyzzz_yyz[i] = ts_yzzz_yyz[i] * fe_0 + ts_xyzzz_yyz[i] * pa_x[i];

        ts_xxyzzz_yzz[i] = ts_yzzz_yzz[i] * fe_0 + ts_xyzzz_yzz[i] * pa_x[i];

        ts_xxyzzz_zzz[i] = ts_xxzzz_zzz[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : IF

    auto ts_xxzzzz_xxx = pbuffer.data(idx_ovl_if + 140);

    auto ts_xxzzzz_xxy = pbuffer.data(idx_ovl_if + 141);

    auto ts_xxzzzz_xxz = pbuffer.data(idx_ovl_if + 142);

    auto ts_xxzzzz_xyy = pbuffer.data(idx_ovl_if + 143);

    auto ts_xxzzzz_xyz = pbuffer.data(idx_ovl_if + 144);

    auto ts_xxzzzz_xzz = pbuffer.data(idx_ovl_if + 145);

    auto ts_xxzzzz_yyy = pbuffer.data(idx_ovl_if + 146);

    auto ts_xxzzzz_yyz = pbuffer.data(idx_ovl_if + 147);

    auto ts_xxzzzz_yzz = pbuffer.data(idx_ovl_if + 148);

    auto ts_xxzzzz_zzz = pbuffer.data(idx_ovl_if + 149);

    #pragma omp simd aligned(pa_x, pa_z, ts_xxzz_xxx, ts_xxzz_xxy, ts_xxzz_xyy, ts_xxzzz_xxx, ts_xxzzz_xxy, ts_xxzzz_xyy, ts_xxzzzz_xxx, ts_xxzzzz_xxy, ts_xxzzzz_xxz, ts_xxzzzz_xyy, ts_xxzzzz_xyz, ts_xxzzzz_xzz, ts_xxzzzz_yyy, ts_xxzzzz_yyz, ts_xxzzzz_yzz, ts_xxzzzz_zzz, ts_xzzzz_xxz, ts_xzzzz_xyz, ts_xzzzz_xz, ts_xzzzz_xzz, ts_xzzzz_yyy, ts_xzzzz_yyz, ts_xzzzz_yz, ts_xzzzz_yzz, ts_xzzzz_zz, ts_xzzzz_zzz, ts_zzzz_xxz, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xxzzzz_xxx[i] = 3.0 * ts_xxzz_xxx[i] * fe_0 + ts_xxzzz_xxx[i] * pa_z[i];

        ts_xxzzzz_xxy[i] = 3.0 * ts_xxzz_xxy[i] * fe_0 + ts_xxzzz_xxy[i] * pa_z[i];

        ts_xxzzzz_xxz[i] = ts_zzzz_xxz[i] * fe_0 + 2.0 * ts_xzzzz_xz[i] * fe_0 + ts_xzzzz_xxz[i] * pa_x[i];

        ts_xxzzzz_xyy[i] = 3.0 * ts_xxzz_xyy[i] * fe_0 + ts_xxzzz_xyy[i] * pa_z[i];

        ts_xxzzzz_xyz[i] = ts_zzzz_xyz[i] * fe_0 + ts_xzzzz_yz[i] * fe_0 + ts_xzzzz_xyz[i] * pa_x[i];

        ts_xxzzzz_xzz[i] = ts_zzzz_xzz[i] * fe_0 + ts_xzzzz_zz[i] * fe_0 + ts_xzzzz_xzz[i] * pa_x[i];

        ts_xxzzzz_yyy[i] = ts_zzzz_yyy[i] * fe_0 + ts_xzzzz_yyy[i] * pa_x[i];

        ts_xxzzzz_yyz[i] = ts_zzzz_yyz[i] * fe_0 + ts_xzzzz_yyz[i] * pa_x[i];

        ts_xxzzzz_yzz[i] = ts_zzzz_yzz[i] * fe_0 + ts_xzzzz_yzz[i] * pa_x[i];

        ts_xxzzzz_zzz[i] = ts_zzzz_zzz[i] * fe_0 + ts_xzzzz_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : IF

    auto ts_xyyyyy_xxx = pbuffer.data(idx_ovl_if + 150);

    auto ts_xyyyyy_xxy = pbuffer.data(idx_ovl_if + 151);

    auto ts_xyyyyy_xxz = pbuffer.data(idx_ovl_if + 152);

    auto ts_xyyyyy_xyy = pbuffer.data(idx_ovl_if + 153);

    auto ts_xyyyyy_xyz = pbuffer.data(idx_ovl_if + 154);

    auto ts_xyyyyy_xzz = pbuffer.data(idx_ovl_if + 155);

    auto ts_xyyyyy_yyy = pbuffer.data(idx_ovl_if + 156);

    auto ts_xyyyyy_yyz = pbuffer.data(idx_ovl_if + 157);

    auto ts_xyyyyy_yzz = pbuffer.data(idx_ovl_if + 158);

    auto ts_xyyyyy_zzz = pbuffer.data(idx_ovl_if + 159);

    #pragma omp simd aligned(pa_x, ts_xyyyyy_xxx, ts_xyyyyy_xxy, ts_xyyyyy_xxz, ts_xyyyyy_xyy, ts_xyyyyy_xyz, ts_xyyyyy_xzz, ts_xyyyyy_yyy, ts_xyyyyy_yyz, ts_xyyyyy_yzz, ts_xyyyyy_zzz, ts_yyyyy_xx, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xz, ts_yyyyy_xzz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyy_zz, ts_yyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyy_xxx[i] = 3.0 * ts_yyyyy_xx[i] * fe_0 + ts_yyyyy_xxx[i] * pa_x[i];

        ts_xyyyyy_xxy[i] = 2.0 * ts_yyyyy_xy[i] * fe_0 + ts_yyyyy_xxy[i] * pa_x[i];

        ts_xyyyyy_xxz[i] = 2.0 * ts_yyyyy_xz[i] * fe_0 + ts_yyyyy_xxz[i] * pa_x[i];

        ts_xyyyyy_xyy[i] = ts_yyyyy_yy[i] * fe_0 + ts_yyyyy_xyy[i] * pa_x[i];

        ts_xyyyyy_xyz[i] = ts_yyyyy_yz[i] * fe_0 + ts_yyyyy_xyz[i] * pa_x[i];

        ts_xyyyyy_xzz[i] = ts_yyyyy_zz[i] * fe_0 + ts_yyyyy_xzz[i] * pa_x[i];

        ts_xyyyyy_yyy[i] = ts_yyyyy_yyy[i] * pa_x[i];

        ts_xyyyyy_yyz[i] = ts_yyyyy_yyz[i] * pa_x[i];

        ts_xyyyyy_yzz[i] = ts_yyyyy_yzz[i] * pa_x[i];

        ts_xyyyyy_zzz[i] = ts_yyyyy_zzz[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : IF

    auto ts_xyyyyz_xxx = pbuffer.data(idx_ovl_if + 160);

    auto ts_xyyyyz_xxy = pbuffer.data(idx_ovl_if + 161);

    auto ts_xyyyyz_xxz = pbuffer.data(idx_ovl_if + 162);

    auto ts_xyyyyz_xyy = pbuffer.data(idx_ovl_if + 163);

    auto ts_xyyyyz_xyz = pbuffer.data(idx_ovl_if + 164);

    auto ts_xyyyyz_xzz = pbuffer.data(idx_ovl_if + 165);

    auto ts_xyyyyz_yyy = pbuffer.data(idx_ovl_if + 166);

    auto ts_xyyyyz_yyz = pbuffer.data(idx_ovl_if + 167);

    auto ts_xyyyyz_yzz = pbuffer.data(idx_ovl_if + 168);

    auto ts_xyyyyz_zzz = pbuffer.data(idx_ovl_if + 169);

    #pragma omp simd aligned(pa_x, pa_z, ts_xyyyy_xxx, ts_xyyyy_xxy, ts_xyyyy_xyy, ts_xyyyyz_xxx, ts_xyyyyz_xxy, ts_xyyyyz_xxz, ts_xyyyyz_xyy, ts_xyyyyz_xyz, ts_xyyyyz_xzz, ts_xyyyyz_yyy, ts_xyyyyz_yyz, ts_xyyyyz_yzz, ts_xyyyyz_zzz, ts_yyyyz_xxz, ts_yyyyz_xyz, ts_yyyyz_xz, ts_yyyyz_xzz, ts_yyyyz_yyy, ts_yyyyz_yyz, ts_yyyyz_yz, ts_yyyyz_yzz, ts_yyyyz_zz, ts_yyyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyyz_xxx[i] = ts_xyyyy_xxx[i] * pa_z[i];

        ts_xyyyyz_xxy[i] = ts_xyyyy_xxy[i] * pa_z[i];

        ts_xyyyyz_xxz[i] = 2.0 * ts_yyyyz_xz[i] * fe_0 + ts_yyyyz_xxz[i] * pa_x[i];

        ts_xyyyyz_xyy[i] = ts_xyyyy_xyy[i] * pa_z[i];

        ts_xyyyyz_xyz[i] = ts_yyyyz_yz[i] * fe_0 + ts_yyyyz_xyz[i] * pa_x[i];

        ts_xyyyyz_xzz[i] = ts_yyyyz_zz[i] * fe_0 + ts_yyyyz_xzz[i] * pa_x[i];

        ts_xyyyyz_yyy[i] = ts_yyyyz_yyy[i] * pa_x[i];

        ts_xyyyyz_yyz[i] = ts_yyyyz_yyz[i] * pa_x[i];

        ts_xyyyyz_yzz[i] = ts_yyyyz_yzz[i] * pa_x[i];

        ts_xyyyyz_zzz[i] = ts_yyyyz_zzz[i] * pa_x[i];
    }

    // Set up 170-180 components of targeted buffer : IF

    auto ts_xyyyzz_xxx = pbuffer.data(idx_ovl_if + 170);

    auto ts_xyyyzz_xxy = pbuffer.data(idx_ovl_if + 171);

    auto ts_xyyyzz_xxz = pbuffer.data(idx_ovl_if + 172);

    auto ts_xyyyzz_xyy = pbuffer.data(idx_ovl_if + 173);

    auto ts_xyyyzz_xyz = pbuffer.data(idx_ovl_if + 174);

    auto ts_xyyyzz_xzz = pbuffer.data(idx_ovl_if + 175);

    auto ts_xyyyzz_yyy = pbuffer.data(idx_ovl_if + 176);

    auto ts_xyyyzz_yyz = pbuffer.data(idx_ovl_if + 177);

    auto ts_xyyyzz_yzz = pbuffer.data(idx_ovl_if + 178);

    auto ts_xyyyzz_zzz = pbuffer.data(idx_ovl_if + 179);

    #pragma omp simd aligned(pa_x, ts_xyyyzz_xxx, ts_xyyyzz_xxy, ts_xyyyzz_xxz, ts_xyyyzz_xyy, ts_xyyyzz_xyz, ts_xyyyzz_xzz, ts_xyyyzz_yyy, ts_xyyyzz_yyz, ts_xyyyzz_yzz, ts_xyyyzz_zzz, ts_yyyzz_xx, ts_yyyzz_xxx, ts_yyyzz_xxy, ts_yyyzz_xxz, ts_yyyzz_xy, ts_yyyzz_xyy, ts_yyyzz_xyz, ts_yyyzz_xz, ts_yyyzz_xzz, ts_yyyzz_yy, ts_yyyzz_yyy, ts_yyyzz_yyz, ts_yyyzz_yz, ts_yyyzz_yzz, ts_yyyzz_zz, ts_yyyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyyzz_xxx[i] = 3.0 * ts_yyyzz_xx[i] * fe_0 + ts_yyyzz_xxx[i] * pa_x[i];

        ts_xyyyzz_xxy[i] = 2.0 * ts_yyyzz_xy[i] * fe_0 + ts_yyyzz_xxy[i] * pa_x[i];

        ts_xyyyzz_xxz[i] = 2.0 * ts_yyyzz_xz[i] * fe_0 + ts_yyyzz_xxz[i] * pa_x[i];

        ts_xyyyzz_xyy[i] = ts_yyyzz_yy[i] * fe_0 + ts_yyyzz_xyy[i] * pa_x[i];

        ts_xyyyzz_xyz[i] = ts_yyyzz_yz[i] * fe_0 + ts_yyyzz_xyz[i] * pa_x[i];

        ts_xyyyzz_xzz[i] = ts_yyyzz_zz[i] * fe_0 + ts_yyyzz_xzz[i] * pa_x[i];

        ts_xyyyzz_yyy[i] = ts_yyyzz_yyy[i] * pa_x[i];

        ts_xyyyzz_yyz[i] = ts_yyyzz_yyz[i] * pa_x[i];

        ts_xyyyzz_yzz[i] = ts_yyyzz_yzz[i] * pa_x[i];

        ts_xyyyzz_zzz[i] = ts_yyyzz_zzz[i] * pa_x[i];
    }

    // Set up 180-190 components of targeted buffer : IF

    auto ts_xyyzzz_xxx = pbuffer.data(idx_ovl_if + 180);

    auto ts_xyyzzz_xxy = pbuffer.data(idx_ovl_if + 181);

    auto ts_xyyzzz_xxz = pbuffer.data(idx_ovl_if + 182);

    auto ts_xyyzzz_xyy = pbuffer.data(idx_ovl_if + 183);

    auto ts_xyyzzz_xyz = pbuffer.data(idx_ovl_if + 184);

    auto ts_xyyzzz_xzz = pbuffer.data(idx_ovl_if + 185);

    auto ts_xyyzzz_yyy = pbuffer.data(idx_ovl_if + 186);

    auto ts_xyyzzz_yyz = pbuffer.data(idx_ovl_if + 187);

    auto ts_xyyzzz_yzz = pbuffer.data(idx_ovl_if + 188);

    auto ts_xyyzzz_zzz = pbuffer.data(idx_ovl_if + 189);

    #pragma omp simd aligned(pa_x, ts_xyyzzz_xxx, ts_xyyzzz_xxy, ts_xyyzzz_xxz, ts_xyyzzz_xyy, ts_xyyzzz_xyz, ts_xyyzzz_xzz, ts_xyyzzz_yyy, ts_xyyzzz_yyz, ts_xyyzzz_yzz, ts_xyyzzz_zzz, ts_yyzzz_xx, ts_yyzzz_xxx, ts_yyzzz_xxy, ts_yyzzz_xxz, ts_yyzzz_xy, ts_yyzzz_xyy, ts_yyzzz_xyz, ts_yyzzz_xz, ts_yyzzz_xzz, ts_yyzzz_yy, ts_yyzzz_yyy, ts_yyzzz_yyz, ts_yyzzz_yz, ts_yyzzz_yzz, ts_yyzzz_zz, ts_yyzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyyzzz_xxx[i] = 3.0 * ts_yyzzz_xx[i] * fe_0 + ts_yyzzz_xxx[i] * pa_x[i];

        ts_xyyzzz_xxy[i] = 2.0 * ts_yyzzz_xy[i] * fe_0 + ts_yyzzz_xxy[i] * pa_x[i];

        ts_xyyzzz_xxz[i] = 2.0 * ts_yyzzz_xz[i] * fe_0 + ts_yyzzz_xxz[i] * pa_x[i];

        ts_xyyzzz_xyy[i] = ts_yyzzz_yy[i] * fe_0 + ts_yyzzz_xyy[i] * pa_x[i];

        ts_xyyzzz_xyz[i] = ts_yyzzz_yz[i] * fe_0 + ts_yyzzz_xyz[i] * pa_x[i];

        ts_xyyzzz_xzz[i] = ts_yyzzz_zz[i] * fe_0 + ts_yyzzz_xzz[i] * pa_x[i];

        ts_xyyzzz_yyy[i] = ts_yyzzz_yyy[i] * pa_x[i];

        ts_xyyzzz_yyz[i] = ts_yyzzz_yyz[i] * pa_x[i];

        ts_xyyzzz_yzz[i] = ts_yyzzz_yzz[i] * pa_x[i];

        ts_xyyzzz_zzz[i] = ts_yyzzz_zzz[i] * pa_x[i];
    }

    // Set up 190-200 components of targeted buffer : IF

    auto ts_xyzzzz_xxx = pbuffer.data(idx_ovl_if + 190);

    auto ts_xyzzzz_xxy = pbuffer.data(idx_ovl_if + 191);

    auto ts_xyzzzz_xxz = pbuffer.data(idx_ovl_if + 192);

    auto ts_xyzzzz_xyy = pbuffer.data(idx_ovl_if + 193);

    auto ts_xyzzzz_xyz = pbuffer.data(idx_ovl_if + 194);

    auto ts_xyzzzz_xzz = pbuffer.data(idx_ovl_if + 195);

    auto ts_xyzzzz_yyy = pbuffer.data(idx_ovl_if + 196);

    auto ts_xyzzzz_yyz = pbuffer.data(idx_ovl_if + 197);

    auto ts_xyzzzz_yzz = pbuffer.data(idx_ovl_if + 198);

    auto ts_xyzzzz_zzz = pbuffer.data(idx_ovl_if + 199);

    #pragma omp simd aligned(pa_x, pa_y, ts_xyzzzz_xxx, ts_xyzzzz_xxy, ts_xyzzzz_xxz, ts_xyzzzz_xyy, ts_xyzzzz_xyz, ts_xyzzzz_xzz, ts_xyzzzz_yyy, ts_xyzzzz_yyz, ts_xyzzzz_yzz, ts_xyzzzz_zzz, ts_xzzzz_xxx, ts_xzzzz_xxz, ts_xzzzz_xzz, ts_yzzzz_xxy, ts_yzzzz_xy, ts_yzzzz_xyy, ts_yzzzz_xyz, ts_yzzzz_yy, ts_yzzzz_yyy, ts_yzzzz_yyz, ts_yzzzz_yz, ts_yzzzz_yzz, ts_yzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xyzzzz_xxx[i] = ts_xzzzz_xxx[i] * pa_y[i];

        ts_xyzzzz_xxy[i] = 2.0 * ts_yzzzz_xy[i] * fe_0 + ts_yzzzz_xxy[i] * pa_x[i];

        ts_xyzzzz_xxz[i] = ts_xzzzz_xxz[i] * pa_y[i];

        ts_xyzzzz_xyy[i] = ts_yzzzz_yy[i] * fe_0 + ts_yzzzz_xyy[i] * pa_x[i];

        ts_xyzzzz_xyz[i] = ts_yzzzz_yz[i] * fe_0 + ts_yzzzz_xyz[i] * pa_x[i];

        ts_xyzzzz_xzz[i] = ts_xzzzz_xzz[i] * pa_y[i];

        ts_xyzzzz_yyy[i] = ts_yzzzz_yyy[i] * pa_x[i];

        ts_xyzzzz_yyz[i] = ts_yzzzz_yyz[i] * pa_x[i];

        ts_xyzzzz_yzz[i] = ts_yzzzz_yzz[i] * pa_x[i];

        ts_xyzzzz_zzz[i] = ts_yzzzz_zzz[i] * pa_x[i];
    }

    // Set up 200-210 components of targeted buffer : IF

    auto ts_xzzzzz_xxx = pbuffer.data(idx_ovl_if + 200);

    auto ts_xzzzzz_xxy = pbuffer.data(idx_ovl_if + 201);

    auto ts_xzzzzz_xxz = pbuffer.data(idx_ovl_if + 202);

    auto ts_xzzzzz_xyy = pbuffer.data(idx_ovl_if + 203);

    auto ts_xzzzzz_xyz = pbuffer.data(idx_ovl_if + 204);

    auto ts_xzzzzz_xzz = pbuffer.data(idx_ovl_if + 205);

    auto ts_xzzzzz_yyy = pbuffer.data(idx_ovl_if + 206);

    auto ts_xzzzzz_yyz = pbuffer.data(idx_ovl_if + 207);

    auto ts_xzzzzz_yzz = pbuffer.data(idx_ovl_if + 208);

    auto ts_xzzzzz_zzz = pbuffer.data(idx_ovl_if + 209);

    #pragma omp simd aligned(pa_x, ts_xzzzzz_xxx, ts_xzzzzz_xxy, ts_xzzzzz_xxz, ts_xzzzzz_xyy, ts_xzzzzz_xyz, ts_xzzzzz_xzz, ts_xzzzzz_yyy, ts_xzzzzz_yyz, ts_xzzzzz_yzz, ts_xzzzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_xzzzzz_xxx[i] = 3.0 * ts_zzzzz_xx[i] * fe_0 + ts_zzzzz_xxx[i] * pa_x[i];

        ts_xzzzzz_xxy[i] = 2.0 * ts_zzzzz_xy[i] * fe_0 + ts_zzzzz_xxy[i] * pa_x[i];

        ts_xzzzzz_xxz[i] = 2.0 * ts_zzzzz_xz[i] * fe_0 + ts_zzzzz_xxz[i] * pa_x[i];

        ts_xzzzzz_xyy[i] = ts_zzzzz_yy[i] * fe_0 + ts_zzzzz_xyy[i] * pa_x[i];

        ts_xzzzzz_xyz[i] = ts_zzzzz_yz[i] * fe_0 + ts_zzzzz_xyz[i] * pa_x[i];

        ts_xzzzzz_xzz[i] = ts_zzzzz_zz[i] * fe_0 + ts_zzzzz_xzz[i] * pa_x[i];

        ts_xzzzzz_yyy[i] = ts_zzzzz_yyy[i] * pa_x[i];

        ts_xzzzzz_yyz[i] = ts_zzzzz_yyz[i] * pa_x[i];

        ts_xzzzzz_yzz[i] = ts_zzzzz_yzz[i] * pa_x[i];

        ts_xzzzzz_zzz[i] = ts_zzzzz_zzz[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : IF

    auto ts_yyyyyy_xxx = pbuffer.data(idx_ovl_if + 210);

    auto ts_yyyyyy_xxy = pbuffer.data(idx_ovl_if + 211);

    auto ts_yyyyyy_xxz = pbuffer.data(idx_ovl_if + 212);

    auto ts_yyyyyy_xyy = pbuffer.data(idx_ovl_if + 213);

    auto ts_yyyyyy_xyz = pbuffer.data(idx_ovl_if + 214);

    auto ts_yyyyyy_xzz = pbuffer.data(idx_ovl_if + 215);

    auto ts_yyyyyy_yyy = pbuffer.data(idx_ovl_if + 216);

    auto ts_yyyyyy_yyz = pbuffer.data(idx_ovl_if + 217);

    auto ts_yyyyyy_yzz = pbuffer.data(idx_ovl_if + 218);

    auto ts_yyyyyy_zzz = pbuffer.data(idx_ovl_if + 219);

    #pragma omp simd aligned(pa_y, ts_yyyy_xxx, ts_yyyy_xxy, ts_yyyy_xxz, ts_yyyy_xyy, ts_yyyy_xyz, ts_yyyy_xzz, ts_yyyy_yyy, ts_yyyy_yyz, ts_yyyy_yzz, ts_yyyy_zzz, ts_yyyyy_xx, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xxz, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_xz, ts_yyyyy_xzz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyy_zz, ts_yyyyy_zzz, ts_yyyyyy_xxx, ts_yyyyyy_xxy, ts_yyyyyy_xxz, ts_yyyyyy_xyy, ts_yyyyyy_xyz, ts_yyyyyy_xzz, ts_yyyyyy_yyy, ts_yyyyyy_yyz, ts_yyyyyy_yzz, ts_yyyyyy_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyy_xxx[i] = 5.0 * ts_yyyy_xxx[i] * fe_0 + ts_yyyyy_xxx[i] * pa_y[i];

        ts_yyyyyy_xxy[i] = 5.0 * ts_yyyy_xxy[i] * fe_0 + ts_yyyyy_xx[i] * fe_0 + ts_yyyyy_xxy[i] * pa_y[i];

        ts_yyyyyy_xxz[i] = 5.0 * ts_yyyy_xxz[i] * fe_0 + ts_yyyyy_xxz[i] * pa_y[i];

        ts_yyyyyy_xyy[i] = 5.0 * ts_yyyy_xyy[i] * fe_0 + 2.0 * ts_yyyyy_xy[i] * fe_0 + ts_yyyyy_xyy[i] * pa_y[i];

        ts_yyyyyy_xyz[i] = 5.0 * ts_yyyy_xyz[i] * fe_0 + ts_yyyyy_xz[i] * fe_0 + ts_yyyyy_xyz[i] * pa_y[i];

        ts_yyyyyy_xzz[i] = 5.0 * ts_yyyy_xzz[i] * fe_0 + ts_yyyyy_xzz[i] * pa_y[i];

        ts_yyyyyy_yyy[i] = 5.0 * ts_yyyy_yyy[i] * fe_0 + 3.0 * ts_yyyyy_yy[i] * fe_0 + ts_yyyyy_yyy[i] * pa_y[i];

        ts_yyyyyy_yyz[i] = 5.0 * ts_yyyy_yyz[i] * fe_0 + 2.0 * ts_yyyyy_yz[i] * fe_0 + ts_yyyyy_yyz[i] * pa_y[i];

        ts_yyyyyy_yzz[i] = 5.0 * ts_yyyy_yzz[i] * fe_0 + ts_yyyyy_zz[i] * fe_0 + ts_yyyyy_yzz[i] * pa_y[i];

        ts_yyyyyy_zzz[i] = 5.0 * ts_yyyy_zzz[i] * fe_0 + ts_yyyyy_zzz[i] * pa_y[i];
    }

    // Set up 220-230 components of targeted buffer : IF

    auto ts_yyyyyz_xxx = pbuffer.data(idx_ovl_if + 220);

    auto ts_yyyyyz_xxy = pbuffer.data(idx_ovl_if + 221);

    auto ts_yyyyyz_xxz = pbuffer.data(idx_ovl_if + 222);

    auto ts_yyyyyz_xyy = pbuffer.data(idx_ovl_if + 223);

    auto ts_yyyyyz_xyz = pbuffer.data(idx_ovl_if + 224);

    auto ts_yyyyyz_xzz = pbuffer.data(idx_ovl_if + 225);

    auto ts_yyyyyz_yyy = pbuffer.data(idx_ovl_if + 226);

    auto ts_yyyyyz_yyz = pbuffer.data(idx_ovl_if + 227);

    auto ts_yyyyyz_yzz = pbuffer.data(idx_ovl_if + 228);

    auto ts_yyyyyz_zzz = pbuffer.data(idx_ovl_if + 229);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyyy_xxx, ts_yyyyy_xxy, ts_yyyyy_xy, ts_yyyyy_xyy, ts_yyyyy_xyz, ts_yyyyy_yy, ts_yyyyy_yyy, ts_yyyyy_yyz, ts_yyyyy_yz, ts_yyyyy_yzz, ts_yyyyyz_xxx, ts_yyyyyz_xxy, ts_yyyyyz_xxz, ts_yyyyyz_xyy, ts_yyyyyz_xyz, ts_yyyyyz_xzz, ts_yyyyyz_yyy, ts_yyyyyz_yyz, ts_yyyyyz_yzz, ts_yyyyyz_zzz, ts_yyyyz_xxz, ts_yyyyz_xzz, ts_yyyyz_zzz, ts_yyyz_xxz, ts_yyyz_xzz, ts_yyyz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyyz_xxx[i] = ts_yyyyy_xxx[i] * pa_z[i];

        ts_yyyyyz_xxy[i] = ts_yyyyy_xxy[i] * pa_z[i];

        ts_yyyyyz_xxz[i] = 4.0 * ts_yyyz_xxz[i] * fe_0 + ts_yyyyz_xxz[i] * pa_y[i];

        ts_yyyyyz_xyy[i] = ts_yyyyy_xyy[i] * pa_z[i];

        ts_yyyyyz_xyz[i] = ts_yyyyy_xy[i] * fe_0 + ts_yyyyy_xyz[i] * pa_z[i];

        ts_yyyyyz_xzz[i] = 4.0 * ts_yyyz_xzz[i] * fe_0 + ts_yyyyz_xzz[i] * pa_y[i];

        ts_yyyyyz_yyy[i] = ts_yyyyy_yyy[i] * pa_z[i];

        ts_yyyyyz_yyz[i] = ts_yyyyy_yy[i] * fe_0 + ts_yyyyy_yyz[i] * pa_z[i];

        ts_yyyyyz_yzz[i] = 2.0 * ts_yyyyy_yz[i] * fe_0 + ts_yyyyy_yzz[i] * pa_z[i];

        ts_yyyyyz_zzz[i] = 4.0 * ts_yyyz_zzz[i] * fe_0 + ts_yyyyz_zzz[i] * pa_y[i];
    }

    // Set up 230-240 components of targeted buffer : IF

    auto ts_yyyyzz_xxx = pbuffer.data(idx_ovl_if + 230);

    auto ts_yyyyzz_xxy = pbuffer.data(idx_ovl_if + 231);

    auto ts_yyyyzz_xxz = pbuffer.data(idx_ovl_if + 232);

    auto ts_yyyyzz_xyy = pbuffer.data(idx_ovl_if + 233);

    auto ts_yyyyzz_xyz = pbuffer.data(idx_ovl_if + 234);

    auto ts_yyyyzz_xzz = pbuffer.data(idx_ovl_if + 235);

    auto ts_yyyyzz_yyy = pbuffer.data(idx_ovl_if + 236);

    auto ts_yyyyzz_yyz = pbuffer.data(idx_ovl_if + 237);

    auto ts_yyyyzz_yzz = pbuffer.data(idx_ovl_if + 238);

    auto ts_yyyyzz_zzz = pbuffer.data(idx_ovl_if + 239);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyy_xxy, ts_yyyy_xyy, ts_yyyy_yyy, ts_yyyyz_xxy, ts_yyyyz_xyy, ts_yyyyz_yyy, ts_yyyyzz_xxx, ts_yyyyzz_xxy, ts_yyyyzz_xxz, ts_yyyyzz_xyy, ts_yyyyzz_xyz, ts_yyyyzz_xzz, ts_yyyyzz_yyy, ts_yyyyzz_yyz, ts_yyyyzz_yzz, ts_yyyyzz_zzz, ts_yyyzz_xxx, ts_yyyzz_xxz, ts_yyyzz_xyz, ts_yyyzz_xz, ts_yyyzz_xzz, ts_yyyzz_yyz, ts_yyyzz_yz, ts_yyyzz_yzz, ts_yyyzz_zz, ts_yyyzz_zzz, ts_yyzz_xxx, ts_yyzz_xxz, ts_yyzz_xyz, ts_yyzz_xzz, ts_yyzz_yyz, ts_yyzz_yzz, ts_yyzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyyzz_xxx[i] = 3.0 * ts_yyzz_xxx[i] * fe_0 + ts_yyyzz_xxx[i] * pa_y[i];

        ts_yyyyzz_xxy[i] = ts_yyyy_xxy[i] * fe_0 + ts_yyyyz_xxy[i] * pa_z[i];

        ts_yyyyzz_xxz[i] = 3.0 * ts_yyzz_xxz[i] * fe_0 + ts_yyyzz_xxz[i] * pa_y[i];

        ts_yyyyzz_xyy[i] = ts_yyyy_xyy[i] * fe_0 + ts_yyyyz_xyy[i] * pa_z[i];

        ts_yyyyzz_xyz[i] = 3.0 * ts_yyzz_xyz[i] * fe_0 + ts_yyyzz_xz[i] * fe_0 + ts_yyyzz_xyz[i] * pa_y[i];

        ts_yyyyzz_xzz[i] = 3.0 * ts_yyzz_xzz[i] * fe_0 + ts_yyyzz_xzz[i] * pa_y[i];

        ts_yyyyzz_yyy[i] = ts_yyyy_yyy[i] * fe_0 + ts_yyyyz_yyy[i] * pa_z[i];

        ts_yyyyzz_yyz[i] = 3.0 * ts_yyzz_yyz[i] * fe_0 + 2.0 * ts_yyyzz_yz[i] * fe_0 + ts_yyyzz_yyz[i] * pa_y[i];

        ts_yyyyzz_yzz[i] = 3.0 * ts_yyzz_yzz[i] * fe_0 + ts_yyyzz_zz[i] * fe_0 + ts_yyyzz_yzz[i] * pa_y[i];

        ts_yyyyzz_zzz[i] = 3.0 * ts_yyzz_zzz[i] * fe_0 + ts_yyyzz_zzz[i] * pa_y[i];
    }

    // Set up 240-250 components of targeted buffer : IF

    auto ts_yyyzzz_xxx = pbuffer.data(idx_ovl_if + 240);

    auto ts_yyyzzz_xxy = pbuffer.data(idx_ovl_if + 241);

    auto ts_yyyzzz_xxz = pbuffer.data(idx_ovl_if + 242);

    auto ts_yyyzzz_xyy = pbuffer.data(idx_ovl_if + 243);

    auto ts_yyyzzz_xyz = pbuffer.data(idx_ovl_if + 244);

    auto ts_yyyzzz_xzz = pbuffer.data(idx_ovl_if + 245);

    auto ts_yyyzzz_yyy = pbuffer.data(idx_ovl_if + 246);

    auto ts_yyyzzz_yyz = pbuffer.data(idx_ovl_if + 247);

    auto ts_yyyzzz_yzz = pbuffer.data(idx_ovl_if + 248);

    auto ts_yyyzzz_zzz = pbuffer.data(idx_ovl_if + 249);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyyz_xxy, ts_yyyz_xyy, ts_yyyz_yyy, ts_yyyzz_xxy, ts_yyyzz_xyy, ts_yyyzz_yyy, ts_yyyzzz_xxx, ts_yyyzzz_xxy, ts_yyyzzz_xxz, ts_yyyzzz_xyy, ts_yyyzzz_xyz, ts_yyyzzz_xzz, ts_yyyzzz_yyy, ts_yyyzzz_yyz, ts_yyyzzz_yzz, ts_yyyzzz_zzz, ts_yyzzz_xxx, ts_yyzzz_xxz, ts_yyzzz_xyz, ts_yyzzz_xz, ts_yyzzz_xzz, ts_yyzzz_yyz, ts_yyzzz_yz, ts_yyzzz_yzz, ts_yyzzz_zz, ts_yyzzz_zzz, ts_yzzz_xxx, ts_yzzz_xxz, ts_yzzz_xyz, ts_yzzz_xzz, ts_yzzz_yyz, ts_yzzz_yzz, ts_yzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyyzzz_xxx[i] = 2.0 * ts_yzzz_xxx[i] * fe_0 + ts_yyzzz_xxx[i] * pa_y[i];

        ts_yyyzzz_xxy[i] = 2.0 * ts_yyyz_xxy[i] * fe_0 + ts_yyyzz_xxy[i] * pa_z[i];

        ts_yyyzzz_xxz[i] = 2.0 * ts_yzzz_xxz[i] * fe_0 + ts_yyzzz_xxz[i] * pa_y[i];

        ts_yyyzzz_xyy[i] = 2.0 * ts_yyyz_xyy[i] * fe_0 + ts_yyyzz_xyy[i] * pa_z[i];

        ts_yyyzzz_xyz[i] = 2.0 * ts_yzzz_xyz[i] * fe_0 + ts_yyzzz_xz[i] * fe_0 + ts_yyzzz_xyz[i] * pa_y[i];

        ts_yyyzzz_xzz[i] = 2.0 * ts_yzzz_xzz[i] * fe_0 + ts_yyzzz_xzz[i] * pa_y[i];

        ts_yyyzzz_yyy[i] = 2.0 * ts_yyyz_yyy[i] * fe_0 + ts_yyyzz_yyy[i] * pa_z[i];

        ts_yyyzzz_yyz[i] = 2.0 * ts_yzzz_yyz[i] * fe_0 + 2.0 * ts_yyzzz_yz[i] * fe_0 + ts_yyzzz_yyz[i] * pa_y[i];

        ts_yyyzzz_yzz[i] = 2.0 * ts_yzzz_yzz[i] * fe_0 + ts_yyzzz_zz[i] * fe_0 + ts_yyzzz_yzz[i] * pa_y[i];

        ts_yyyzzz_zzz[i] = 2.0 * ts_yzzz_zzz[i] * fe_0 + ts_yyzzz_zzz[i] * pa_y[i];
    }

    // Set up 250-260 components of targeted buffer : IF

    auto ts_yyzzzz_xxx = pbuffer.data(idx_ovl_if + 250);

    auto ts_yyzzzz_xxy = pbuffer.data(idx_ovl_if + 251);

    auto ts_yyzzzz_xxz = pbuffer.data(idx_ovl_if + 252);

    auto ts_yyzzzz_xyy = pbuffer.data(idx_ovl_if + 253);

    auto ts_yyzzzz_xyz = pbuffer.data(idx_ovl_if + 254);

    auto ts_yyzzzz_xzz = pbuffer.data(idx_ovl_if + 255);

    auto ts_yyzzzz_yyy = pbuffer.data(idx_ovl_if + 256);

    auto ts_yyzzzz_yyz = pbuffer.data(idx_ovl_if + 257);

    auto ts_yyzzzz_yzz = pbuffer.data(idx_ovl_if + 258);

    auto ts_yyzzzz_zzz = pbuffer.data(idx_ovl_if + 259);

    #pragma omp simd aligned(pa_y, pa_z, ts_yyzz_xxy, ts_yyzz_xyy, ts_yyzz_yyy, ts_yyzzz_xxy, ts_yyzzz_xyy, ts_yyzzz_yyy, ts_yyzzzz_xxx, ts_yyzzzz_xxy, ts_yyzzzz_xxz, ts_yyzzzz_xyy, ts_yyzzzz_xyz, ts_yyzzzz_xzz, ts_yyzzzz_yyy, ts_yyzzzz_yyz, ts_yyzzzz_yzz, ts_yyzzzz_zzz, ts_yzzzz_xxx, ts_yzzzz_xxz, ts_yzzzz_xyz, ts_yzzzz_xz, ts_yzzzz_xzz, ts_yzzzz_yyz, ts_yzzzz_yz, ts_yzzzz_yzz, ts_yzzzz_zz, ts_yzzzz_zzz, ts_zzzz_xxx, ts_zzzz_xxz, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yyzzzz_xxx[i] = ts_zzzz_xxx[i] * fe_0 + ts_yzzzz_xxx[i] * pa_y[i];

        ts_yyzzzz_xxy[i] = 3.0 * ts_yyzz_xxy[i] * fe_0 + ts_yyzzz_xxy[i] * pa_z[i];

        ts_yyzzzz_xxz[i] = ts_zzzz_xxz[i] * fe_0 + ts_yzzzz_xxz[i] * pa_y[i];

        ts_yyzzzz_xyy[i] = 3.0 * ts_yyzz_xyy[i] * fe_0 + ts_yyzzz_xyy[i] * pa_z[i];

        ts_yyzzzz_xyz[i] = ts_zzzz_xyz[i] * fe_0 + ts_yzzzz_xz[i] * fe_0 + ts_yzzzz_xyz[i] * pa_y[i];

        ts_yyzzzz_xzz[i] = ts_zzzz_xzz[i] * fe_0 + ts_yzzzz_xzz[i] * pa_y[i];

        ts_yyzzzz_yyy[i] = 3.0 * ts_yyzz_yyy[i] * fe_0 + ts_yyzzz_yyy[i] * pa_z[i];

        ts_yyzzzz_yyz[i] = ts_zzzz_yyz[i] * fe_0 + 2.0 * ts_yzzzz_yz[i] * fe_0 + ts_yzzzz_yyz[i] * pa_y[i];

        ts_yyzzzz_yzz[i] = ts_zzzz_yzz[i] * fe_0 + ts_yzzzz_zz[i] * fe_0 + ts_yzzzz_yzz[i] * pa_y[i];

        ts_yyzzzz_zzz[i] = ts_zzzz_zzz[i] * fe_0 + ts_yzzzz_zzz[i] * pa_y[i];
    }

    // Set up 260-270 components of targeted buffer : IF

    auto ts_yzzzzz_xxx = pbuffer.data(idx_ovl_if + 260);

    auto ts_yzzzzz_xxy = pbuffer.data(idx_ovl_if + 261);

    auto ts_yzzzzz_xxz = pbuffer.data(idx_ovl_if + 262);

    auto ts_yzzzzz_xyy = pbuffer.data(idx_ovl_if + 263);

    auto ts_yzzzzz_xyz = pbuffer.data(idx_ovl_if + 264);

    auto ts_yzzzzz_xzz = pbuffer.data(idx_ovl_if + 265);

    auto ts_yzzzzz_yyy = pbuffer.data(idx_ovl_if + 266);

    auto ts_yzzzzz_yyz = pbuffer.data(idx_ovl_if + 267);

    auto ts_yzzzzz_yzz = pbuffer.data(idx_ovl_if + 268);

    auto ts_yzzzzz_zzz = pbuffer.data(idx_ovl_if + 269);

    #pragma omp simd aligned(pa_y, ts_yzzzzz_xxx, ts_yzzzzz_xxy, ts_yzzzzz_xxz, ts_yzzzzz_xyy, ts_yzzzzz_xyz, ts_yzzzzz_xzz, ts_yzzzzz_yyy, ts_yzzzzz_yyz, ts_yzzzzz_yzz, ts_yzzzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_yzzzzz_xxx[i] = ts_zzzzz_xxx[i] * pa_y[i];

        ts_yzzzzz_xxy[i] = ts_zzzzz_xx[i] * fe_0 + ts_zzzzz_xxy[i] * pa_y[i];

        ts_yzzzzz_xxz[i] = ts_zzzzz_xxz[i] * pa_y[i];

        ts_yzzzzz_xyy[i] = 2.0 * ts_zzzzz_xy[i] * fe_0 + ts_zzzzz_xyy[i] * pa_y[i];

        ts_yzzzzz_xyz[i] = ts_zzzzz_xz[i] * fe_0 + ts_zzzzz_xyz[i] * pa_y[i];

        ts_yzzzzz_xzz[i] = ts_zzzzz_xzz[i] * pa_y[i];

        ts_yzzzzz_yyy[i] = 3.0 * ts_zzzzz_yy[i] * fe_0 + ts_zzzzz_yyy[i] * pa_y[i];

        ts_yzzzzz_yyz[i] = 2.0 * ts_zzzzz_yz[i] * fe_0 + ts_zzzzz_yyz[i] * pa_y[i];

        ts_yzzzzz_yzz[i] = ts_zzzzz_zz[i] * fe_0 + ts_zzzzz_yzz[i] * pa_y[i];

        ts_yzzzzz_zzz[i] = ts_zzzzz_zzz[i] * pa_y[i];
    }

    // Set up 270-280 components of targeted buffer : IF

    auto ts_zzzzzz_xxx = pbuffer.data(idx_ovl_if + 270);

    auto ts_zzzzzz_xxy = pbuffer.data(idx_ovl_if + 271);

    auto ts_zzzzzz_xxz = pbuffer.data(idx_ovl_if + 272);

    auto ts_zzzzzz_xyy = pbuffer.data(idx_ovl_if + 273);

    auto ts_zzzzzz_xyz = pbuffer.data(idx_ovl_if + 274);

    auto ts_zzzzzz_xzz = pbuffer.data(idx_ovl_if + 275);

    auto ts_zzzzzz_yyy = pbuffer.data(idx_ovl_if + 276);

    auto ts_zzzzzz_yyz = pbuffer.data(idx_ovl_if + 277);

    auto ts_zzzzzz_yzz = pbuffer.data(idx_ovl_if + 278);

    auto ts_zzzzzz_zzz = pbuffer.data(idx_ovl_if + 279);

    #pragma omp simd aligned(pa_z, ts_zzzz_xxx, ts_zzzz_xxy, ts_zzzz_xxz, ts_zzzz_xyy, ts_zzzz_xyz, ts_zzzz_xzz, ts_zzzz_yyy, ts_zzzz_yyz, ts_zzzz_yzz, ts_zzzz_zzz, ts_zzzzz_xx, ts_zzzzz_xxx, ts_zzzzz_xxy, ts_zzzzz_xxz, ts_zzzzz_xy, ts_zzzzz_xyy, ts_zzzzz_xyz, ts_zzzzz_xz, ts_zzzzz_xzz, ts_zzzzz_yy, ts_zzzzz_yyy, ts_zzzzz_yyz, ts_zzzzz_yz, ts_zzzzz_yzz, ts_zzzzz_zz, ts_zzzzz_zzz, ts_zzzzzz_xxx, ts_zzzzzz_xxy, ts_zzzzzz_xxz, ts_zzzzzz_xyy, ts_zzzzzz_xyz, ts_zzzzzz_xzz, ts_zzzzzz_yyy, ts_zzzzzz_yyz, ts_zzzzzz_yzz, ts_zzzzzz_zzz, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ts_zzzzzz_xxx[i] = 5.0 * ts_zzzz_xxx[i] * fe_0 + ts_zzzzz_xxx[i] * pa_z[i];

        ts_zzzzzz_xxy[i] = 5.0 * ts_zzzz_xxy[i] * fe_0 + ts_zzzzz_xxy[i] * pa_z[i];

        ts_zzzzzz_xxz[i] = 5.0 * ts_zzzz_xxz[i] * fe_0 + ts_zzzzz_xx[i] * fe_0 + ts_zzzzz_xxz[i] * pa_z[i];

        ts_zzzzzz_xyy[i] = 5.0 * ts_zzzz_xyy[i] * fe_0 + ts_zzzzz_xyy[i] * pa_z[i];

        ts_zzzzzz_xyz[i] = 5.0 * ts_zzzz_xyz[i] * fe_0 + ts_zzzzz_xy[i] * fe_0 + ts_zzzzz_xyz[i] * pa_z[i];

        ts_zzzzzz_xzz[i] = 5.0 * ts_zzzz_xzz[i] * fe_0 + 2.0 * ts_zzzzz_xz[i] * fe_0 + ts_zzzzz_xzz[i] * pa_z[i];

        ts_zzzzzz_yyy[i] = 5.0 * ts_zzzz_yyy[i] * fe_0 + ts_zzzzz_yyy[i] * pa_z[i];

        ts_zzzzzz_yyz[i] = 5.0 * ts_zzzz_yyz[i] * fe_0 + ts_zzzzz_yy[i] * fe_0 + ts_zzzzz_yyz[i] * pa_z[i];

        ts_zzzzzz_yzz[i] = 5.0 * ts_zzzz_yzz[i] * fe_0 + 2.0 * ts_zzzzz_yz[i] * fe_0 + ts_zzzzz_yzz[i] * pa_z[i];

        ts_zzzzzz_zzz[i] = 5.0 * ts_zzzz_zzz[i] * fe_0 + 3.0 * ts_zzzzz_zz[i] * fe_0 + ts_zzzzz_zzz[i] * pa_z[i];
    }

}

} // ovlrec namespace

