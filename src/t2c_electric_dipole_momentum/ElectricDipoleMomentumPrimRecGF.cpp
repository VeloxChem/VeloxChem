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

#include "ElectricDipoleMomentumPrimRecGF.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_gf(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_gf,
                                      const size_t              idx_dip_df,
                                      const size_t              idx_dip_fd,
                                      const size_t              idx_ovl_ff,
                                      const size_t              idx_dip_ff,
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

    // Set up components of auxiliary buffer : DF

    auto tr_x_xx_xxx = pbuffer.data(idx_dip_df);

    auto tr_x_xx_xxy = pbuffer.data(idx_dip_df + 1);

    auto tr_x_xx_xxz = pbuffer.data(idx_dip_df + 2);

    auto tr_x_xx_xyy = pbuffer.data(idx_dip_df + 3);

    auto tr_x_xx_xyz = pbuffer.data(idx_dip_df + 4);

    auto tr_x_xx_xzz = pbuffer.data(idx_dip_df + 5);

    auto tr_x_xx_yyy = pbuffer.data(idx_dip_df + 6);

    auto tr_x_xx_yyz = pbuffer.data(idx_dip_df + 7);

    auto tr_x_xx_yzz = pbuffer.data(idx_dip_df + 8);

    auto tr_x_xx_zzz = pbuffer.data(idx_dip_df + 9);

    auto tr_x_xy_xxx = pbuffer.data(idx_dip_df + 10);

    auto tr_x_xy_xxz = pbuffer.data(idx_dip_df + 12);

    auto tr_x_xy_xzz = pbuffer.data(idx_dip_df + 15);

    auto tr_x_xz_xxx = pbuffer.data(idx_dip_df + 20);

    auto tr_x_xz_xxy = pbuffer.data(idx_dip_df + 21);

    auto tr_x_xz_xxz = pbuffer.data(idx_dip_df + 22);

    auto tr_x_xz_xyy = pbuffer.data(idx_dip_df + 23);

    auto tr_x_xz_xzz = pbuffer.data(idx_dip_df + 25);

    auto tr_x_yy_xxx = pbuffer.data(idx_dip_df + 30);

    auto tr_x_yy_xxy = pbuffer.data(idx_dip_df + 31);

    auto tr_x_yy_xxz = pbuffer.data(idx_dip_df + 32);

    auto tr_x_yy_xyy = pbuffer.data(idx_dip_df + 33);

    auto tr_x_yy_xyz = pbuffer.data(idx_dip_df + 34);

    auto tr_x_yy_xzz = pbuffer.data(idx_dip_df + 35);

    auto tr_x_yy_yyy = pbuffer.data(idx_dip_df + 36);

    auto tr_x_yy_yyz = pbuffer.data(idx_dip_df + 37);

    auto tr_x_yy_yzz = pbuffer.data(idx_dip_df + 38);

    auto tr_x_yy_zzz = pbuffer.data(idx_dip_df + 39);

    auto tr_x_yz_xxz = pbuffer.data(idx_dip_df + 42);

    auto tr_x_yz_xzz = pbuffer.data(idx_dip_df + 45);

    auto tr_x_yz_zzz = pbuffer.data(idx_dip_df + 49);

    auto tr_x_zz_xxx = pbuffer.data(idx_dip_df + 50);

    auto tr_x_zz_xxy = pbuffer.data(idx_dip_df + 51);

    auto tr_x_zz_xxz = pbuffer.data(idx_dip_df + 52);

    auto tr_x_zz_xyy = pbuffer.data(idx_dip_df + 53);

    auto tr_x_zz_xyz = pbuffer.data(idx_dip_df + 54);

    auto tr_x_zz_xzz = pbuffer.data(idx_dip_df + 55);

    auto tr_x_zz_yyy = pbuffer.data(idx_dip_df + 56);

    auto tr_x_zz_yyz = pbuffer.data(idx_dip_df + 57);

    auto tr_x_zz_yzz = pbuffer.data(idx_dip_df + 58);

    auto tr_x_zz_zzz = pbuffer.data(idx_dip_df + 59);

    auto tr_y_xx_xxx = pbuffer.data(idx_dip_df + 60);

    auto tr_y_xx_xxy = pbuffer.data(idx_dip_df + 61);

    auto tr_y_xx_xxz = pbuffer.data(idx_dip_df + 62);

    auto tr_y_xx_xyy = pbuffer.data(idx_dip_df + 63);

    auto tr_y_xx_xyz = pbuffer.data(idx_dip_df + 64);

    auto tr_y_xx_xzz = pbuffer.data(idx_dip_df + 65);

    auto tr_y_xx_yyy = pbuffer.data(idx_dip_df + 66);

    auto tr_y_xx_yyz = pbuffer.data(idx_dip_df + 67);

    auto tr_y_xx_yzz = pbuffer.data(idx_dip_df + 68);

    auto tr_y_xx_zzz = pbuffer.data(idx_dip_df + 69);

    auto tr_y_xy_xxy = pbuffer.data(idx_dip_df + 71);

    auto tr_y_xy_xyy = pbuffer.data(idx_dip_df + 73);

    auto tr_y_xy_xyz = pbuffer.data(idx_dip_df + 74);

    auto tr_y_xy_yyy = pbuffer.data(idx_dip_df + 76);

    auto tr_y_xy_yyz = pbuffer.data(idx_dip_df + 77);

    auto tr_y_xy_yzz = pbuffer.data(idx_dip_df + 78);

    auto tr_y_xy_zzz = pbuffer.data(idx_dip_df + 79);

    auto tr_y_xz_yyz = pbuffer.data(idx_dip_df + 87);

    auto tr_y_xz_yzz = pbuffer.data(idx_dip_df + 88);

    auto tr_y_xz_zzz = pbuffer.data(idx_dip_df + 89);

    auto tr_y_yy_xxx = pbuffer.data(idx_dip_df + 90);

    auto tr_y_yy_xxy = pbuffer.data(idx_dip_df + 91);

    auto tr_y_yy_xxz = pbuffer.data(idx_dip_df + 92);

    auto tr_y_yy_xyy = pbuffer.data(idx_dip_df + 93);

    auto tr_y_yy_xyz = pbuffer.data(idx_dip_df + 94);

    auto tr_y_yy_xzz = pbuffer.data(idx_dip_df + 95);

    auto tr_y_yy_yyy = pbuffer.data(idx_dip_df + 96);

    auto tr_y_yy_yyz = pbuffer.data(idx_dip_df + 97);

    auto tr_y_yy_yzz = pbuffer.data(idx_dip_df + 98);

    auto tr_y_yy_zzz = pbuffer.data(idx_dip_df + 99);

    auto tr_y_yz_xxy = pbuffer.data(idx_dip_df + 101);

    auto tr_y_yz_xyy = pbuffer.data(idx_dip_df + 103);

    auto tr_y_yz_yyy = pbuffer.data(idx_dip_df + 106);

    auto tr_y_yz_yyz = pbuffer.data(idx_dip_df + 107);

    auto tr_y_yz_yzz = pbuffer.data(idx_dip_df + 108);

    auto tr_y_yz_zzz = pbuffer.data(idx_dip_df + 109);

    auto tr_y_zz_xxx = pbuffer.data(idx_dip_df + 110);

    auto tr_y_zz_xxy = pbuffer.data(idx_dip_df + 111);

    auto tr_y_zz_xxz = pbuffer.data(idx_dip_df + 112);

    auto tr_y_zz_xyy = pbuffer.data(idx_dip_df + 113);

    auto tr_y_zz_xyz = pbuffer.data(idx_dip_df + 114);

    auto tr_y_zz_xzz = pbuffer.data(idx_dip_df + 115);

    auto tr_y_zz_yyy = pbuffer.data(idx_dip_df + 116);

    auto tr_y_zz_yyz = pbuffer.data(idx_dip_df + 117);

    auto tr_y_zz_yzz = pbuffer.data(idx_dip_df + 118);

    auto tr_y_zz_zzz = pbuffer.data(idx_dip_df + 119);

    auto tr_z_xx_xxx = pbuffer.data(idx_dip_df + 120);

    auto tr_z_xx_xxy = pbuffer.data(idx_dip_df + 121);

    auto tr_z_xx_xxz = pbuffer.data(idx_dip_df + 122);

    auto tr_z_xx_xyy = pbuffer.data(idx_dip_df + 123);

    auto tr_z_xx_xyz = pbuffer.data(idx_dip_df + 124);

    auto tr_z_xx_xzz = pbuffer.data(idx_dip_df + 125);

    auto tr_z_xx_yyy = pbuffer.data(idx_dip_df + 126);

    auto tr_z_xx_yyz = pbuffer.data(idx_dip_df + 127);

    auto tr_z_xx_yzz = pbuffer.data(idx_dip_df + 128);

    auto tr_z_xx_zzz = pbuffer.data(idx_dip_df + 129);

    auto tr_z_xy_yyy = pbuffer.data(idx_dip_df + 136);

    auto tr_z_xy_yyz = pbuffer.data(idx_dip_df + 137);

    auto tr_z_xy_yzz = pbuffer.data(idx_dip_df + 138);

    auto tr_z_xz_xxz = pbuffer.data(idx_dip_df + 142);

    auto tr_z_xz_xyz = pbuffer.data(idx_dip_df + 144);

    auto tr_z_xz_xzz = pbuffer.data(idx_dip_df + 145);

    auto tr_z_xz_yyy = pbuffer.data(idx_dip_df + 146);

    auto tr_z_xz_yyz = pbuffer.data(idx_dip_df + 147);

    auto tr_z_xz_yzz = pbuffer.data(idx_dip_df + 148);

    auto tr_z_xz_zzz = pbuffer.data(idx_dip_df + 149);

    auto tr_z_yy_xxx = pbuffer.data(idx_dip_df + 150);

    auto tr_z_yy_xxy = pbuffer.data(idx_dip_df + 151);

    auto tr_z_yy_xxz = pbuffer.data(idx_dip_df + 152);

    auto tr_z_yy_xyy = pbuffer.data(idx_dip_df + 153);

    auto tr_z_yy_xyz = pbuffer.data(idx_dip_df + 154);

    auto tr_z_yy_xzz = pbuffer.data(idx_dip_df + 155);

    auto tr_z_yy_yyy = pbuffer.data(idx_dip_df + 156);

    auto tr_z_yy_yyz = pbuffer.data(idx_dip_df + 157);

    auto tr_z_yy_yzz = pbuffer.data(idx_dip_df + 158);

    auto tr_z_yy_zzz = pbuffer.data(idx_dip_df + 159);

    auto tr_z_yz_xxx = pbuffer.data(idx_dip_df + 160);

    auto tr_z_yz_xxz = pbuffer.data(idx_dip_df + 162);

    auto tr_z_yz_xyz = pbuffer.data(idx_dip_df + 164);

    auto tr_z_yz_xzz = pbuffer.data(idx_dip_df + 165);

    auto tr_z_yz_yyy = pbuffer.data(idx_dip_df + 166);

    auto tr_z_yz_yyz = pbuffer.data(idx_dip_df + 167);

    auto tr_z_yz_yzz = pbuffer.data(idx_dip_df + 168);

    auto tr_z_yz_zzz = pbuffer.data(idx_dip_df + 169);

    auto tr_z_zz_xxx = pbuffer.data(idx_dip_df + 170);

    auto tr_z_zz_xxy = pbuffer.data(idx_dip_df + 171);

    auto tr_z_zz_xxz = pbuffer.data(idx_dip_df + 172);

    auto tr_z_zz_xyy = pbuffer.data(idx_dip_df + 173);

    auto tr_z_zz_xyz = pbuffer.data(idx_dip_df + 174);

    auto tr_z_zz_xzz = pbuffer.data(idx_dip_df + 175);

    auto tr_z_zz_yyy = pbuffer.data(idx_dip_df + 176);

    auto tr_z_zz_yyz = pbuffer.data(idx_dip_df + 177);

    auto tr_z_zz_yzz = pbuffer.data(idx_dip_df + 178);

    auto tr_z_zz_zzz = pbuffer.data(idx_dip_df + 179);

    // Set up components of auxiliary buffer : FD

    auto tr_x_xxx_xx = pbuffer.data(idx_dip_fd);

    auto tr_x_xxx_xy = pbuffer.data(idx_dip_fd + 1);

    auto tr_x_xxx_xz = pbuffer.data(idx_dip_fd + 2);

    auto tr_x_xxx_yy = pbuffer.data(idx_dip_fd + 3);

    auto tr_x_xxx_yz = pbuffer.data(idx_dip_fd + 4);

    auto tr_x_xxx_zz = pbuffer.data(idx_dip_fd + 5);

    auto tr_x_xxy_xx = pbuffer.data(idx_dip_fd + 6);

    auto tr_x_xxy_xy = pbuffer.data(idx_dip_fd + 7);

    auto tr_x_xxy_xz = pbuffer.data(idx_dip_fd + 8);

    auto tr_x_xxz_xx = pbuffer.data(idx_dip_fd + 12);

    auto tr_x_xxz_xy = pbuffer.data(idx_dip_fd + 13);

    auto tr_x_xxz_xz = pbuffer.data(idx_dip_fd + 14);

    auto tr_x_xxz_yz = pbuffer.data(idx_dip_fd + 16);

    auto tr_x_xxz_zz = pbuffer.data(idx_dip_fd + 17);

    auto tr_x_xyy_xy = pbuffer.data(idx_dip_fd + 19);

    auto tr_x_xzz_xx = pbuffer.data(idx_dip_fd + 30);

    auto tr_x_xzz_xy = pbuffer.data(idx_dip_fd + 31);

    auto tr_x_xzz_xz = pbuffer.data(idx_dip_fd + 32);

    auto tr_x_yyy_xx = pbuffer.data(idx_dip_fd + 36);

    auto tr_x_yyy_xy = pbuffer.data(idx_dip_fd + 37);

    auto tr_x_yyy_xz = pbuffer.data(idx_dip_fd + 38);

    auto tr_x_yyy_yy = pbuffer.data(idx_dip_fd + 39);

    auto tr_x_yyy_yz = pbuffer.data(idx_dip_fd + 40);

    auto tr_x_yyy_zz = pbuffer.data(idx_dip_fd + 41);

    auto tr_x_yzz_xz = pbuffer.data(idx_dip_fd + 50);

    auto tr_x_yzz_yz = pbuffer.data(idx_dip_fd + 52);

    auto tr_x_yzz_zz = pbuffer.data(idx_dip_fd + 53);

    auto tr_x_zzz_xx = pbuffer.data(idx_dip_fd + 54);

    auto tr_x_zzz_xy = pbuffer.data(idx_dip_fd + 55);

    auto tr_x_zzz_xz = pbuffer.data(idx_dip_fd + 56);

    auto tr_x_zzz_yy = pbuffer.data(idx_dip_fd + 57);

    auto tr_x_zzz_yz = pbuffer.data(idx_dip_fd + 58);

    auto tr_x_zzz_zz = pbuffer.data(idx_dip_fd + 59);

    auto tr_y_xxx_xx = pbuffer.data(idx_dip_fd + 60);

    auto tr_y_xxx_xy = pbuffer.data(idx_dip_fd + 61);

    auto tr_y_xxx_xz = pbuffer.data(idx_dip_fd + 62);

    auto tr_y_xxx_yy = pbuffer.data(idx_dip_fd + 63);

    auto tr_y_xxx_yz = pbuffer.data(idx_dip_fd + 64);

    auto tr_y_xxx_zz = pbuffer.data(idx_dip_fd + 65);

    auto tr_y_xxy_xy = pbuffer.data(idx_dip_fd + 67);

    auto tr_y_xxy_yy = pbuffer.data(idx_dip_fd + 69);

    auto tr_y_xxy_yz = pbuffer.data(idx_dip_fd + 70);

    auto tr_y_xyy_xx = pbuffer.data(idx_dip_fd + 78);

    auto tr_y_xyy_xy = pbuffer.data(idx_dip_fd + 79);

    auto tr_y_xyy_xz = pbuffer.data(idx_dip_fd + 80);

    auto tr_y_xyy_yy = pbuffer.data(idx_dip_fd + 81);

    auto tr_y_xyy_yz = pbuffer.data(idx_dip_fd + 82);

    auto tr_y_xyy_zz = pbuffer.data(idx_dip_fd + 83);

    auto tr_y_xzz_xz = pbuffer.data(idx_dip_fd + 92);

    auto tr_y_xzz_yz = pbuffer.data(idx_dip_fd + 94);

    auto tr_y_xzz_zz = pbuffer.data(idx_dip_fd + 95);

    auto tr_y_yyy_xx = pbuffer.data(idx_dip_fd + 96);

    auto tr_y_yyy_xy = pbuffer.data(idx_dip_fd + 97);

    auto tr_y_yyy_xz = pbuffer.data(idx_dip_fd + 98);

    auto tr_y_yyy_yy = pbuffer.data(idx_dip_fd + 99);

    auto tr_y_yyy_yz = pbuffer.data(idx_dip_fd + 100);

    auto tr_y_yyy_zz = pbuffer.data(idx_dip_fd + 101);

    auto tr_y_yyz_xy = pbuffer.data(idx_dip_fd + 103);

    auto tr_y_yyz_xz = pbuffer.data(idx_dip_fd + 104);

    auto tr_y_yyz_yy = pbuffer.data(idx_dip_fd + 105);

    auto tr_y_yyz_yz = pbuffer.data(idx_dip_fd + 106);

    auto tr_y_yyz_zz = pbuffer.data(idx_dip_fd + 107);

    auto tr_y_yzz_xx = pbuffer.data(idx_dip_fd + 108);

    auto tr_y_yzz_xy = pbuffer.data(idx_dip_fd + 109);

    auto tr_y_yzz_xz = pbuffer.data(idx_dip_fd + 110);

    auto tr_y_yzz_yy = pbuffer.data(idx_dip_fd + 111);

    auto tr_y_yzz_yz = pbuffer.data(idx_dip_fd + 112);

    auto tr_y_yzz_zz = pbuffer.data(idx_dip_fd + 113);

    auto tr_y_zzz_xx = pbuffer.data(idx_dip_fd + 114);

    auto tr_y_zzz_xy = pbuffer.data(idx_dip_fd + 115);

    auto tr_y_zzz_xz = pbuffer.data(idx_dip_fd + 116);

    auto tr_y_zzz_yy = pbuffer.data(idx_dip_fd + 117);

    auto tr_y_zzz_yz = pbuffer.data(idx_dip_fd + 118);

    auto tr_y_zzz_zz = pbuffer.data(idx_dip_fd + 119);

    auto tr_z_xxx_xx = pbuffer.data(idx_dip_fd + 120);

    auto tr_z_xxx_xy = pbuffer.data(idx_dip_fd + 121);

    auto tr_z_xxx_xz = pbuffer.data(idx_dip_fd + 122);

    auto tr_z_xxx_yy = pbuffer.data(idx_dip_fd + 123);

    auto tr_z_xxx_yz = pbuffer.data(idx_dip_fd + 124);

    auto tr_z_xxx_zz = pbuffer.data(idx_dip_fd + 125);

    auto tr_z_xxz_xx = pbuffer.data(idx_dip_fd + 132);

    auto tr_z_xxz_xy = pbuffer.data(idx_dip_fd + 133);

    auto tr_z_xxz_xz = pbuffer.data(idx_dip_fd + 134);

    auto tr_z_xxz_yz = pbuffer.data(idx_dip_fd + 136);

    auto tr_z_xxz_zz = pbuffer.data(idx_dip_fd + 137);

    auto tr_z_xyy_xy = pbuffer.data(idx_dip_fd + 139);

    auto tr_z_xyy_yy = pbuffer.data(idx_dip_fd + 141);

    auto tr_z_xyy_yz = pbuffer.data(idx_dip_fd + 142);

    auto tr_z_xzz_xx = pbuffer.data(idx_dip_fd + 150);

    auto tr_z_xzz_xy = pbuffer.data(idx_dip_fd + 151);

    auto tr_z_xzz_xz = pbuffer.data(idx_dip_fd + 152);

    auto tr_z_xzz_yy = pbuffer.data(idx_dip_fd + 153);

    auto tr_z_xzz_yz = pbuffer.data(idx_dip_fd + 154);

    auto tr_z_xzz_zz = pbuffer.data(idx_dip_fd + 155);

    auto tr_z_yyy_xx = pbuffer.data(idx_dip_fd + 156);

    auto tr_z_yyy_xy = pbuffer.data(idx_dip_fd + 157);

    auto tr_z_yyy_xz = pbuffer.data(idx_dip_fd + 158);

    auto tr_z_yyy_yy = pbuffer.data(idx_dip_fd + 159);

    auto tr_z_yyy_yz = pbuffer.data(idx_dip_fd + 160);

    auto tr_z_yyy_zz = pbuffer.data(idx_dip_fd + 161);

    auto tr_z_yyz_xx = pbuffer.data(idx_dip_fd + 162);

    auto tr_z_yyz_xy = pbuffer.data(idx_dip_fd + 163);

    auto tr_z_yyz_xz = pbuffer.data(idx_dip_fd + 164);

    auto tr_z_yyz_yy = pbuffer.data(idx_dip_fd + 165);

    auto tr_z_yyz_yz = pbuffer.data(idx_dip_fd + 166);

    auto tr_z_yyz_zz = pbuffer.data(idx_dip_fd + 167);

    auto tr_z_yzz_xx = pbuffer.data(idx_dip_fd + 168);

    auto tr_z_yzz_xy = pbuffer.data(idx_dip_fd + 169);

    auto tr_z_yzz_xz = pbuffer.data(idx_dip_fd + 170);

    auto tr_z_yzz_yy = pbuffer.data(idx_dip_fd + 171);

    auto tr_z_yzz_yz = pbuffer.data(idx_dip_fd + 172);

    auto tr_z_yzz_zz = pbuffer.data(idx_dip_fd + 173);

    auto tr_z_zzz_xx = pbuffer.data(idx_dip_fd + 174);

    auto tr_z_zzz_xy = pbuffer.data(idx_dip_fd + 175);

    auto tr_z_zzz_xz = pbuffer.data(idx_dip_fd + 176);

    auto tr_z_zzz_yy = pbuffer.data(idx_dip_fd + 177);

    auto tr_z_zzz_yz = pbuffer.data(idx_dip_fd + 178);

    auto tr_z_zzz_zz = pbuffer.data(idx_dip_fd + 179);

    // Set up components of auxiliary buffer : FF

    auto ts_xxx_xxx = pbuffer.data(idx_ovl_ff);

    auto ts_xxx_xxy = pbuffer.data(idx_ovl_ff + 1);

    auto ts_xxx_xxz = pbuffer.data(idx_ovl_ff + 2);

    auto ts_xxx_xyy = pbuffer.data(idx_ovl_ff + 3);

    auto ts_xxx_xyz = pbuffer.data(idx_ovl_ff + 4);

    auto ts_xxx_xzz = pbuffer.data(idx_ovl_ff + 5);

    auto ts_xxx_yyy = pbuffer.data(idx_ovl_ff + 6);

    auto ts_xxx_yyz = pbuffer.data(idx_ovl_ff + 7);

    auto ts_xxx_yzz = pbuffer.data(idx_ovl_ff + 8);

    auto ts_xxx_zzz = pbuffer.data(idx_ovl_ff + 9);

    auto ts_xxz_xxz = pbuffer.data(idx_ovl_ff + 22);

    auto ts_xxz_xzz = pbuffer.data(idx_ovl_ff + 25);

    auto ts_xyy_yyy = pbuffer.data(idx_ovl_ff + 36);

    auto ts_xyy_yyz = pbuffer.data(idx_ovl_ff + 37);

    auto ts_xyy_yzz = pbuffer.data(idx_ovl_ff + 38);

    auto ts_xzz_yyz = pbuffer.data(idx_ovl_ff + 57);

    auto ts_xzz_yzz = pbuffer.data(idx_ovl_ff + 58);

    auto ts_xzz_zzz = pbuffer.data(idx_ovl_ff + 59);

    auto ts_yyy_xxx = pbuffer.data(idx_ovl_ff + 60);

    auto ts_yyy_xxy = pbuffer.data(idx_ovl_ff + 61);

    auto ts_yyy_xxz = pbuffer.data(idx_ovl_ff + 62);

    auto ts_yyy_xyy = pbuffer.data(idx_ovl_ff + 63);

    auto ts_yyy_xyz = pbuffer.data(idx_ovl_ff + 64);

    auto ts_yyy_xzz = pbuffer.data(idx_ovl_ff + 65);

    auto ts_yyy_yyy = pbuffer.data(idx_ovl_ff + 66);

    auto ts_yyy_yyz = pbuffer.data(idx_ovl_ff + 67);

    auto ts_yyy_yzz = pbuffer.data(idx_ovl_ff + 68);

    auto ts_yyy_zzz = pbuffer.data(idx_ovl_ff + 69);

    auto ts_yyz_yyz = pbuffer.data(idx_ovl_ff + 77);

    auto ts_yyz_yzz = pbuffer.data(idx_ovl_ff + 78);

    auto ts_yyz_zzz = pbuffer.data(idx_ovl_ff + 79);

    auto ts_yzz_xxz = pbuffer.data(idx_ovl_ff + 82);

    auto ts_yzz_xzz = pbuffer.data(idx_ovl_ff + 85);

    auto ts_yzz_yyy = pbuffer.data(idx_ovl_ff + 86);

    auto ts_yzz_yyz = pbuffer.data(idx_ovl_ff + 87);

    auto ts_yzz_yzz = pbuffer.data(idx_ovl_ff + 88);

    auto ts_yzz_zzz = pbuffer.data(idx_ovl_ff + 89);

    auto ts_zzz_xxx = pbuffer.data(idx_ovl_ff + 90);

    auto ts_zzz_xxy = pbuffer.data(idx_ovl_ff + 91);

    auto ts_zzz_xxz = pbuffer.data(idx_ovl_ff + 92);

    auto ts_zzz_xyy = pbuffer.data(idx_ovl_ff + 93);

    auto ts_zzz_xyz = pbuffer.data(idx_ovl_ff + 94);

    auto ts_zzz_xzz = pbuffer.data(idx_ovl_ff + 95);

    auto ts_zzz_yyy = pbuffer.data(idx_ovl_ff + 96);

    auto ts_zzz_yyz = pbuffer.data(idx_ovl_ff + 97);

    auto ts_zzz_yzz = pbuffer.data(idx_ovl_ff + 98);

    auto ts_zzz_zzz = pbuffer.data(idx_ovl_ff + 99);

    // Set up components of auxiliary buffer : FF

    auto tr_x_xxx_xxx = pbuffer.data(idx_dip_ff);

    auto tr_x_xxx_xxy = pbuffer.data(idx_dip_ff + 1);

    auto tr_x_xxx_xxz = pbuffer.data(idx_dip_ff + 2);

    auto tr_x_xxx_xyy = pbuffer.data(idx_dip_ff + 3);

    auto tr_x_xxx_xyz = pbuffer.data(idx_dip_ff + 4);

    auto tr_x_xxx_xzz = pbuffer.data(idx_dip_ff + 5);

    auto tr_x_xxx_yyy = pbuffer.data(idx_dip_ff + 6);

    auto tr_x_xxx_yyz = pbuffer.data(idx_dip_ff + 7);

    auto tr_x_xxx_yzz = pbuffer.data(idx_dip_ff + 8);

    auto tr_x_xxx_zzz = pbuffer.data(idx_dip_ff + 9);

    auto tr_x_xxy_xxx = pbuffer.data(idx_dip_ff + 10);

    auto tr_x_xxy_xxy = pbuffer.data(idx_dip_ff + 11);

    auto tr_x_xxy_xxz = pbuffer.data(idx_dip_ff + 12);

    auto tr_x_xxy_xyy = pbuffer.data(idx_dip_ff + 13);

    auto tr_x_xxy_xyz = pbuffer.data(idx_dip_ff + 14);

    auto tr_x_xxy_xzz = pbuffer.data(idx_dip_ff + 15);

    auto tr_x_xxy_yyy = pbuffer.data(idx_dip_ff + 16);

    auto tr_x_xxy_zzz = pbuffer.data(idx_dip_ff + 19);

    auto tr_x_xxz_xxx = pbuffer.data(idx_dip_ff + 20);

    auto tr_x_xxz_xxy = pbuffer.data(idx_dip_ff + 21);

    auto tr_x_xxz_xxz = pbuffer.data(idx_dip_ff + 22);

    auto tr_x_xxz_xyy = pbuffer.data(idx_dip_ff + 23);

    auto tr_x_xxz_xyz = pbuffer.data(idx_dip_ff + 24);

    auto tr_x_xxz_xzz = pbuffer.data(idx_dip_ff + 25);

    auto tr_x_xxz_yyy = pbuffer.data(idx_dip_ff + 26);

    auto tr_x_xxz_yyz = pbuffer.data(idx_dip_ff + 27);

    auto tr_x_xxz_yzz = pbuffer.data(idx_dip_ff + 28);

    auto tr_x_xxz_zzz = pbuffer.data(idx_dip_ff + 29);

    auto tr_x_xyy_xxx = pbuffer.data(idx_dip_ff + 30);

    auto tr_x_xyy_xxy = pbuffer.data(idx_dip_ff + 31);

    auto tr_x_xyy_xxz = pbuffer.data(idx_dip_ff + 32);

    auto tr_x_xyy_xyy = pbuffer.data(idx_dip_ff + 33);

    auto tr_x_xyy_xyz = pbuffer.data(idx_dip_ff + 34);

    auto tr_x_xyy_xzz = pbuffer.data(idx_dip_ff + 35);

    auto tr_x_xyy_yyy = pbuffer.data(idx_dip_ff + 36);

    auto tr_x_xyy_yyz = pbuffer.data(idx_dip_ff + 37);

    auto tr_x_xyy_yzz = pbuffer.data(idx_dip_ff + 38);

    auto tr_x_xyz_xxz = pbuffer.data(idx_dip_ff + 42);

    auto tr_x_xyz_xzz = pbuffer.data(idx_dip_ff + 45);

    auto tr_x_xzz_xxx = pbuffer.data(idx_dip_ff + 50);

    auto tr_x_xzz_xxy = pbuffer.data(idx_dip_ff + 51);

    auto tr_x_xzz_xxz = pbuffer.data(idx_dip_ff + 52);

    auto tr_x_xzz_xyy = pbuffer.data(idx_dip_ff + 53);

    auto tr_x_xzz_xyz = pbuffer.data(idx_dip_ff + 54);

    auto tr_x_xzz_xzz = pbuffer.data(idx_dip_ff + 55);

    auto tr_x_xzz_yyz = pbuffer.data(idx_dip_ff + 57);

    auto tr_x_xzz_yzz = pbuffer.data(idx_dip_ff + 58);

    auto tr_x_xzz_zzz = pbuffer.data(idx_dip_ff + 59);

    auto tr_x_yyy_xxx = pbuffer.data(idx_dip_ff + 60);

    auto tr_x_yyy_xxy = pbuffer.data(idx_dip_ff + 61);

    auto tr_x_yyy_xxz = pbuffer.data(idx_dip_ff + 62);

    auto tr_x_yyy_xyy = pbuffer.data(idx_dip_ff + 63);

    auto tr_x_yyy_xyz = pbuffer.data(idx_dip_ff + 64);

    auto tr_x_yyy_xzz = pbuffer.data(idx_dip_ff + 65);

    auto tr_x_yyy_yyy = pbuffer.data(idx_dip_ff + 66);

    auto tr_x_yyy_yyz = pbuffer.data(idx_dip_ff + 67);

    auto tr_x_yyy_yzz = pbuffer.data(idx_dip_ff + 68);

    auto tr_x_yyy_zzz = pbuffer.data(idx_dip_ff + 69);

    auto tr_x_yyz_xxy = pbuffer.data(idx_dip_ff + 71);

    auto tr_x_yyz_xxz = pbuffer.data(idx_dip_ff + 72);

    auto tr_x_yyz_xyy = pbuffer.data(idx_dip_ff + 73);

    auto tr_x_yyz_xzz = pbuffer.data(idx_dip_ff + 75);

    auto tr_x_yyz_yyy = pbuffer.data(idx_dip_ff + 76);

    auto tr_x_yyz_yyz = pbuffer.data(idx_dip_ff + 77);

    auto tr_x_yyz_yzz = pbuffer.data(idx_dip_ff + 78);

    auto tr_x_yyz_zzz = pbuffer.data(idx_dip_ff + 79);

    auto tr_x_yzz_xxx = pbuffer.data(idx_dip_ff + 80);

    auto tr_x_yzz_xxz = pbuffer.data(idx_dip_ff + 82);

    auto tr_x_yzz_xyz = pbuffer.data(idx_dip_ff + 84);

    auto tr_x_yzz_xzz = pbuffer.data(idx_dip_ff + 85);

    auto tr_x_yzz_yyy = pbuffer.data(idx_dip_ff + 86);

    auto tr_x_yzz_yyz = pbuffer.data(idx_dip_ff + 87);

    auto tr_x_yzz_yzz = pbuffer.data(idx_dip_ff + 88);

    auto tr_x_yzz_zzz = pbuffer.data(idx_dip_ff + 89);

    auto tr_x_zzz_xxx = pbuffer.data(idx_dip_ff + 90);

    auto tr_x_zzz_xxy = pbuffer.data(idx_dip_ff + 91);

    auto tr_x_zzz_xxz = pbuffer.data(idx_dip_ff + 92);

    auto tr_x_zzz_xyy = pbuffer.data(idx_dip_ff + 93);

    auto tr_x_zzz_xyz = pbuffer.data(idx_dip_ff + 94);

    auto tr_x_zzz_xzz = pbuffer.data(idx_dip_ff + 95);

    auto tr_x_zzz_yyy = pbuffer.data(idx_dip_ff + 96);

    auto tr_x_zzz_yyz = pbuffer.data(idx_dip_ff + 97);

    auto tr_x_zzz_yzz = pbuffer.data(idx_dip_ff + 98);

    auto tr_x_zzz_zzz = pbuffer.data(idx_dip_ff + 99);

    auto tr_y_xxx_xxx = pbuffer.data(idx_dip_ff + 100);

    auto tr_y_xxx_xxy = pbuffer.data(idx_dip_ff + 101);

    auto tr_y_xxx_xxz = pbuffer.data(idx_dip_ff + 102);

    auto tr_y_xxx_xyy = pbuffer.data(idx_dip_ff + 103);

    auto tr_y_xxx_xyz = pbuffer.data(idx_dip_ff + 104);

    auto tr_y_xxx_xzz = pbuffer.data(idx_dip_ff + 105);

    auto tr_y_xxx_yyy = pbuffer.data(idx_dip_ff + 106);

    auto tr_y_xxx_yyz = pbuffer.data(idx_dip_ff + 107);

    auto tr_y_xxx_yzz = pbuffer.data(idx_dip_ff + 108);

    auto tr_y_xxx_zzz = pbuffer.data(idx_dip_ff + 109);

    auto tr_y_xxy_xxx = pbuffer.data(idx_dip_ff + 110);

    auto tr_y_xxy_xxy = pbuffer.data(idx_dip_ff + 111);

    auto tr_y_xxy_xyy = pbuffer.data(idx_dip_ff + 113);

    auto tr_y_xxy_xyz = pbuffer.data(idx_dip_ff + 114);

    auto tr_y_xxy_yyy = pbuffer.data(idx_dip_ff + 116);

    auto tr_y_xxy_yyz = pbuffer.data(idx_dip_ff + 117);

    auto tr_y_xxy_yzz = pbuffer.data(idx_dip_ff + 118);

    auto tr_y_xxy_zzz = pbuffer.data(idx_dip_ff + 119);

    auto tr_y_xxz_xxx = pbuffer.data(idx_dip_ff + 120);

    auto tr_y_xxz_xxy = pbuffer.data(idx_dip_ff + 121);

    auto tr_y_xxz_xxz = pbuffer.data(idx_dip_ff + 122);

    auto tr_y_xxz_xyy = pbuffer.data(idx_dip_ff + 123);

    auto tr_y_xxz_xzz = pbuffer.data(idx_dip_ff + 125);

    auto tr_y_xxz_yyz = pbuffer.data(idx_dip_ff + 127);

    auto tr_y_xxz_yzz = pbuffer.data(idx_dip_ff + 128);

    auto tr_y_xxz_zzz = pbuffer.data(idx_dip_ff + 129);

    auto tr_y_xyy_xxx = pbuffer.data(idx_dip_ff + 130);

    auto tr_y_xyy_xxy = pbuffer.data(idx_dip_ff + 131);

    auto tr_y_xyy_xxz = pbuffer.data(idx_dip_ff + 132);

    auto tr_y_xyy_xyy = pbuffer.data(idx_dip_ff + 133);

    auto tr_y_xyy_xyz = pbuffer.data(idx_dip_ff + 134);

    auto tr_y_xyy_xzz = pbuffer.data(idx_dip_ff + 135);

    auto tr_y_xyy_yyy = pbuffer.data(idx_dip_ff + 136);

    auto tr_y_xyy_yyz = pbuffer.data(idx_dip_ff + 137);

    auto tr_y_xyy_yzz = pbuffer.data(idx_dip_ff + 138);

    auto tr_y_xyy_zzz = pbuffer.data(idx_dip_ff + 139);

    auto tr_y_xyz_yyz = pbuffer.data(idx_dip_ff + 147);

    auto tr_y_xyz_yzz = pbuffer.data(idx_dip_ff + 148);

    auto tr_y_xyz_zzz = pbuffer.data(idx_dip_ff + 149);

    auto tr_y_xzz_xxz = pbuffer.data(idx_dip_ff + 152);

    auto tr_y_xzz_xyz = pbuffer.data(idx_dip_ff + 154);

    auto tr_y_xzz_xzz = pbuffer.data(idx_dip_ff + 155);

    auto tr_y_xzz_yyy = pbuffer.data(idx_dip_ff + 156);

    auto tr_y_xzz_yyz = pbuffer.data(idx_dip_ff + 157);

    auto tr_y_xzz_yzz = pbuffer.data(idx_dip_ff + 158);

    auto tr_y_xzz_zzz = pbuffer.data(idx_dip_ff + 159);

    auto tr_y_yyy_xxx = pbuffer.data(idx_dip_ff + 160);

    auto tr_y_yyy_xxy = pbuffer.data(idx_dip_ff + 161);

    auto tr_y_yyy_xxz = pbuffer.data(idx_dip_ff + 162);

    auto tr_y_yyy_xyy = pbuffer.data(idx_dip_ff + 163);

    auto tr_y_yyy_xyz = pbuffer.data(idx_dip_ff + 164);

    auto tr_y_yyy_xzz = pbuffer.data(idx_dip_ff + 165);

    auto tr_y_yyy_yyy = pbuffer.data(idx_dip_ff + 166);

    auto tr_y_yyy_yyz = pbuffer.data(idx_dip_ff + 167);

    auto tr_y_yyy_yzz = pbuffer.data(idx_dip_ff + 168);

    auto tr_y_yyy_zzz = pbuffer.data(idx_dip_ff + 169);

    auto tr_y_yyz_xxx = pbuffer.data(idx_dip_ff + 170);

    auto tr_y_yyz_xxy = pbuffer.data(idx_dip_ff + 171);

    auto tr_y_yyz_xxz = pbuffer.data(idx_dip_ff + 172);

    auto tr_y_yyz_xyy = pbuffer.data(idx_dip_ff + 173);

    auto tr_y_yyz_xyz = pbuffer.data(idx_dip_ff + 174);

    auto tr_y_yyz_xzz = pbuffer.data(idx_dip_ff + 175);

    auto tr_y_yyz_yyy = pbuffer.data(idx_dip_ff + 176);

    auto tr_y_yyz_yyz = pbuffer.data(idx_dip_ff + 177);

    auto tr_y_yyz_yzz = pbuffer.data(idx_dip_ff + 178);

    auto tr_y_yyz_zzz = pbuffer.data(idx_dip_ff + 179);

    auto tr_y_yzz_xxx = pbuffer.data(idx_dip_ff + 180);

    auto tr_y_yzz_xxy = pbuffer.data(idx_dip_ff + 181);

    auto tr_y_yzz_xxz = pbuffer.data(idx_dip_ff + 182);

    auto tr_y_yzz_xyy = pbuffer.data(idx_dip_ff + 183);

    auto tr_y_yzz_xyz = pbuffer.data(idx_dip_ff + 184);

    auto tr_y_yzz_xzz = pbuffer.data(idx_dip_ff + 185);

    auto tr_y_yzz_yyy = pbuffer.data(idx_dip_ff + 186);

    auto tr_y_yzz_yyz = pbuffer.data(idx_dip_ff + 187);

    auto tr_y_yzz_yzz = pbuffer.data(idx_dip_ff + 188);

    auto tr_y_yzz_zzz = pbuffer.data(idx_dip_ff + 189);

    auto tr_y_zzz_xxx = pbuffer.data(idx_dip_ff + 190);

    auto tr_y_zzz_xxy = pbuffer.data(idx_dip_ff + 191);

    auto tr_y_zzz_xxz = pbuffer.data(idx_dip_ff + 192);

    auto tr_y_zzz_xyy = pbuffer.data(idx_dip_ff + 193);

    auto tr_y_zzz_xyz = pbuffer.data(idx_dip_ff + 194);

    auto tr_y_zzz_xzz = pbuffer.data(idx_dip_ff + 195);

    auto tr_y_zzz_yyy = pbuffer.data(idx_dip_ff + 196);

    auto tr_y_zzz_yyz = pbuffer.data(idx_dip_ff + 197);

    auto tr_y_zzz_yzz = pbuffer.data(idx_dip_ff + 198);

    auto tr_y_zzz_zzz = pbuffer.data(idx_dip_ff + 199);

    auto tr_z_xxx_xxx = pbuffer.data(idx_dip_ff + 200);

    auto tr_z_xxx_xxy = pbuffer.data(idx_dip_ff + 201);

    auto tr_z_xxx_xxz = pbuffer.data(idx_dip_ff + 202);

    auto tr_z_xxx_xyy = pbuffer.data(idx_dip_ff + 203);

    auto tr_z_xxx_xyz = pbuffer.data(idx_dip_ff + 204);

    auto tr_z_xxx_xzz = pbuffer.data(idx_dip_ff + 205);

    auto tr_z_xxx_yyy = pbuffer.data(idx_dip_ff + 206);

    auto tr_z_xxx_yyz = pbuffer.data(idx_dip_ff + 207);

    auto tr_z_xxx_yzz = pbuffer.data(idx_dip_ff + 208);

    auto tr_z_xxx_zzz = pbuffer.data(idx_dip_ff + 209);

    auto tr_z_xxy_xxx = pbuffer.data(idx_dip_ff + 210);

    auto tr_z_xxy_xxz = pbuffer.data(idx_dip_ff + 212);

    auto tr_z_xxy_xzz = pbuffer.data(idx_dip_ff + 215);

    auto tr_z_xxy_yyy = pbuffer.data(idx_dip_ff + 216);

    auto tr_z_xxy_yyz = pbuffer.data(idx_dip_ff + 217);

    auto tr_z_xxy_yzz = pbuffer.data(idx_dip_ff + 218);

    auto tr_z_xxz_xxx = pbuffer.data(idx_dip_ff + 220);

    auto tr_z_xxz_xxy = pbuffer.data(idx_dip_ff + 221);

    auto tr_z_xxz_xxz = pbuffer.data(idx_dip_ff + 222);

    auto tr_z_xxz_xyy = pbuffer.data(idx_dip_ff + 223);

    auto tr_z_xxz_xyz = pbuffer.data(idx_dip_ff + 224);

    auto tr_z_xxz_xzz = pbuffer.data(idx_dip_ff + 225);

    auto tr_z_xxz_yyy = pbuffer.data(idx_dip_ff + 226);

    auto tr_z_xxz_yyz = pbuffer.data(idx_dip_ff + 227);

    auto tr_z_xxz_yzz = pbuffer.data(idx_dip_ff + 228);

    auto tr_z_xxz_zzz = pbuffer.data(idx_dip_ff + 229);

    auto tr_z_xyy_xxy = pbuffer.data(idx_dip_ff + 231);

    auto tr_z_xyy_xyy = pbuffer.data(idx_dip_ff + 233);

    auto tr_z_xyy_xyz = pbuffer.data(idx_dip_ff + 234);

    auto tr_z_xyy_yyy = pbuffer.data(idx_dip_ff + 236);

    auto tr_z_xyy_yyz = pbuffer.data(idx_dip_ff + 237);

    auto tr_z_xyy_yzz = pbuffer.data(idx_dip_ff + 238);

    auto tr_z_xyy_zzz = pbuffer.data(idx_dip_ff + 239);

    auto tr_z_xyz_yyy = pbuffer.data(idx_dip_ff + 246);

    auto tr_z_xyz_yyz = pbuffer.data(idx_dip_ff + 247);

    auto tr_z_xyz_yzz = pbuffer.data(idx_dip_ff + 248);

    auto tr_z_xzz_xxx = pbuffer.data(idx_dip_ff + 250);

    auto tr_z_xzz_xxy = pbuffer.data(idx_dip_ff + 251);

    auto tr_z_xzz_xxz = pbuffer.data(idx_dip_ff + 252);

    auto tr_z_xzz_xyy = pbuffer.data(idx_dip_ff + 253);

    auto tr_z_xzz_xyz = pbuffer.data(idx_dip_ff + 254);

    auto tr_z_xzz_xzz = pbuffer.data(idx_dip_ff + 255);

    auto tr_z_xzz_yyy = pbuffer.data(idx_dip_ff + 256);

    auto tr_z_xzz_yyz = pbuffer.data(idx_dip_ff + 257);

    auto tr_z_xzz_yzz = pbuffer.data(idx_dip_ff + 258);

    auto tr_z_xzz_zzz = pbuffer.data(idx_dip_ff + 259);

    auto tr_z_yyy_xxx = pbuffer.data(idx_dip_ff + 260);

    auto tr_z_yyy_xxy = pbuffer.data(idx_dip_ff + 261);

    auto tr_z_yyy_xxz = pbuffer.data(idx_dip_ff + 262);

    auto tr_z_yyy_xyy = pbuffer.data(idx_dip_ff + 263);

    auto tr_z_yyy_xyz = pbuffer.data(idx_dip_ff + 264);

    auto tr_z_yyy_xzz = pbuffer.data(idx_dip_ff + 265);

    auto tr_z_yyy_yyy = pbuffer.data(idx_dip_ff + 266);

    auto tr_z_yyy_yyz = pbuffer.data(idx_dip_ff + 267);

    auto tr_z_yyy_yzz = pbuffer.data(idx_dip_ff + 268);

    auto tr_z_yyy_zzz = pbuffer.data(idx_dip_ff + 269);

    auto tr_z_yyz_xxx = pbuffer.data(idx_dip_ff + 270);

    auto tr_z_yyz_xxy = pbuffer.data(idx_dip_ff + 271);

    auto tr_z_yyz_xxz = pbuffer.data(idx_dip_ff + 272);

    auto tr_z_yyz_xyy = pbuffer.data(idx_dip_ff + 273);

    auto tr_z_yyz_xyz = pbuffer.data(idx_dip_ff + 274);

    auto tr_z_yyz_xzz = pbuffer.data(idx_dip_ff + 275);

    auto tr_z_yyz_yyy = pbuffer.data(idx_dip_ff + 276);

    auto tr_z_yyz_yyz = pbuffer.data(idx_dip_ff + 277);

    auto tr_z_yyz_yzz = pbuffer.data(idx_dip_ff + 278);

    auto tr_z_yyz_zzz = pbuffer.data(idx_dip_ff + 279);

    auto tr_z_yzz_xxx = pbuffer.data(idx_dip_ff + 280);

    auto tr_z_yzz_xxy = pbuffer.data(idx_dip_ff + 281);

    auto tr_z_yzz_xxz = pbuffer.data(idx_dip_ff + 282);

    auto tr_z_yzz_xyy = pbuffer.data(idx_dip_ff + 283);

    auto tr_z_yzz_xyz = pbuffer.data(idx_dip_ff + 284);

    auto tr_z_yzz_xzz = pbuffer.data(idx_dip_ff + 285);

    auto tr_z_yzz_yyy = pbuffer.data(idx_dip_ff + 286);

    auto tr_z_yzz_yyz = pbuffer.data(idx_dip_ff + 287);

    auto tr_z_yzz_yzz = pbuffer.data(idx_dip_ff + 288);

    auto tr_z_yzz_zzz = pbuffer.data(idx_dip_ff + 289);

    auto tr_z_zzz_xxx = pbuffer.data(idx_dip_ff + 290);

    auto tr_z_zzz_xxy = pbuffer.data(idx_dip_ff + 291);

    auto tr_z_zzz_xxz = pbuffer.data(idx_dip_ff + 292);

    auto tr_z_zzz_xyy = pbuffer.data(idx_dip_ff + 293);

    auto tr_z_zzz_xyz = pbuffer.data(idx_dip_ff + 294);

    auto tr_z_zzz_xzz = pbuffer.data(idx_dip_ff + 295);

    auto tr_z_zzz_yyy = pbuffer.data(idx_dip_ff + 296);

    auto tr_z_zzz_yyz = pbuffer.data(idx_dip_ff + 297);

    auto tr_z_zzz_yzz = pbuffer.data(idx_dip_ff + 298);

    auto tr_z_zzz_zzz = pbuffer.data(idx_dip_ff + 299);

    // Set up 0-10 components of targeted buffer : GF

    auto tr_x_xxxx_xxx = pbuffer.data(idx_dip_gf);

    auto tr_x_xxxx_xxy = pbuffer.data(idx_dip_gf + 1);

    auto tr_x_xxxx_xxz = pbuffer.data(idx_dip_gf + 2);

    auto tr_x_xxxx_xyy = pbuffer.data(idx_dip_gf + 3);

    auto tr_x_xxxx_xyz = pbuffer.data(idx_dip_gf + 4);

    auto tr_x_xxxx_xzz = pbuffer.data(idx_dip_gf + 5);

    auto tr_x_xxxx_yyy = pbuffer.data(idx_dip_gf + 6);

    auto tr_x_xxxx_yyz = pbuffer.data(idx_dip_gf + 7);

    auto tr_x_xxxx_yzz = pbuffer.data(idx_dip_gf + 8);

    auto tr_x_xxxx_zzz = pbuffer.data(idx_dip_gf + 9);

#pragma omp simd aligned(pa_x,              \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xzz,   \
                             tr_x_xx_yyy,   \
                             tr_x_xx_yyz,   \
                             tr_x_xx_yzz,   \
                             tr_x_xx_zzz,   \
                             tr_x_xxx_xx,   \
                             tr_x_xxx_xxx,  \
                             tr_x_xxx_xxy,  \
                             tr_x_xxx_xxz,  \
                             tr_x_xxx_xy,   \
                             tr_x_xxx_xyy,  \
                             tr_x_xxx_xyz,  \
                             tr_x_xxx_xz,   \
                             tr_x_xxx_xzz,  \
                             tr_x_xxx_yy,   \
                             tr_x_xxx_yyy,  \
                             tr_x_xxx_yyz,  \
                             tr_x_xxx_yz,   \
                             tr_x_xxx_yzz,  \
                             tr_x_xxx_zz,   \
                             tr_x_xxx_zzz,  \
                             tr_x_xxxx_xxx, \
                             tr_x_xxxx_xxy, \
                             tr_x_xxxx_xxz, \
                             tr_x_xxxx_xyy, \
                             tr_x_xxxx_xyz, \
                             tr_x_xxxx_xzz, \
                             tr_x_xxxx_yyy, \
                             tr_x_xxxx_yyz, \
                             tr_x_xxxx_yzz, \
                             tr_x_xxxx_zzz, \
                             ts_xxx_xxx,    \
                             ts_xxx_xxy,    \
                             ts_xxx_xxz,    \
                             ts_xxx_xyy,    \
                             ts_xxx_xyz,    \
                             ts_xxx_xzz,    \
                             ts_xxx_yyy,    \
                             ts_xxx_yyz,    \
                             ts_xxx_yzz,    \
                             ts_xxx_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_xxx[i] = 3.0 * tr_x_xx_xxx[i] * fe_0 + 3.0 * tr_x_xxx_xx[i] * fe_0 + ts_xxx_xxx[i] * fe_0 + tr_x_xxx_xxx[i] * pa_x[i];

        tr_x_xxxx_xxy[i] = 3.0 * tr_x_xx_xxy[i] * fe_0 + 2.0 * tr_x_xxx_xy[i] * fe_0 + ts_xxx_xxy[i] * fe_0 + tr_x_xxx_xxy[i] * pa_x[i];

        tr_x_xxxx_xxz[i] = 3.0 * tr_x_xx_xxz[i] * fe_0 + 2.0 * tr_x_xxx_xz[i] * fe_0 + ts_xxx_xxz[i] * fe_0 + tr_x_xxx_xxz[i] * pa_x[i];

        tr_x_xxxx_xyy[i] = 3.0 * tr_x_xx_xyy[i] * fe_0 + tr_x_xxx_yy[i] * fe_0 + ts_xxx_xyy[i] * fe_0 + tr_x_xxx_xyy[i] * pa_x[i];

        tr_x_xxxx_xyz[i] = 3.0 * tr_x_xx_xyz[i] * fe_0 + tr_x_xxx_yz[i] * fe_0 + ts_xxx_xyz[i] * fe_0 + tr_x_xxx_xyz[i] * pa_x[i];

        tr_x_xxxx_xzz[i] = 3.0 * tr_x_xx_xzz[i] * fe_0 + tr_x_xxx_zz[i] * fe_0 + ts_xxx_xzz[i] * fe_0 + tr_x_xxx_xzz[i] * pa_x[i];

        tr_x_xxxx_yyy[i] = 3.0 * tr_x_xx_yyy[i] * fe_0 + ts_xxx_yyy[i] * fe_0 + tr_x_xxx_yyy[i] * pa_x[i];

        tr_x_xxxx_yyz[i] = 3.0 * tr_x_xx_yyz[i] * fe_0 + ts_xxx_yyz[i] * fe_0 + tr_x_xxx_yyz[i] * pa_x[i];

        tr_x_xxxx_yzz[i] = 3.0 * tr_x_xx_yzz[i] * fe_0 + ts_xxx_yzz[i] * fe_0 + tr_x_xxx_yzz[i] * pa_x[i];

        tr_x_xxxx_zzz[i] = 3.0 * tr_x_xx_zzz[i] * fe_0 + ts_xxx_zzz[i] * fe_0 + tr_x_xxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto tr_x_xxxy_xxx = pbuffer.data(idx_dip_gf + 10);

    auto tr_x_xxxy_xxy = pbuffer.data(idx_dip_gf + 11);

    auto tr_x_xxxy_xxz = pbuffer.data(idx_dip_gf + 12);

    auto tr_x_xxxy_xyy = pbuffer.data(idx_dip_gf + 13);

    auto tr_x_xxxy_xyz = pbuffer.data(idx_dip_gf + 14);

    auto tr_x_xxxy_xzz = pbuffer.data(idx_dip_gf + 15);

    auto tr_x_xxxy_yyy = pbuffer.data(idx_dip_gf + 16);

    auto tr_x_xxxy_yyz = pbuffer.data(idx_dip_gf + 17);

    auto tr_x_xxxy_yzz = pbuffer.data(idx_dip_gf + 18);

    auto tr_x_xxxy_zzz = pbuffer.data(idx_dip_gf + 19);

#pragma omp simd aligned(pa_y,              \
                             tr_x_xxx_xx,   \
                             tr_x_xxx_xxx,  \
                             tr_x_xxx_xxy,  \
                             tr_x_xxx_xxz,  \
                             tr_x_xxx_xy,   \
                             tr_x_xxx_xyy,  \
                             tr_x_xxx_xyz,  \
                             tr_x_xxx_xz,   \
                             tr_x_xxx_xzz,  \
                             tr_x_xxx_yy,   \
                             tr_x_xxx_yyy,  \
                             tr_x_xxx_yyz,  \
                             tr_x_xxx_yz,   \
                             tr_x_xxx_yzz,  \
                             tr_x_xxx_zz,   \
                             tr_x_xxx_zzz,  \
                             tr_x_xxxy_xxx, \
                             tr_x_xxxy_xxy, \
                             tr_x_xxxy_xxz, \
                             tr_x_xxxy_xyy, \
                             tr_x_xxxy_xyz, \
                             tr_x_xxxy_xzz, \
                             tr_x_xxxy_yyy, \
                             tr_x_xxxy_yyz, \
                             tr_x_xxxy_yzz, \
                             tr_x_xxxy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_xxx[i] = tr_x_xxx_xxx[i] * pa_y[i];

        tr_x_xxxy_xxy[i] = tr_x_xxx_xx[i] * fe_0 + tr_x_xxx_xxy[i] * pa_y[i];

        tr_x_xxxy_xxz[i] = tr_x_xxx_xxz[i] * pa_y[i];

        tr_x_xxxy_xyy[i] = 2.0 * tr_x_xxx_xy[i] * fe_0 + tr_x_xxx_xyy[i] * pa_y[i];

        tr_x_xxxy_xyz[i] = tr_x_xxx_xz[i] * fe_0 + tr_x_xxx_xyz[i] * pa_y[i];

        tr_x_xxxy_xzz[i] = tr_x_xxx_xzz[i] * pa_y[i];

        tr_x_xxxy_yyy[i] = 3.0 * tr_x_xxx_yy[i] * fe_0 + tr_x_xxx_yyy[i] * pa_y[i];

        tr_x_xxxy_yyz[i] = 2.0 * tr_x_xxx_yz[i] * fe_0 + tr_x_xxx_yyz[i] * pa_y[i];

        tr_x_xxxy_yzz[i] = tr_x_xxx_zz[i] * fe_0 + tr_x_xxx_yzz[i] * pa_y[i];

        tr_x_xxxy_zzz[i] = tr_x_xxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto tr_x_xxxz_xxx = pbuffer.data(idx_dip_gf + 20);

    auto tr_x_xxxz_xxy = pbuffer.data(idx_dip_gf + 21);

    auto tr_x_xxxz_xxz = pbuffer.data(idx_dip_gf + 22);

    auto tr_x_xxxz_xyy = pbuffer.data(idx_dip_gf + 23);

    auto tr_x_xxxz_xyz = pbuffer.data(idx_dip_gf + 24);

    auto tr_x_xxxz_xzz = pbuffer.data(idx_dip_gf + 25);

    auto tr_x_xxxz_yyy = pbuffer.data(idx_dip_gf + 26);

    auto tr_x_xxxz_yyz = pbuffer.data(idx_dip_gf + 27);

    auto tr_x_xxxz_yzz = pbuffer.data(idx_dip_gf + 28);

    auto tr_x_xxxz_zzz = pbuffer.data(idx_dip_gf + 29);

#pragma omp simd aligned(pa_z,              \
                             tr_x_xxx_xx,   \
                             tr_x_xxx_xxx,  \
                             tr_x_xxx_xxy,  \
                             tr_x_xxx_xxz,  \
                             tr_x_xxx_xy,   \
                             tr_x_xxx_xyy,  \
                             tr_x_xxx_xyz,  \
                             tr_x_xxx_xz,   \
                             tr_x_xxx_xzz,  \
                             tr_x_xxx_yy,   \
                             tr_x_xxx_yyy,  \
                             tr_x_xxx_yyz,  \
                             tr_x_xxx_yz,   \
                             tr_x_xxx_yzz,  \
                             tr_x_xxx_zz,   \
                             tr_x_xxx_zzz,  \
                             tr_x_xxxz_xxx, \
                             tr_x_xxxz_xxy, \
                             tr_x_xxxz_xxz, \
                             tr_x_xxxz_xyy, \
                             tr_x_xxxz_xyz, \
                             tr_x_xxxz_xzz, \
                             tr_x_xxxz_yyy, \
                             tr_x_xxxz_yyz, \
                             tr_x_xxxz_yzz, \
                             tr_x_xxxz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_xxx[i] = tr_x_xxx_xxx[i] * pa_z[i];

        tr_x_xxxz_xxy[i] = tr_x_xxx_xxy[i] * pa_z[i];

        tr_x_xxxz_xxz[i] = tr_x_xxx_xx[i] * fe_0 + tr_x_xxx_xxz[i] * pa_z[i];

        tr_x_xxxz_xyy[i] = tr_x_xxx_xyy[i] * pa_z[i];

        tr_x_xxxz_xyz[i] = tr_x_xxx_xy[i] * fe_0 + tr_x_xxx_xyz[i] * pa_z[i];

        tr_x_xxxz_xzz[i] = 2.0 * tr_x_xxx_xz[i] * fe_0 + tr_x_xxx_xzz[i] * pa_z[i];

        tr_x_xxxz_yyy[i] = tr_x_xxx_yyy[i] * pa_z[i];

        tr_x_xxxz_yyz[i] = tr_x_xxx_yy[i] * fe_0 + tr_x_xxx_yyz[i] * pa_z[i];

        tr_x_xxxz_yzz[i] = 2.0 * tr_x_xxx_yz[i] * fe_0 + tr_x_xxx_yzz[i] * pa_z[i];

        tr_x_xxxz_zzz[i] = 3.0 * tr_x_xxx_zz[i] * fe_0 + tr_x_xxx_zzz[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto tr_x_xxyy_xxx = pbuffer.data(idx_dip_gf + 30);

    auto tr_x_xxyy_xxy = pbuffer.data(idx_dip_gf + 31);

    auto tr_x_xxyy_xxz = pbuffer.data(idx_dip_gf + 32);

    auto tr_x_xxyy_xyy = pbuffer.data(idx_dip_gf + 33);

    auto tr_x_xxyy_xyz = pbuffer.data(idx_dip_gf + 34);

    auto tr_x_xxyy_xzz = pbuffer.data(idx_dip_gf + 35);

    auto tr_x_xxyy_yyy = pbuffer.data(idx_dip_gf + 36);

    auto tr_x_xxyy_yyz = pbuffer.data(idx_dip_gf + 37);

    auto tr_x_xxyy_yzz = pbuffer.data(idx_dip_gf + 38);

    auto tr_x_xxyy_zzz = pbuffer.data(idx_dip_gf + 39);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xzz,   \
                             tr_x_xx_zzz,   \
                             tr_x_xxy_xx,   \
                             tr_x_xxy_xxx,  \
                             tr_x_xxy_xxy,  \
                             tr_x_xxy_xxz,  \
                             tr_x_xxy_xy,   \
                             tr_x_xxy_xyy,  \
                             tr_x_xxy_xyz,  \
                             tr_x_xxy_xz,   \
                             tr_x_xxy_xzz,  \
                             tr_x_xxy_zzz,  \
                             tr_x_xxyy_xxx, \
                             tr_x_xxyy_xxy, \
                             tr_x_xxyy_xxz, \
                             tr_x_xxyy_xyy, \
                             tr_x_xxyy_xyz, \
                             tr_x_xxyy_xzz, \
                             tr_x_xxyy_yyy, \
                             tr_x_xxyy_yyz, \
                             tr_x_xxyy_yzz, \
                             tr_x_xxyy_zzz, \
                             tr_x_xyy_yyy,  \
                             tr_x_xyy_yyz,  \
                             tr_x_xyy_yzz,  \
                             tr_x_yy_yyy,   \
                             tr_x_yy_yyz,   \
                             tr_x_yy_yzz,   \
                             ts_xyy_yyy,    \
                             ts_xyy_yyz,    \
                             ts_xyy_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_xxx[i] = tr_x_xx_xxx[i] * fe_0 + tr_x_xxy_xxx[i] * pa_y[i];

        tr_x_xxyy_xxy[i] = tr_x_xx_xxy[i] * fe_0 + tr_x_xxy_xx[i] * fe_0 + tr_x_xxy_xxy[i] * pa_y[i];

        tr_x_xxyy_xxz[i] = tr_x_xx_xxz[i] * fe_0 + tr_x_xxy_xxz[i] * pa_y[i];

        tr_x_xxyy_xyy[i] = tr_x_xx_xyy[i] * fe_0 + 2.0 * tr_x_xxy_xy[i] * fe_0 + tr_x_xxy_xyy[i] * pa_y[i];

        tr_x_xxyy_xyz[i] = tr_x_xx_xyz[i] * fe_0 + tr_x_xxy_xz[i] * fe_0 + tr_x_xxy_xyz[i] * pa_y[i];

        tr_x_xxyy_xzz[i] = tr_x_xx_xzz[i] * fe_0 + tr_x_xxy_xzz[i] * pa_y[i];

        tr_x_xxyy_yyy[i] = tr_x_yy_yyy[i] * fe_0 + ts_xyy_yyy[i] * fe_0 + tr_x_xyy_yyy[i] * pa_x[i];

        tr_x_xxyy_yyz[i] = tr_x_yy_yyz[i] * fe_0 + ts_xyy_yyz[i] * fe_0 + tr_x_xyy_yyz[i] * pa_x[i];

        tr_x_xxyy_yzz[i] = tr_x_yy_yzz[i] * fe_0 + ts_xyy_yzz[i] * fe_0 + tr_x_xyy_yzz[i] * pa_x[i];

        tr_x_xxyy_zzz[i] = tr_x_xx_zzz[i] * fe_0 + tr_x_xxy_zzz[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto tr_x_xxyz_xxx = pbuffer.data(idx_dip_gf + 40);

    auto tr_x_xxyz_xxy = pbuffer.data(idx_dip_gf + 41);

    auto tr_x_xxyz_xxz = pbuffer.data(idx_dip_gf + 42);

    auto tr_x_xxyz_xyy = pbuffer.data(idx_dip_gf + 43);

    auto tr_x_xxyz_xyz = pbuffer.data(idx_dip_gf + 44);

    auto tr_x_xxyz_xzz = pbuffer.data(idx_dip_gf + 45);

    auto tr_x_xxyz_yyy = pbuffer.data(idx_dip_gf + 46);

    auto tr_x_xxyz_yyz = pbuffer.data(idx_dip_gf + 47);

    auto tr_x_xxyz_yzz = pbuffer.data(idx_dip_gf + 48);

    auto tr_x_xxyz_zzz = pbuffer.data(idx_dip_gf + 49);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_xxy_xxy,  \
                             tr_x_xxy_xyy,  \
                             tr_x_xxy_yyy,  \
                             tr_x_xxyz_xxx, \
                             tr_x_xxyz_xxy, \
                             tr_x_xxyz_xxz, \
                             tr_x_xxyz_xyy, \
                             tr_x_xxyz_xyz, \
                             tr_x_xxyz_xzz, \
                             tr_x_xxyz_yyy, \
                             tr_x_xxyz_yyz, \
                             tr_x_xxyz_yzz, \
                             tr_x_xxyz_zzz, \
                             tr_x_xxz_xxx,  \
                             tr_x_xxz_xxz,  \
                             tr_x_xxz_xyz,  \
                             tr_x_xxz_xz,   \
                             tr_x_xxz_xzz,  \
                             tr_x_xxz_yyz,  \
                             tr_x_xxz_yz,   \
                             tr_x_xxz_yzz,  \
                             tr_x_xxz_zz,   \
                             tr_x_xxz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyz_xxx[i] = tr_x_xxz_xxx[i] * pa_y[i];

        tr_x_xxyz_xxy[i] = tr_x_xxy_xxy[i] * pa_z[i];

        tr_x_xxyz_xxz[i] = tr_x_xxz_xxz[i] * pa_y[i];

        tr_x_xxyz_xyy[i] = tr_x_xxy_xyy[i] * pa_z[i];

        tr_x_xxyz_xyz[i] = tr_x_xxz_xz[i] * fe_0 + tr_x_xxz_xyz[i] * pa_y[i];

        tr_x_xxyz_xzz[i] = tr_x_xxz_xzz[i] * pa_y[i];

        tr_x_xxyz_yyy[i] = tr_x_xxy_yyy[i] * pa_z[i];

        tr_x_xxyz_yyz[i] = 2.0 * tr_x_xxz_yz[i] * fe_0 + tr_x_xxz_yyz[i] * pa_y[i];

        tr_x_xxyz_yzz[i] = tr_x_xxz_zz[i] * fe_0 + tr_x_xxz_yzz[i] * pa_y[i];

        tr_x_xxyz_zzz[i] = tr_x_xxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto tr_x_xxzz_xxx = pbuffer.data(idx_dip_gf + 50);

    auto tr_x_xxzz_xxy = pbuffer.data(idx_dip_gf + 51);

    auto tr_x_xxzz_xxz = pbuffer.data(idx_dip_gf + 52);

    auto tr_x_xxzz_xyy = pbuffer.data(idx_dip_gf + 53);

    auto tr_x_xxzz_xyz = pbuffer.data(idx_dip_gf + 54);

    auto tr_x_xxzz_xzz = pbuffer.data(idx_dip_gf + 55);

    auto tr_x_xxzz_yyy = pbuffer.data(idx_dip_gf + 56);

    auto tr_x_xxzz_yyz = pbuffer.data(idx_dip_gf + 57);

    auto tr_x_xxzz_yzz = pbuffer.data(idx_dip_gf + 58);

    auto tr_x_xxzz_zzz = pbuffer.data(idx_dip_gf + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xx_xxx,   \
                             tr_x_xx_xxy,   \
                             tr_x_xx_xxz,   \
                             tr_x_xx_xyy,   \
                             tr_x_xx_xyz,   \
                             tr_x_xx_xzz,   \
                             tr_x_xx_yyy,   \
                             tr_x_xxz_xx,   \
                             tr_x_xxz_xxx,  \
                             tr_x_xxz_xxy,  \
                             tr_x_xxz_xxz,  \
                             tr_x_xxz_xy,   \
                             tr_x_xxz_xyy,  \
                             tr_x_xxz_xyz,  \
                             tr_x_xxz_xz,   \
                             tr_x_xxz_xzz,  \
                             tr_x_xxz_yyy,  \
                             tr_x_xxzz_xxx, \
                             tr_x_xxzz_xxy, \
                             tr_x_xxzz_xxz, \
                             tr_x_xxzz_xyy, \
                             tr_x_xxzz_xyz, \
                             tr_x_xxzz_xzz, \
                             tr_x_xxzz_yyy, \
                             tr_x_xxzz_yyz, \
                             tr_x_xxzz_yzz, \
                             tr_x_xxzz_zzz, \
                             tr_x_xzz_yyz,  \
                             tr_x_xzz_yzz,  \
                             tr_x_xzz_zzz,  \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yzz,   \
                             tr_x_zz_zzz,   \
                             ts_xzz_yyz,    \
                             ts_xzz_yzz,    \
                             ts_xzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_xxx[i] = tr_x_xx_xxx[i] * fe_0 + tr_x_xxz_xxx[i] * pa_z[i];

        tr_x_xxzz_xxy[i] = tr_x_xx_xxy[i] * fe_0 + tr_x_xxz_xxy[i] * pa_z[i];

        tr_x_xxzz_xxz[i] = tr_x_xx_xxz[i] * fe_0 + tr_x_xxz_xx[i] * fe_0 + tr_x_xxz_xxz[i] * pa_z[i];

        tr_x_xxzz_xyy[i] = tr_x_xx_xyy[i] * fe_0 + tr_x_xxz_xyy[i] * pa_z[i];

        tr_x_xxzz_xyz[i] = tr_x_xx_xyz[i] * fe_0 + tr_x_xxz_xy[i] * fe_0 + tr_x_xxz_xyz[i] * pa_z[i];

        tr_x_xxzz_xzz[i] = tr_x_xx_xzz[i] * fe_0 + 2.0 * tr_x_xxz_xz[i] * fe_0 + tr_x_xxz_xzz[i] * pa_z[i];

        tr_x_xxzz_yyy[i] = tr_x_xx_yyy[i] * fe_0 + tr_x_xxz_yyy[i] * pa_z[i];

        tr_x_xxzz_yyz[i] = tr_x_zz_yyz[i] * fe_0 + ts_xzz_yyz[i] * fe_0 + tr_x_xzz_yyz[i] * pa_x[i];

        tr_x_xxzz_yzz[i] = tr_x_zz_yzz[i] * fe_0 + ts_xzz_yzz[i] * fe_0 + tr_x_xzz_yzz[i] * pa_x[i];

        tr_x_xxzz_zzz[i] = tr_x_zz_zzz[i] * fe_0 + ts_xzz_zzz[i] * fe_0 + tr_x_xzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto tr_x_xyyy_xxx = pbuffer.data(idx_dip_gf + 60);

    auto tr_x_xyyy_xxy = pbuffer.data(idx_dip_gf + 61);

    auto tr_x_xyyy_xxz = pbuffer.data(idx_dip_gf + 62);

    auto tr_x_xyyy_xyy = pbuffer.data(idx_dip_gf + 63);

    auto tr_x_xyyy_xyz = pbuffer.data(idx_dip_gf + 64);

    auto tr_x_xyyy_xzz = pbuffer.data(idx_dip_gf + 65);

    auto tr_x_xyyy_yyy = pbuffer.data(idx_dip_gf + 66);

    auto tr_x_xyyy_yyz = pbuffer.data(idx_dip_gf + 67);

    auto tr_x_xyyy_yzz = pbuffer.data(idx_dip_gf + 68);

    auto tr_x_xyyy_zzz = pbuffer.data(idx_dip_gf + 69);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xy_xxx,   \
                             tr_x_xy_xxz,   \
                             tr_x_xy_xzz,   \
                             tr_x_xyy_xxx,  \
                             tr_x_xyy_xxz,  \
                             tr_x_xyy_xzz,  \
                             tr_x_xyyy_xxx, \
                             tr_x_xyyy_xxy, \
                             tr_x_xyyy_xxz, \
                             tr_x_xyyy_xyy, \
                             tr_x_xyyy_xyz, \
                             tr_x_xyyy_xzz, \
                             tr_x_xyyy_yyy, \
                             tr_x_xyyy_yyz, \
                             tr_x_xyyy_yzz, \
                             tr_x_xyyy_zzz, \
                             tr_x_yyy_xxy,  \
                             tr_x_yyy_xy,   \
                             tr_x_yyy_xyy,  \
                             tr_x_yyy_xyz,  \
                             tr_x_yyy_yy,   \
                             tr_x_yyy_yyy,  \
                             tr_x_yyy_yyz,  \
                             tr_x_yyy_yz,   \
                             tr_x_yyy_yzz,  \
                             tr_x_yyy_zzz,  \
                             ts_yyy_xxy,    \
                             ts_yyy_xyy,    \
                             ts_yyy_xyz,    \
                             ts_yyy_yyy,    \
                             ts_yyy_yyz,    \
                             ts_yyy_yzz,    \
                             ts_yyy_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_xxx[i] = 2.0 * tr_x_xy_xxx[i] * fe_0 + tr_x_xyy_xxx[i] * pa_y[i];

        tr_x_xyyy_xxy[i] = 2.0 * tr_x_yyy_xy[i] * fe_0 + ts_yyy_xxy[i] * fe_0 + tr_x_yyy_xxy[i] * pa_x[i];

        tr_x_xyyy_xxz[i] = 2.0 * tr_x_xy_xxz[i] * fe_0 + tr_x_xyy_xxz[i] * pa_y[i];

        tr_x_xyyy_xyy[i] = tr_x_yyy_yy[i] * fe_0 + ts_yyy_xyy[i] * fe_0 + tr_x_yyy_xyy[i] * pa_x[i];

        tr_x_xyyy_xyz[i] = tr_x_yyy_yz[i] * fe_0 + ts_yyy_xyz[i] * fe_0 + tr_x_yyy_xyz[i] * pa_x[i];

        tr_x_xyyy_xzz[i] = 2.0 * tr_x_xy_xzz[i] * fe_0 + tr_x_xyy_xzz[i] * pa_y[i];

        tr_x_xyyy_yyy[i] = ts_yyy_yyy[i] * fe_0 + tr_x_yyy_yyy[i] * pa_x[i];

        tr_x_xyyy_yyz[i] = ts_yyy_yyz[i] * fe_0 + tr_x_yyy_yyz[i] * pa_x[i];

        tr_x_xyyy_yzz[i] = ts_yyy_yzz[i] * fe_0 + tr_x_yyy_yzz[i] * pa_x[i];

        tr_x_xyyy_zzz[i] = ts_yyy_zzz[i] * fe_0 + tr_x_yyy_zzz[i] * pa_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto tr_x_xyyz_xxx = pbuffer.data(idx_dip_gf + 70);

    auto tr_x_xyyz_xxy = pbuffer.data(idx_dip_gf + 71);

    auto tr_x_xyyz_xxz = pbuffer.data(idx_dip_gf + 72);

    auto tr_x_xyyz_xyy = pbuffer.data(idx_dip_gf + 73);

    auto tr_x_xyyz_xyz = pbuffer.data(idx_dip_gf + 74);

    auto tr_x_xyyz_xzz = pbuffer.data(idx_dip_gf + 75);

    auto tr_x_xyyz_yyy = pbuffer.data(idx_dip_gf + 76);

    auto tr_x_xyyz_yyz = pbuffer.data(idx_dip_gf + 77);

    auto tr_x_xyyz_yzz = pbuffer.data(idx_dip_gf + 78);

    auto tr_x_xyyz_zzz = pbuffer.data(idx_dip_gf + 79);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             tr_x_xyy_xxx,  \
                             tr_x_xyy_xxy,  \
                             tr_x_xyy_xy,   \
                             tr_x_xyy_xyy,  \
                             tr_x_xyy_xyz,  \
                             tr_x_xyy_yyy,  \
                             tr_x_xyyz_xxx, \
                             tr_x_xyyz_xxy, \
                             tr_x_xyyz_xxz, \
                             tr_x_xyyz_xyy, \
                             tr_x_xyyz_xyz, \
                             tr_x_xyyz_xzz, \
                             tr_x_xyyz_yyy, \
                             tr_x_xyyz_yyz, \
                             tr_x_xyyz_yzz, \
                             tr_x_xyyz_zzz, \
                             tr_x_xyz_xxz,  \
                             tr_x_xyz_xzz,  \
                             tr_x_xz_xxz,   \
                             tr_x_xz_xzz,   \
                             tr_x_yyz_yyz,  \
                             tr_x_yyz_yzz,  \
                             tr_x_yyz_zzz,  \
                             ts_yyz_yyz,    \
                             ts_yyz_yzz,    \
                             ts_yyz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_xxx[i] = tr_x_xyy_xxx[i] * pa_z[i];

        tr_x_xyyz_xxy[i] = tr_x_xyy_xxy[i] * pa_z[i];

        tr_x_xyyz_xxz[i] = tr_x_xz_xxz[i] * fe_0 + tr_x_xyz_xxz[i] * pa_y[i];

        tr_x_xyyz_xyy[i] = tr_x_xyy_xyy[i] * pa_z[i];

        tr_x_xyyz_xyz[i] = tr_x_xyy_xy[i] * fe_0 + tr_x_xyy_xyz[i] * pa_z[i];

        tr_x_xyyz_xzz[i] = tr_x_xz_xzz[i] * fe_0 + tr_x_xyz_xzz[i] * pa_y[i];

        tr_x_xyyz_yyy[i] = tr_x_xyy_yyy[i] * pa_z[i];

        tr_x_xyyz_yyz[i] = ts_yyz_yyz[i] * fe_0 + tr_x_yyz_yyz[i] * pa_x[i];

        tr_x_xyyz_yzz[i] = ts_yyz_yzz[i] * fe_0 + tr_x_yyz_yzz[i] * pa_x[i];

        tr_x_xyyz_zzz[i] = ts_yyz_zzz[i] * fe_0 + tr_x_yyz_zzz[i] * pa_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto tr_x_xyzz_xxx = pbuffer.data(idx_dip_gf + 80);

    auto tr_x_xyzz_xxy = pbuffer.data(idx_dip_gf + 81);

    auto tr_x_xyzz_xxz = pbuffer.data(idx_dip_gf + 82);

    auto tr_x_xyzz_xyy = pbuffer.data(idx_dip_gf + 83);

    auto tr_x_xyzz_xyz = pbuffer.data(idx_dip_gf + 84);

    auto tr_x_xyzz_xzz = pbuffer.data(idx_dip_gf + 85);

    auto tr_x_xyzz_yyy = pbuffer.data(idx_dip_gf + 86);

    auto tr_x_xyzz_yyz = pbuffer.data(idx_dip_gf + 87);

    auto tr_x_xyzz_yzz = pbuffer.data(idx_dip_gf + 88);

    auto tr_x_xyzz_zzz = pbuffer.data(idx_dip_gf + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_x_xyzz_xxx, \
                             tr_x_xyzz_xxy, \
                             tr_x_xyzz_xxz, \
                             tr_x_xyzz_xyy, \
                             tr_x_xyzz_xyz, \
                             tr_x_xyzz_xzz, \
                             tr_x_xyzz_yyy, \
                             tr_x_xyzz_yyz, \
                             tr_x_xyzz_yzz, \
                             tr_x_xyzz_zzz, \
                             tr_x_xzz_xx,   \
                             tr_x_xzz_xxx,  \
                             tr_x_xzz_xxy,  \
                             tr_x_xzz_xxz,  \
                             tr_x_xzz_xy,   \
                             tr_x_xzz_xyy,  \
                             tr_x_xzz_xyz,  \
                             tr_x_xzz_xz,   \
                             tr_x_xzz_xzz,  \
                             tr_x_xzz_zzz,  \
                             tr_x_yzz_yyy,  \
                             tr_x_yzz_yyz,  \
                             tr_x_yzz_yzz,  \
                             ts_yzz_yyy,    \
                             ts_yzz_yyz,    \
                             ts_yzz_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_xxx[i] = tr_x_xzz_xxx[i] * pa_y[i];

        tr_x_xyzz_xxy[i] = tr_x_xzz_xx[i] * fe_0 + tr_x_xzz_xxy[i] * pa_y[i];

        tr_x_xyzz_xxz[i] = tr_x_xzz_xxz[i] * pa_y[i];

        tr_x_xyzz_xyy[i] = 2.0 * tr_x_xzz_xy[i] * fe_0 + tr_x_xzz_xyy[i] * pa_y[i];

        tr_x_xyzz_xyz[i] = tr_x_xzz_xz[i] * fe_0 + tr_x_xzz_xyz[i] * pa_y[i];

        tr_x_xyzz_xzz[i] = tr_x_xzz_xzz[i] * pa_y[i];

        tr_x_xyzz_yyy[i] = ts_yzz_yyy[i] * fe_0 + tr_x_yzz_yyy[i] * pa_x[i];

        tr_x_xyzz_yyz[i] = ts_yzz_yyz[i] * fe_0 + tr_x_yzz_yyz[i] * pa_x[i];

        tr_x_xyzz_yzz[i] = ts_yzz_yzz[i] * fe_0 + tr_x_yzz_yzz[i] * pa_x[i];

        tr_x_xyzz_zzz[i] = tr_x_xzz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto tr_x_xzzz_xxx = pbuffer.data(idx_dip_gf + 90);

    auto tr_x_xzzz_xxy = pbuffer.data(idx_dip_gf + 91);

    auto tr_x_xzzz_xxz = pbuffer.data(idx_dip_gf + 92);

    auto tr_x_xzzz_xyy = pbuffer.data(idx_dip_gf + 93);

    auto tr_x_xzzz_xyz = pbuffer.data(idx_dip_gf + 94);

    auto tr_x_xzzz_xzz = pbuffer.data(idx_dip_gf + 95);

    auto tr_x_xzzz_yyy = pbuffer.data(idx_dip_gf + 96);

    auto tr_x_xzzz_yyz = pbuffer.data(idx_dip_gf + 97);

    auto tr_x_xzzz_yzz = pbuffer.data(idx_dip_gf + 98);

    auto tr_x_xzzz_zzz = pbuffer.data(idx_dip_gf + 99);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_x_xz_xxx,   \
                             tr_x_xz_xxy,   \
                             tr_x_xz_xyy,   \
                             tr_x_xzz_xxx,  \
                             tr_x_xzz_xxy,  \
                             tr_x_xzz_xyy,  \
                             tr_x_xzzz_xxx, \
                             tr_x_xzzz_xxy, \
                             tr_x_xzzz_xxz, \
                             tr_x_xzzz_xyy, \
                             tr_x_xzzz_xyz, \
                             tr_x_xzzz_xzz, \
                             tr_x_xzzz_yyy, \
                             tr_x_xzzz_yyz, \
                             tr_x_xzzz_yzz, \
                             tr_x_xzzz_zzz, \
                             tr_x_zzz_xxz,  \
                             tr_x_zzz_xyz,  \
                             tr_x_zzz_xz,   \
                             tr_x_zzz_xzz,  \
                             tr_x_zzz_yyy,  \
                             tr_x_zzz_yyz,  \
                             tr_x_zzz_yz,   \
                             tr_x_zzz_yzz,  \
                             tr_x_zzz_zz,   \
                             tr_x_zzz_zzz,  \
                             ts_zzz_xxz,    \
                             ts_zzz_xyz,    \
                             ts_zzz_xzz,    \
                             ts_zzz_yyy,    \
                             ts_zzz_yyz,    \
                             ts_zzz_yzz,    \
                             ts_zzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_xxx[i] = 2.0 * tr_x_xz_xxx[i] * fe_0 + tr_x_xzz_xxx[i] * pa_z[i];

        tr_x_xzzz_xxy[i] = 2.0 * tr_x_xz_xxy[i] * fe_0 + tr_x_xzz_xxy[i] * pa_z[i];

        tr_x_xzzz_xxz[i] = 2.0 * tr_x_zzz_xz[i] * fe_0 + ts_zzz_xxz[i] * fe_0 + tr_x_zzz_xxz[i] * pa_x[i];

        tr_x_xzzz_xyy[i] = 2.0 * tr_x_xz_xyy[i] * fe_0 + tr_x_xzz_xyy[i] * pa_z[i];

        tr_x_xzzz_xyz[i] = tr_x_zzz_yz[i] * fe_0 + ts_zzz_xyz[i] * fe_0 + tr_x_zzz_xyz[i] * pa_x[i];

        tr_x_xzzz_xzz[i] = tr_x_zzz_zz[i] * fe_0 + ts_zzz_xzz[i] * fe_0 + tr_x_zzz_xzz[i] * pa_x[i];

        tr_x_xzzz_yyy[i] = ts_zzz_yyy[i] * fe_0 + tr_x_zzz_yyy[i] * pa_x[i];

        tr_x_xzzz_yyz[i] = ts_zzz_yyz[i] * fe_0 + tr_x_zzz_yyz[i] * pa_x[i];

        tr_x_xzzz_yzz[i] = ts_zzz_yzz[i] * fe_0 + tr_x_zzz_yzz[i] * pa_x[i];

        tr_x_xzzz_zzz[i] = ts_zzz_zzz[i] * fe_0 + tr_x_zzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto tr_x_yyyy_xxx = pbuffer.data(idx_dip_gf + 100);

    auto tr_x_yyyy_xxy = pbuffer.data(idx_dip_gf + 101);

    auto tr_x_yyyy_xxz = pbuffer.data(idx_dip_gf + 102);

    auto tr_x_yyyy_xyy = pbuffer.data(idx_dip_gf + 103);

    auto tr_x_yyyy_xyz = pbuffer.data(idx_dip_gf + 104);

    auto tr_x_yyyy_xzz = pbuffer.data(idx_dip_gf + 105);

    auto tr_x_yyyy_yyy = pbuffer.data(idx_dip_gf + 106);

    auto tr_x_yyyy_yyz = pbuffer.data(idx_dip_gf + 107);

    auto tr_x_yyyy_yzz = pbuffer.data(idx_dip_gf + 108);

    auto tr_x_yyyy_zzz = pbuffer.data(idx_dip_gf + 109);

#pragma omp simd aligned(pa_y,              \
                             tr_x_yy_xxx,   \
                             tr_x_yy_xxy,   \
                             tr_x_yy_xxz,   \
                             tr_x_yy_xyy,   \
                             tr_x_yy_xyz,   \
                             tr_x_yy_xzz,   \
                             tr_x_yy_yyy,   \
                             tr_x_yy_yyz,   \
                             tr_x_yy_yzz,   \
                             tr_x_yy_zzz,   \
                             tr_x_yyy_xx,   \
                             tr_x_yyy_xxx,  \
                             tr_x_yyy_xxy,  \
                             tr_x_yyy_xxz,  \
                             tr_x_yyy_xy,   \
                             tr_x_yyy_xyy,  \
                             tr_x_yyy_xyz,  \
                             tr_x_yyy_xz,   \
                             tr_x_yyy_xzz,  \
                             tr_x_yyy_yy,   \
                             tr_x_yyy_yyy,  \
                             tr_x_yyy_yyz,  \
                             tr_x_yyy_yz,   \
                             tr_x_yyy_yzz,  \
                             tr_x_yyy_zz,   \
                             tr_x_yyy_zzz,  \
                             tr_x_yyyy_xxx, \
                             tr_x_yyyy_xxy, \
                             tr_x_yyyy_xxz, \
                             tr_x_yyyy_xyy, \
                             tr_x_yyyy_xyz, \
                             tr_x_yyyy_xzz, \
                             tr_x_yyyy_yyy, \
                             tr_x_yyyy_yyz, \
                             tr_x_yyyy_yzz, \
                             tr_x_yyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_xxx[i] = 3.0 * tr_x_yy_xxx[i] * fe_0 + tr_x_yyy_xxx[i] * pa_y[i];

        tr_x_yyyy_xxy[i] = 3.0 * tr_x_yy_xxy[i] * fe_0 + tr_x_yyy_xx[i] * fe_0 + tr_x_yyy_xxy[i] * pa_y[i];

        tr_x_yyyy_xxz[i] = 3.0 * tr_x_yy_xxz[i] * fe_0 + tr_x_yyy_xxz[i] * pa_y[i];

        tr_x_yyyy_xyy[i] = 3.0 * tr_x_yy_xyy[i] * fe_0 + 2.0 * tr_x_yyy_xy[i] * fe_0 + tr_x_yyy_xyy[i] * pa_y[i];

        tr_x_yyyy_xyz[i] = 3.0 * tr_x_yy_xyz[i] * fe_0 + tr_x_yyy_xz[i] * fe_0 + tr_x_yyy_xyz[i] * pa_y[i];

        tr_x_yyyy_xzz[i] = 3.0 * tr_x_yy_xzz[i] * fe_0 + tr_x_yyy_xzz[i] * pa_y[i];

        tr_x_yyyy_yyy[i] = 3.0 * tr_x_yy_yyy[i] * fe_0 + 3.0 * tr_x_yyy_yy[i] * fe_0 + tr_x_yyy_yyy[i] * pa_y[i];

        tr_x_yyyy_yyz[i] = 3.0 * tr_x_yy_yyz[i] * fe_0 + 2.0 * tr_x_yyy_yz[i] * fe_0 + tr_x_yyy_yyz[i] * pa_y[i];

        tr_x_yyyy_yzz[i] = 3.0 * tr_x_yy_yzz[i] * fe_0 + tr_x_yyy_zz[i] * fe_0 + tr_x_yyy_yzz[i] * pa_y[i];

        tr_x_yyyy_zzz[i] = 3.0 * tr_x_yy_zzz[i] * fe_0 + tr_x_yyy_zzz[i] * pa_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto tr_x_yyyz_xxx = pbuffer.data(idx_dip_gf + 110);

    auto tr_x_yyyz_xxy = pbuffer.data(idx_dip_gf + 111);

    auto tr_x_yyyz_xxz = pbuffer.data(idx_dip_gf + 112);

    auto tr_x_yyyz_xyy = pbuffer.data(idx_dip_gf + 113);

    auto tr_x_yyyz_xyz = pbuffer.data(idx_dip_gf + 114);

    auto tr_x_yyyz_xzz = pbuffer.data(idx_dip_gf + 115);

    auto tr_x_yyyz_yyy = pbuffer.data(idx_dip_gf + 116);

    auto tr_x_yyyz_yyz = pbuffer.data(idx_dip_gf + 117);

    auto tr_x_yyyz_yzz = pbuffer.data(idx_dip_gf + 118);

    auto tr_x_yyyz_zzz = pbuffer.data(idx_dip_gf + 119);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yyy_xxx,  \
                             tr_x_yyy_xxy,  \
                             tr_x_yyy_xy,   \
                             tr_x_yyy_xyy,  \
                             tr_x_yyy_xyz,  \
                             tr_x_yyy_yy,   \
                             tr_x_yyy_yyy,  \
                             tr_x_yyy_yyz,  \
                             tr_x_yyy_yz,   \
                             tr_x_yyy_yzz,  \
                             tr_x_yyyz_xxx, \
                             tr_x_yyyz_xxy, \
                             tr_x_yyyz_xxz, \
                             tr_x_yyyz_xyy, \
                             tr_x_yyyz_xyz, \
                             tr_x_yyyz_xzz, \
                             tr_x_yyyz_yyy, \
                             tr_x_yyyz_yyz, \
                             tr_x_yyyz_yzz, \
                             tr_x_yyyz_zzz, \
                             tr_x_yyz_xxz,  \
                             tr_x_yyz_xzz,  \
                             tr_x_yyz_zzz,  \
                             tr_x_yz_xxz,   \
                             tr_x_yz_xzz,   \
                             tr_x_yz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_xxx[i] = tr_x_yyy_xxx[i] * pa_z[i];

        tr_x_yyyz_xxy[i] = tr_x_yyy_xxy[i] * pa_z[i];

        tr_x_yyyz_xxz[i] = 2.0 * tr_x_yz_xxz[i] * fe_0 + tr_x_yyz_xxz[i] * pa_y[i];

        tr_x_yyyz_xyy[i] = tr_x_yyy_xyy[i] * pa_z[i];

        tr_x_yyyz_xyz[i] = tr_x_yyy_xy[i] * fe_0 + tr_x_yyy_xyz[i] * pa_z[i];

        tr_x_yyyz_xzz[i] = 2.0 * tr_x_yz_xzz[i] * fe_0 + tr_x_yyz_xzz[i] * pa_y[i];

        tr_x_yyyz_yyy[i] = tr_x_yyy_yyy[i] * pa_z[i];

        tr_x_yyyz_yyz[i] = tr_x_yyy_yy[i] * fe_0 + tr_x_yyy_yyz[i] * pa_z[i];

        tr_x_yyyz_yzz[i] = 2.0 * tr_x_yyy_yz[i] * fe_0 + tr_x_yyy_yzz[i] * pa_z[i];

        tr_x_yyyz_zzz[i] = 2.0 * tr_x_yz_zzz[i] * fe_0 + tr_x_yyz_zzz[i] * pa_y[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto tr_x_yyzz_xxx = pbuffer.data(idx_dip_gf + 120);

    auto tr_x_yyzz_xxy = pbuffer.data(idx_dip_gf + 121);

    auto tr_x_yyzz_xxz = pbuffer.data(idx_dip_gf + 122);

    auto tr_x_yyzz_xyy = pbuffer.data(idx_dip_gf + 123);

    auto tr_x_yyzz_xyz = pbuffer.data(idx_dip_gf + 124);

    auto tr_x_yyzz_xzz = pbuffer.data(idx_dip_gf + 125);

    auto tr_x_yyzz_yyy = pbuffer.data(idx_dip_gf + 126);

    auto tr_x_yyzz_yyz = pbuffer.data(idx_dip_gf + 127);

    auto tr_x_yyzz_yzz = pbuffer.data(idx_dip_gf + 128);

    auto tr_x_yyzz_zzz = pbuffer.data(idx_dip_gf + 129);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_x_yy_xxy,   \
                             tr_x_yy_xyy,   \
                             tr_x_yy_yyy,   \
                             tr_x_yyz_xxy,  \
                             tr_x_yyz_xyy,  \
                             tr_x_yyz_yyy,  \
                             tr_x_yyzz_xxx, \
                             tr_x_yyzz_xxy, \
                             tr_x_yyzz_xxz, \
                             tr_x_yyzz_xyy, \
                             tr_x_yyzz_xyz, \
                             tr_x_yyzz_xzz, \
                             tr_x_yyzz_yyy, \
                             tr_x_yyzz_yyz, \
                             tr_x_yyzz_yzz, \
                             tr_x_yyzz_zzz, \
                             tr_x_yzz_xxx,  \
                             tr_x_yzz_xxz,  \
                             tr_x_yzz_xyz,  \
                             tr_x_yzz_xz,   \
                             tr_x_yzz_xzz,  \
                             tr_x_yzz_yyz,  \
                             tr_x_yzz_yz,   \
                             tr_x_yzz_yzz,  \
                             tr_x_yzz_zz,   \
                             tr_x_yzz_zzz,  \
                             tr_x_zz_xxx,   \
                             tr_x_zz_xxz,   \
                             tr_x_zz_xyz,   \
                             tr_x_zz_xzz,   \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yzz,   \
                             tr_x_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_xxx[i] = tr_x_zz_xxx[i] * fe_0 + tr_x_yzz_xxx[i] * pa_y[i];

        tr_x_yyzz_xxy[i] = tr_x_yy_xxy[i] * fe_0 + tr_x_yyz_xxy[i] * pa_z[i];

        tr_x_yyzz_xxz[i] = tr_x_zz_xxz[i] * fe_0 + tr_x_yzz_xxz[i] * pa_y[i];

        tr_x_yyzz_xyy[i] = tr_x_yy_xyy[i] * fe_0 + tr_x_yyz_xyy[i] * pa_z[i];

        tr_x_yyzz_xyz[i] = tr_x_zz_xyz[i] * fe_0 + tr_x_yzz_xz[i] * fe_0 + tr_x_yzz_xyz[i] * pa_y[i];

        tr_x_yyzz_xzz[i] = tr_x_zz_xzz[i] * fe_0 + tr_x_yzz_xzz[i] * pa_y[i];

        tr_x_yyzz_yyy[i] = tr_x_yy_yyy[i] * fe_0 + tr_x_yyz_yyy[i] * pa_z[i];

        tr_x_yyzz_yyz[i] = tr_x_zz_yyz[i] * fe_0 + 2.0 * tr_x_yzz_yz[i] * fe_0 + tr_x_yzz_yyz[i] * pa_y[i];

        tr_x_yyzz_yzz[i] = tr_x_zz_yzz[i] * fe_0 + tr_x_yzz_zz[i] * fe_0 + tr_x_yzz_yzz[i] * pa_y[i];

        tr_x_yyzz_zzz[i] = tr_x_zz_zzz[i] * fe_0 + tr_x_yzz_zzz[i] * pa_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto tr_x_yzzz_xxx = pbuffer.data(idx_dip_gf + 130);

    auto tr_x_yzzz_xxy = pbuffer.data(idx_dip_gf + 131);

    auto tr_x_yzzz_xxz = pbuffer.data(idx_dip_gf + 132);

    auto tr_x_yzzz_xyy = pbuffer.data(idx_dip_gf + 133);

    auto tr_x_yzzz_xyz = pbuffer.data(idx_dip_gf + 134);

    auto tr_x_yzzz_xzz = pbuffer.data(idx_dip_gf + 135);

    auto tr_x_yzzz_yyy = pbuffer.data(idx_dip_gf + 136);

    auto tr_x_yzzz_yyz = pbuffer.data(idx_dip_gf + 137);

    auto tr_x_yzzz_yzz = pbuffer.data(idx_dip_gf + 138);

    auto tr_x_yzzz_zzz = pbuffer.data(idx_dip_gf + 139);

#pragma omp simd aligned(pa_y,              \
                             tr_x_yzzz_xxx, \
                             tr_x_yzzz_xxy, \
                             tr_x_yzzz_xxz, \
                             tr_x_yzzz_xyy, \
                             tr_x_yzzz_xyz, \
                             tr_x_yzzz_xzz, \
                             tr_x_yzzz_yyy, \
                             tr_x_yzzz_yyz, \
                             tr_x_yzzz_yzz, \
                             tr_x_yzzz_zzz, \
                             tr_x_zzz_xx,   \
                             tr_x_zzz_xxx,  \
                             tr_x_zzz_xxy,  \
                             tr_x_zzz_xxz,  \
                             tr_x_zzz_xy,   \
                             tr_x_zzz_xyy,  \
                             tr_x_zzz_xyz,  \
                             tr_x_zzz_xz,   \
                             tr_x_zzz_xzz,  \
                             tr_x_zzz_yy,   \
                             tr_x_zzz_yyy,  \
                             tr_x_zzz_yyz,  \
                             tr_x_zzz_yz,   \
                             tr_x_zzz_yzz,  \
                             tr_x_zzz_zz,   \
                             tr_x_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_xxx[i] = tr_x_zzz_xxx[i] * pa_y[i];

        tr_x_yzzz_xxy[i] = tr_x_zzz_xx[i] * fe_0 + tr_x_zzz_xxy[i] * pa_y[i];

        tr_x_yzzz_xxz[i] = tr_x_zzz_xxz[i] * pa_y[i];

        tr_x_yzzz_xyy[i] = 2.0 * tr_x_zzz_xy[i] * fe_0 + tr_x_zzz_xyy[i] * pa_y[i];

        tr_x_yzzz_xyz[i] = tr_x_zzz_xz[i] * fe_0 + tr_x_zzz_xyz[i] * pa_y[i];

        tr_x_yzzz_xzz[i] = tr_x_zzz_xzz[i] * pa_y[i];

        tr_x_yzzz_yyy[i] = 3.0 * tr_x_zzz_yy[i] * fe_0 + tr_x_zzz_yyy[i] * pa_y[i];

        tr_x_yzzz_yyz[i] = 2.0 * tr_x_zzz_yz[i] * fe_0 + tr_x_zzz_yyz[i] * pa_y[i];

        tr_x_yzzz_yzz[i] = tr_x_zzz_zz[i] * fe_0 + tr_x_zzz_yzz[i] * pa_y[i];

        tr_x_yzzz_zzz[i] = tr_x_zzz_zzz[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto tr_x_zzzz_xxx = pbuffer.data(idx_dip_gf + 140);

    auto tr_x_zzzz_xxy = pbuffer.data(idx_dip_gf + 141);

    auto tr_x_zzzz_xxz = pbuffer.data(idx_dip_gf + 142);

    auto tr_x_zzzz_xyy = pbuffer.data(idx_dip_gf + 143);

    auto tr_x_zzzz_xyz = pbuffer.data(idx_dip_gf + 144);

    auto tr_x_zzzz_xzz = pbuffer.data(idx_dip_gf + 145);

    auto tr_x_zzzz_yyy = pbuffer.data(idx_dip_gf + 146);

    auto tr_x_zzzz_yyz = pbuffer.data(idx_dip_gf + 147);

    auto tr_x_zzzz_yzz = pbuffer.data(idx_dip_gf + 148);

    auto tr_x_zzzz_zzz = pbuffer.data(idx_dip_gf + 149);

#pragma omp simd aligned(pa_z,              \
                             tr_x_zz_xxx,   \
                             tr_x_zz_xxy,   \
                             tr_x_zz_xxz,   \
                             tr_x_zz_xyy,   \
                             tr_x_zz_xyz,   \
                             tr_x_zz_xzz,   \
                             tr_x_zz_yyy,   \
                             tr_x_zz_yyz,   \
                             tr_x_zz_yzz,   \
                             tr_x_zz_zzz,   \
                             tr_x_zzz_xx,   \
                             tr_x_zzz_xxx,  \
                             tr_x_zzz_xxy,  \
                             tr_x_zzz_xxz,  \
                             tr_x_zzz_xy,   \
                             tr_x_zzz_xyy,  \
                             tr_x_zzz_xyz,  \
                             tr_x_zzz_xz,   \
                             tr_x_zzz_xzz,  \
                             tr_x_zzz_yy,   \
                             tr_x_zzz_yyy,  \
                             tr_x_zzz_yyz,  \
                             tr_x_zzz_yz,   \
                             tr_x_zzz_yzz,  \
                             tr_x_zzz_zz,   \
                             tr_x_zzz_zzz,  \
                             tr_x_zzzz_xxx, \
                             tr_x_zzzz_xxy, \
                             tr_x_zzzz_xxz, \
                             tr_x_zzzz_xyy, \
                             tr_x_zzzz_xyz, \
                             tr_x_zzzz_xzz, \
                             tr_x_zzzz_yyy, \
                             tr_x_zzzz_yyz, \
                             tr_x_zzzz_yzz, \
                             tr_x_zzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_xxx[i] = 3.0 * tr_x_zz_xxx[i] * fe_0 + tr_x_zzz_xxx[i] * pa_z[i];

        tr_x_zzzz_xxy[i] = 3.0 * tr_x_zz_xxy[i] * fe_0 + tr_x_zzz_xxy[i] * pa_z[i];

        tr_x_zzzz_xxz[i] = 3.0 * tr_x_zz_xxz[i] * fe_0 + tr_x_zzz_xx[i] * fe_0 + tr_x_zzz_xxz[i] * pa_z[i];

        tr_x_zzzz_xyy[i] = 3.0 * tr_x_zz_xyy[i] * fe_0 + tr_x_zzz_xyy[i] * pa_z[i];

        tr_x_zzzz_xyz[i] = 3.0 * tr_x_zz_xyz[i] * fe_0 + tr_x_zzz_xy[i] * fe_0 + tr_x_zzz_xyz[i] * pa_z[i];

        tr_x_zzzz_xzz[i] = 3.0 * tr_x_zz_xzz[i] * fe_0 + 2.0 * tr_x_zzz_xz[i] * fe_0 + tr_x_zzz_xzz[i] * pa_z[i];

        tr_x_zzzz_yyy[i] = 3.0 * tr_x_zz_yyy[i] * fe_0 + tr_x_zzz_yyy[i] * pa_z[i];

        tr_x_zzzz_yyz[i] = 3.0 * tr_x_zz_yyz[i] * fe_0 + tr_x_zzz_yy[i] * fe_0 + tr_x_zzz_yyz[i] * pa_z[i];

        tr_x_zzzz_yzz[i] = 3.0 * tr_x_zz_yzz[i] * fe_0 + 2.0 * tr_x_zzz_yz[i] * fe_0 + tr_x_zzz_yzz[i] * pa_z[i];

        tr_x_zzzz_zzz[i] = 3.0 * tr_x_zz_zzz[i] * fe_0 + 3.0 * tr_x_zzz_zz[i] * fe_0 + tr_x_zzz_zzz[i] * pa_z[i];
    }

    // Set up 150-160 components of targeted buffer : GF

    auto tr_y_xxxx_xxx = pbuffer.data(idx_dip_gf + 150);

    auto tr_y_xxxx_xxy = pbuffer.data(idx_dip_gf + 151);

    auto tr_y_xxxx_xxz = pbuffer.data(idx_dip_gf + 152);

    auto tr_y_xxxx_xyy = pbuffer.data(idx_dip_gf + 153);

    auto tr_y_xxxx_xyz = pbuffer.data(idx_dip_gf + 154);

    auto tr_y_xxxx_xzz = pbuffer.data(idx_dip_gf + 155);

    auto tr_y_xxxx_yyy = pbuffer.data(idx_dip_gf + 156);

    auto tr_y_xxxx_yyz = pbuffer.data(idx_dip_gf + 157);

    auto tr_y_xxxx_yzz = pbuffer.data(idx_dip_gf + 158);

    auto tr_y_xxxx_zzz = pbuffer.data(idx_dip_gf + 159);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xx_xxx,   \
                             tr_y_xx_xxy,   \
                             tr_y_xx_xxz,   \
                             tr_y_xx_xyy,   \
                             tr_y_xx_xyz,   \
                             tr_y_xx_xzz,   \
                             tr_y_xx_yyy,   \
                             tr_y_xx_yyz,   \
                             tr_y_xx_yzz,   \
                             tr_y_xx_zzz,   \
                             tr_y_xxx_xx,   \
                             tr_y_xxx_xxx,  \
                             tr_y_xxx_xxy,  \
                             tr_y_xxx_xxz,  \
                             tr_y_xxx_xy,   \
                             tr_y_xxx_xyy,  \
                             tr_y_xxx_xyz,  \
                             tr_y_xxx_xz,   \
                             tr_y_xxx_xzz,  \
                             tr_y_xxx_yy,   \
                             tr_y_xxx_yyy,  \
                             tr_y_xxx_yyz,  \
                             tr_y_xxx_yz,   \
                             tr_y_xxx_yzz,  \
                             tr_y_xxx_zz,   \
                             tr_y_xxx_zzz,  \
                             tr_y_xxxx_xxx, \
                             tr_y_xxxx_xxy, \
                             tr_y_xxxx_xxz, \
                             tr_y_xxxx_xyy, \
                             tr_y_xxxx_xyz, \
                             tr_y_xxxx_xzz, \
                             tr_y_xxxx_yyy, \
                             tr_y_xxxx_yyz, \
                             tr_y_xxxx_yzz, \
                             tr_y_xxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_xxx[i] = 3.0 * tr_y_xx_xxx[i] * fe_0 + 3.0 * tr_y_xxx_xx[i] * fe_0 + tr_y_xxx_xxx[i] * pa_x[i];

        tr_y_xxxx_xxy[i] = 3.0 * tr_y_xx_xxy[i] * fe_0 + 2.0 * tr_y_xxx_xy[i] * fe_0 + tr_y_xxx_xxy[i] * pa_x[i];

        tr_y_xxxx_xxz[i] = 3.0 * tr_y_xx_xxz[i] * fe_0 + 2.0 * tr_y_xxx_xz[i] * fe_0 + tr_y_xxx_xxz[i] * pa_x[i];

        tr_y_xxxx_xyy[i] = 3.0 * tr_y_xx_xyy[i] * fe_0 + tr_y_xxx_yy[i] * fe_0 + tr_y_xxx_xyy[i] * pa_x[i];

        tr_y_xxxx_xyz[i] = 3.0 * tr_y_xx_xyz[i] * fe_0 + tr_y_xxx_yz[i] * fe_0 + tr_y_xxx_xyz[i] * pa_x[i];

        tr_y_xxxx_xzz[i] = 3.0 * tr_y_xx_xzz[i] * fe_0 + tr_y_xxx_zz[i] * fe_0 + tr_y_xxx_xzz[i] * pa_x[i];

        tr_y_xxxx_yyy[i] = 3.0 * tr_y_xx_yyy[i] * fe_0 + tr_y_xxx_yyy[i] * pa_x[i];

        tr_y_xxxx_yyz[i] = 3.0 * tr_y_xx_yyz[i] * fe_0 + tr_y_xxx_yyz[i] * pa_x[i];

        tr_y_xxxx_yzz[i] = 3.0 * tr_y_xx_yzz[i] * fe_0 + tr_y_xxx_yzz[i] * pa_x[i];

        tr_y_xxxx_zzz[i] = 3.0 * tr_y_xx_zzz[i] * fe_0 + tr_y_xxx_zzz[i] * pa_x[i];
    }

    // Set up 160-170 components of targeted buffer : GF

    auto tr_y_xxxy_xxx = pbuffer.data(idx_dip_gf + 160);

    auto tr_y_xxxy_xxy = pbuffer.data(idx_dip_gf + 161);

    auto tr_y_xxxy_xxz = pbuffer.data(idx_dip_gf + 162);

    auto tr_y_xxxy_xyy = pbuffer.data(idx_dip_gf + 163);

    auto tr_y_xxxy_xyz = pbuffer.data(idx_dip_gf + 164);

    auto tr_y_xxxy_xzz = pbuffer.data(idx_dip_gf + 165);

    auto tr_y_xxxy_yyy = pbuffer.data(idx_dip_gf + 166);

    auto tr_y_xxxy_yyz = pbuffer.data(idx_dip_gf + 167);

    auto tr_y_xxxy_yzz = pbuffer.data(idx_dip_gf + 168);

    auto tr_y_xxxy_zzz = pbuffer.data(idx_dip_gf + 169);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_y_xxx_xxx,  \
                             tr_y_xxx_xxz,  \
                             tr_y_xxx_xzz,  \
                             tr_y_xxxy_xxx, \
                             tr_y_xxxy_xxy, \
                             tr_y_xxxy_xxz, \
                             tr_y_xxxy_xyy, \
                             tr_y_xxxy_xyz, \
                             tr_y_xxxy_xzz, \
                             tr_y_xxxy_yyy, \
                             tr_y_xxxy_yyz, \
                             tr_y_xxxy_yzz, \
                             tr_y_xxxy_zzz, \
                             tr_y_xxy_xxy,  \
                             tr_y_xxy_xy,   \
                             tr_y_xxy_xyy,  \
                             tr_y_xxy_xyz,  \
                             tr_y_xxy_yy,   \
                             tr_y_xxy_yyy,  \
                             tr_y_xxy_yyz,  \
                             tr_y_xxy_yz,   \
                             tr_y_xxy_yzz,  \
                             tr_y_xxy_zzz,  \
                             tr_y_xy_xxy,   \
                             tr_y_xy_xyy,   \
                             tr_y_xy_xyz,   \
                             tr_y_xy_yyy,   \
                             tr_y_xy_yyz,   \
                             tr_y_xy_yzz,   \
                             tr_y_xy_zzz,   \
                             ts_xxx_xxx,    \
                             ts_xxx_xxz,    \
                             ts_xxx_xzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_xxx[i] = ts_xxx_xxx[i] * fe_0 + tr_y_xxx_xxx[i] * pa_y[i];

        tr_y_xxxy_xxy[i] = 2.0 * tr_y_xy_xxy[i] * fe_0 + 2.0 * tr_y_xxy_xy[i] * fe_0 + tr_y_xxy_xxy[i] * pa_x[i];

        tr_y_xxxy_xxz[i] = ts_xxx_xxz[i] * fe_0 + tr_y_xxx_xxz[i] * pa_y[i];

        tr_y_xxxy_xyy[i] = 2.0 * tr_y_xy_xyy[i] * fe_0 + tr_y_xxy_yy[i] * fe_0 + tr_y_xxy_xyy[i] * pa_x[i];

        tr_y_xxxy_xyz[i] = 2.0 * tr_y_xy_xyz[i] * fe_0 + tr_y_xxy_yz[i] * fe_0 + tr_y_xxy_xyz[i] * pa_x[i];

        tr_y_xxxy_xzz[i] = ts_xxx_xzz[i] * fe_0 + tr_y_xxx_xzz[i] * pa_y[i];

        tr_y_xxxy_yyy[i] = 2.0 * tr_y_xy_yyy[i] * fe_0 + tr_y_xxy_yyy[i] * pa_x[i];

        tr_y_xxxy_yyz[i] = 2.0 * tr_y_xy_yyz[i] * fe_0 + tr_y_xxy_yyz[i] * pa_x[i];

        tr_y_xxxy_yzz[i] = 2.0 * tr_y_xy_yzz[i] * fe_0 + tr_y_xxy_yzz[i] * pa_x[i];

        tr_y_xxxy_zzz[i] = 2.0 * tr_y_xy_zzz[i] * fe_0 + tr_y_xxy_zzz[i] * pa_x[i];
    }

    // Set up 170-180 components of targeted buffer : GF

    auto tr_y_xxxz_xxx = pbuffer.data(idx_dip_gf + 170);

    auto tr_y_xxxz_xxy = pbuffer.data(idx_dip_gf + 171);

    auto tr_y_xxxz_xxz = pbuffer.data(idx_dip_gf + 172);

    auto tr_y_xxxz_xyy = pbuffer.data(idx_dip_gf + 173);

    auto tr_y_xxxz_xyz = pbuffer.data(idx_dip_gf + 174);

    auto tr_y_xxxz_xzz = pbuffer.data(idx_dip_gf + 175);

    auto tr_y_xxxz_yyy = pbuffer.data(idx_dip_gf + 176);

    auto tr_y_xxxz_yyz = pbuffer.data(idx_dip_gf + 177);

    auto tr_y_xxxz_yzz = pbuffer.data(idx_dip_gf + 178);

    auto tr_y_xxxz_zzz = pbuffer.data(idx_dip_gf + 179);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xxx_xx,   \
                             tr_y_xxx_xxx,  \
                             tr_y_xxx_xxy,  \
                             tr_y_xxx_xxz,  \
                             tr_y_xxx_xy,   \
                             tr_y_xxx_xyy,  \
                             tr_y_xxx_xyz,  \
                             tr_y_xxx_xz,   \
                             tr_y_xxx_xzz,  \
                             tr_y_xxx_yyy,  \
                             tr_y_xxxz_xxx, \
                             tr_y_xxxz_xxy, \
                             tr_y_xxxz_xxz, \
                             tr_y_xxxz_xyy, \
                             tr_y_xxxz_xyz, \
                             tr_y_xxxz_xzz, \
                             tr_y_xxxz_yyy, \
                             tr_y_xxxz_yyz, \
                             tr_y_xxxz_yzz, \
                             tr_y_xxxz_zzz, \
                             tr_y_xxz_yyz,  \
                             tr_y_xxz_yzz,  \
                             tr_y_xxz_zzz,  \
                             tr_y_xz_yyz,   \
                             tr_y_xz_yzz,   \
                             tr_y_xz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_xxx[i] = tr_y_xxx_xxx[i] * pa_z[i];

        tr_y_xxxz_xxy[i] = tr_y_xxx_xxy[i] * pa_z[i];

        tr_y_xxxz_xxz[i] = tr_y_xxx_xx[i] * fe_0 + tr_y_xxx_xxz[i] * pa_z[i];

        tr_y_xxxz_xyy[i] = tr_y_xxx_xyy[i] * pa_z[i];

        tr_y_xxxz_xyz[i] = tr_y_xxx_xy[i] * fe_0 + tr_y_xxx_xyz[i] * pa_z[i];

        tr_y_xxxz_xzz[i] = 2.0 * tr_y_xxx_xz[i] * fe_0 + tr_y_xxx_xzz[i] * pa_z[i];

        tr_y_xxxz_yyy[i] = tr_y_xxx_yyy[i] * pa_z[i];

        tr_y_xxxz_yyz[i] = 2.0 * tr_y_xz_yyz[i] * fe_0 + tr_y_xxz_yyz[i] * pa_x[i];

        tr_y_xxxz_yzz[i] = 2.0 * tr_y_xz_yzz[i] * fe_0 + tr_y_xxz_yzz[i] * pa_x[i];

        tr_y_xxxz_zzz[i] = 2.0 * tr_y_xz_zzz[i] * fe_0 + tr_y_xxz_zzz[i] * pa_x[i];
    }

    // Set up 180-190 components of targeted buffer : GF

    auto tr_y_xxyy_xxx = pbuffer.data(idx_dip_gf + 180);

    auto tr_y_xxyy_xxy = pbuffer.data(idx_dip_gf + 181);

    auto tr_y_xxyy_xxz = pbuffer.data(idx_dip_gf + 182);

    auto tr_y_xxyy_xyy = pbuffer.data(idx_dip_gf + 183);

    auto tr_y_xxyy_xyz = pbuffer.data(idx_dip_gf + 184);

    auto tr_y_xxyy_xzz = pbuffer.data(idx_dip_gf + 185);

    auto tr_y_xxyy_yyy = pbuffer.data(idx_dip_gf + 186);

    auto tr_y_xxyy_yyz = pbuffer.data(idx_dip_gf + 187);

    auto tr_y_xxyy_yzz = pbuffer.data(idx_dip_gf + 188);

    auto tr_y_xxyy_zzz = pbuffer.data(idx_dip_gf + 189);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xxyy_xxx, \
                             tr_y_xxyy_xxy, \
                             tr_y_xxyy_xxz, \
                             tr_y_xxyy_xyy, \
                             tr_y_xxyy_xyz, \
                             tr_y_xxyy_xzz, \
                             tr_y_xxyy_yyy, \
                             tr_y_xxyy_yyz, \
                             tr_y_xxyy_yzz, \
                             tr_y_xxyy_zzz, \
                             tr_y_xyy_xx,   \
                             tr_y_xyy_xxx,  \
                             tr_y_xyy_xxy,  \
                             tr_y_xyy_xxz,  \
                             tr_y_xyy_xy,   \
                             tr_y_xyy_xyy,  \
                             tr_y_xyy_xyz,  \
                             tr_y_xyy_xz,   \
                             tr_y_xyy_xzz,  \
                             tr_y_xyy_yy,   \
                             tr_y_xyy_yyy,  \
                             tr_y_xyy_yyz,  \
                             tr_y_xyy_yz,   \
                             tr_y_xyy_yzz,  \
                             tr_y_xyy_zz,   \
                             tr_y_xyy_zzz,  \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xxz,   \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyz,   \
                             tr_y_yy_xzz,   \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yzz,   \
                             tr_y_yy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_xxx[i] = tr_y_yy_xxx[i] * fe_0 + 3.0 * tr_y_xyy_xx[i] * fe_0 + tr_y_xyy_xxx[i] * pa_x[i];

        tr_y_xxyy_xxy[i] = tr_y_yy_xxy[i] * fe_0 + 2.0 * tr_y_xyy_xy[i] * fe_0 + tr_y_xyy_xxy[i] * pa_x[i];

        tr_y_xxyy_xxz[i] = tr_y_yy_xxz[i] * fe_0 + 2.0 * tr_y_xyy_xz[i] * fe_0 + tr_y_xyy_xxz[i] * pa_x[i];

        tr_y_xxyy_xyy[i] = tr_y_yy_xyy[i] * fe_0 + tr_y_xyy_yy[i] * fe_0 + tr_y_xyy_xyy[i] * pa_x[i];

        tr_y_xxyy_xyz[i] = tr_y_yy_xyz[i] * fe_0 + tr_y_xyy_yz[i] * fe_0 + tr_y_xyy_xyz[i] * pa_x[i];

        tr_y_xxyy_xzz[i] = tr_y_yy_xzz[i] * fe_0 + tr_y_xyy_zz[i] * fe_0 + tr_y_xyy_xzz[i] * pa_x[i];

        tr_y_xxyy_yyy[i] = tr_y_yy_yyy[i] * fe_0 + tr_y_xyy_yyy[i] * pa_x[i];

        tr_y_xxyy_yyz[i] = tr_y_yy_yyz[i] * fe_0 + tr_y_xyy_yyz[i] * pa_x[i];

        tr_y_xxyy_yzz[i] = tr_y_yy_yzz[i] * fe_0 + tr_y_xyy_yzz[i] * pa_x[i];

        tr_y_xxyy_zzz[i] = tr_y_yy_zzz[i] * fe_0 + tr_y_xyy_zzz[i] * pa_x[i];
    }

    // Set up 190-200 components of targeted buffer : GF

    auto tr_y_xxyz_xxx = pbuffer.data(idx_dip_gf + 190);

    auto tr_y_xxyz_xxy = pbuffer.data(idx_dip_gf + 191);

    auto tr_y_xxyz_xxz = pbuffer.data(idx_dip_gf + 192);

    auto tr_y_xxyz_xyy = pbuffer.data(idx_dip_gf + 193);

    auto tr_y_xxyz_xyz = pbuffer.data(idx_dip_gf + 194);

    auto tr_y_xxyz_xzz = pbuffer.data(idx_dip_gf + 195);

    auto tr_y_xxyz_yyy = pbuffer.data(idx_dip_gf + 196);

    auto tr_y_xxyz_yyz = pbuffer.data(idx_dip_gf + 197);

    auto tr_y_xxyz_yzz = pbuffer.data(idx_dip_gf + 198);

    auto tr_y_xxyz_zzz = pbuffer.data(idx_dip_gf + 199);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             tr_y_xxy_xxx,  \
                             tr_y_xxy_xxy,  \
                             tr_y_xxy_xy,   \
                             tr_y_xxy_xyy,  \
                             tr_y_xxy_xyz,  \
                             tr_y_xxy_yyy,  \
                             tr_y_xxyz_xxx, \
                             tr_y_xxyz_xxy, \
                             tr_y_xxyz_xxz, \
                             tr_y_xxyz_xyy, \
                             tr_y_xxyz_xyz, \
                             tr_y_xxyz_xzz, \
                             tr_y_xxyz_yyy, \
                             tr_y_xxyz_yyz, \
                             tr_y_xxyz_yzz, \
                             tr_y_xxyz_zzz, \
                             tr_y_xxz_xxz,  \
                             tr_y_xxz_xzz,  \
                             tr_y_xyz_yyz,  \
                             tr_y_xyz_yzz,  \
                             tr_y_xyz_zzz,  \
                             tr_y_yz_yyz,   \
                             tr_y_yz_yzz,   \
                             tr_y_yz_zzz,   \
                             ts_xxz_xxz,    \
                             ts_xxz_xzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_xxx[i] = tr_y_xxy_xxx[i] * pa_z[i];

        tr_y_xxyz_xxy[i] = tr_y_xxy_xxy[i] * pa_z[i];

        tr_y_xxyz_xxz[i] = ts_xxz_xxz[i] * fe_0 + tr_y_xxz_xxz[i] * pa_y[i];

        tr_y_xxyz_xyy[i] = tr_y_xxy_xyy[i] * pa_z[i];

        tr_y_xxyz_xyz[i] = tr_y_xxy_xy[i] * fe_0 + tr_y_xxy_xyz[i] * pa_z[i];

        tr_y_xxyz_xzz[i] = ts_xxz_xzz[i] * fe_0 + tr_y_xxz_xzz[i] * pa_y[i];

        tr_y_xxyz_yyy[i] = tr_y_xxy_yyy[i] * pa_z[i];

        tr_y_xxyz_yyz[i] = tr_y_yz_yyz[i] * fe_0 + tr_y_xyz_yyz[i] * pa_x[i];

        tr_y_xxyz_yzz[i] = tr_y_yz_yzz[i] * fe_0 + tr_y_xyz_yzz[i] * pa_x[i];

        tr_y_xxyz_zzz[i] = tr_y_yz_zzz[i] * fe_0 + tr_y_xyz_zzz[i] * pa_x[i];
    }

    // Set up 200-210 components of targeted buffer : GF

    auto tr_y_xxzz_xxx = pbuffer.data(idx_dip_gf + 200);

    auto tr_y_xxzz_xxy = pbuffer.data(idx_dip_gf + 201);

    auto tr_y_xxzz_xxz = pbuffer.data(idx_dip_gf + 202);

    auto tr_y_xxzz_xyy = pbuffer.data(idx_dip_gf + 203);

    auto tr_y_xxzz_xyz = pbuffer.data(idx_dip_gf + 204);

    auto tr_y_xxzz_xzz = pbuffer.data(idx_dip_gf + 205);

    auto tr_y_xxzz_yyy = pbuffer.data(idx_dip_gf + 206);

    auto tr_y_xxzz_yyz = pbuffer.data(idx_dip_gf + 207);

    auto tr_y_xxzz_yzz = pbuffer.data(idx_dip_gf + 208);

    auto tr_y_xxzz_zzz = pbuffer.data(idx_dip_gf + 209);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xx_xxx,   \
                             tr_y_xx_xxy,   \
                             tr_y_xx_xyy,   \
                             tr_y_xxz_xxx,  \
                             tr_y_xxz_xxy,  \
                             tr_y_xxz_xyy,  \
                             tr_y_xxzz_xxx, \
                             tr_y_xxzz_xxy, \
                             tr_y_xxzz_xxz, \
                             tr_y_xxzz_xyy, \
                             tr_y_xxzz_xyz, \
                             tr_y_xxzz_xzz, \
                             tr_y_xxzz_yyy, \
                             tr_y_xxzz_yyz, \
                             tr_y_xxzz_yzz, \
                             tr_y_xxzz_zzz, \
                             tr_y_xzz_xxz,  \
                             tr_y_xzz_xyz,  \
                             tr_y_xzz_xz,   \
                             tr_y_xzz_xzz,  \
                             tr_y_xzz_yyy,  \
                             tr_y_xzz_yyz,  \
                             tr_y_xzz_yz,   \
                             tr_y_xzz_yzz,  \
                             tr_y_xzz_zz,   \
                             tr_y_xzz_zzz,  \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xyz,   \
                             tr_y_zz_xzz,   \
                             tr_y_zz_yyy,   \
                             tr_y_zz_yyz,   \
                             tr_y_zz_yzz,   \
                             tr_y_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_xxx[i] = tr_y_xx_xxx[i] * fe_0 + tr_y_xxz_xxx[i] * pa_z[i];

        tr_y_xxzz_xxy[i] = tr_y_xx_xxy[i] * fe_0 + tr_y_xxz_xxy[i] * pa_z[i];

        tr_y_xxzz_xxz[i] = tr_y_zz_xxz[i] * fe_0 + 2.0 * tr_y_xzz_xz[i] * fe_0 + tr_y_xzz_xxz[i] * pa_x[i];

        tr_y_xxzz_xyy[i] = tr_y_xx_xyy[i] * fe_0 + tr_y_xxz_xyy[i] * pa_z[i];

        tr_y_xxzz_xyz[i] = tr_y_zz_xyz[i] * fe_0 + tr_y_xzz_yz[i] * fe_0 + tr_y_xzz_xyz[i] * pa_x[i];

        tr_y_xxzz_xzz[i] = tr_y_zz_xzz[i] * fe_0 + tr_y_xzz_zz[i] * fe_0 + tr_y_xzz_xzz[i] * pa_x[i];

        tr_y_xxzz_yyy[i] = tr_y_zz_yyy[i] * fe_0 + tr_y_xzz_yyy[i] * pa_x[i];

        tr_y_xxzz_yyz[i] = tr_y_zz_yyz[i] * fe_0 + tr_y_xzz_yyz[i] * pa_x[i];

        tr_y_xxzz_yzz[i] = tr_y_zz_yzz[i] * fe_0 + tr_y_xzz_yzz[i] * pa_x[i];

        tr_y_xxzz_zzz[i] = tr_y_zz_zzz[i] * fe_0 + tr_y_xzz_zzz[i] * pa_x[i];
    }

    // Set up 210-220 components of targeted buffer : GF

    auto tr_y_xyyy_xxx = pbuffer.data(idx_dip_gf + 210);

    auto tr_y_xyyy_xxy = pbuffer.data(idx_dip_gf + 211);

    auto tr_y_xyyy_xxz = pbuffer.data(idx_dip_gf + 212);

    auto tr_y_xyyy_xyy = pbuffer.data(idx_dip_gf + 213);

    auto tr_y_xyyy_xyz = pbuffer.data(idx_dip_gf + 214);

    auto tr_y_xyyy_xzz = pbuffer.data(idx_dip_gf + 215);

    auto tr_y_xyyy_yyy = pbuffer.data(idx_dip_gf + 216);

    auto tr_y_xyyy_yyz = pbuffer.data(idx_dip_gf + 217);

    auto tr_y_xyyy_yzz = pbuffer.data(idx_dip_gf + 218);

    auto tr_y_xyyy_zzz = pbuffer.data(idx_dip_gf + 219);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xyyy_xxx, \
                             tr_y_xyyy_xxy, \
                             tr_y_xyyy_xxz, \
                             tr_y_xyyy_xyy, \
                             tr_y_xyyy_xyz, \
                             tr_y_xyyy_xzz, \
                             tr_y_xyyy_yyy, \
                             tr_y_xyyy_yyz, \
                             tr_y_xyyy_yzz, \
                             tr_y_xyyy_zzz, \
                             tr_y_yyy_xx,   \
                             tr_y_yyy_xxx,  \
                             tr_y_yyy_xxy,  \
                             tr_y_yyy_xxz,  \
                             tr_y_yyy_xy,   \
                             tr_y_yyy_xyy,  \
                             tr_y_yyy_xyz,  \
                             tr_y_yyy_xz,   \
                             tr_y_yyy_xzz,  \
                             tr_y_yyy_yy,   \
                             tr_y_yyy_yyy,  \
                             tr_y_yyy_yyz,  \
                             tr_y_yyy_yz,   \
                             tr_y_yyy_yzz,  \
                             tr_y_yyy_zz,   \
                             tr_y_yyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_xxx[i] = 3.0 * tr_y_yyy_xx[i] * fe_0 + tr_y_yyy_xxx[i] * pa_x[i];

        tr_y_xyyy_xxy[i] = 2.0 * tr_y_yyy_xy[i] * fe_0 + tr_y_yyy_xxy[i] * pa_x[i];

        tr_y_xyyy_xxz[i] = 2.0 * tr_y_yyy_xz[i] * fe_0 + tr_y_yyy_xxz[i] * pa_x[i];

        tr_y_xyyy_xyy[i] = tr_y_yyy_yy[i] * fe_0 + tr_y_yyy_xyy[i] * pa_x[i];

        tr_y_xyyy_xyz[i] = tr_y_yyy_yz[i] * fe_0 + tr_y_yyy_xyz[i] * pa_x[i];

        tr_y_xyyy_xzz[i] = tr_y_yyy_zz[i] * fe_0 + tr_y_yyy_xzz[i] * pa_x[i];

        tr_y_xyyy_yyy[i] = tr_y_yyy_yyy[i] * pa_x[i];

        tr_y_xyyy_yyz[i] = tr_y_yyy_yyz[i] * pa_x[i];

        tr_y_xyyy_yzz[i] = tr_y_yyy_yzz[i] * pa_x[i];

        tr_y_xyyy_zzz[i] = tr_y_yyy_zzz[i] * pa_x[i];
    }

    // Set up 220-230 components of targeted buffer : GF

    auto tr_y_xyyz_xxx = pbuffer.data(idx_dip_gf + 220);

    auto tr_y_xyyz_xxy = pbuffer.data(idx_dip_gf + 221);

    auto tr_y_xyyz_xxz = pbuffer.data(idx_dip_gf + 222);

    auto tr_y_xyyz_xyy = pbuffer.data(idx_dip_gf + 223);

    auto tr_y_xyyz_xyz = pbuffer.data(idx_dip_gf + 224);

    auto tr_y_xyyz_xzz = pbuffer.data(idx_dip_gf + 225);

    auto tr_y_xyyz_yyy = pbuffer.data(idx_dip_gf + 226);

    auto tr_y_xyyz_yyz = pbuffer.data(idx_dip_gf + 227);

    auto tr_y_xyyz_yzz = pbuffer.data(idx_dip_gf + 228);

    auto tr_y_xyyz_zzz = pbuffer.data(idx_dip_gf + 229);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_y_xyy_xxx,  \
                             tr_y_xyy_xxy,  \
                             tr_y_xyy_xyy,  \
                             tr_y_xyyz_xxx, \
                             tr_y_xyyz_xxy, \
                             tr_y_xyyz_xxz, \
                             tr_y_xyyz_xyy, \
                             tr_y_xyyz_xyz, \
                             tr_y_xyyz_xzz, \
                             tr_y_xyyz_yyy, \
                             tr_y_xyyz_yyz, \
                             tr_y_xyyz_yzz, \
                             tr_y_xyyz_zzz, \
                             tr_y_yyz_xxz,  \
                             tr_y_yyz_xyz,  \
                             tr_y_yyz_xz,   \
                             tr_y_yyz_xzz,  \
                             tr_y_yyz_yyy,  \
                             tr_y_yyz_yyz,  \
                             tr_y_yyz_yz,   \
                             tr_y_yyz_yzz,  \
                             tr_y_yyz_zz,   \
                             tr_y_yyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyz_xxx[i] = tr_y_xyy_xxx[i] * pa_z[i];

        tr_y_xyyz_xxy[i] = tr_y_xyy_xxy[i] * pa_z[i];

        tr_y_xyyz_xxz[i] = 2.0 * tr_y_yyz_xz[i] * fe_0 + tr_y_yyz_xxz[i] * pa_x[i];

        tr_y_xyyz_xyy[i] = tr_y_xyy_xyy[i] * pa_z[i];

        tr_y_xyyz_xyz[i] = tr_y_yyz_yz[i] * fe_0 + tr_y_yyz_xyz[i] * pa_x[i];

        tr_y_xyyz_xzz[i] = tr_y_yyz_zz[i] * fe_0 + tr_y_yyz_xzz[i] * pa_x[i];

        tr_y_xyyz_yyy[i] = tr_y_yyz_yyy[i] * pa_x[i];

        tr_y_xyyz_yyz[i] = tr_y_yyz_yyz[i] * pa_x[i];

        tr_y_xyyz_yzz[i] = tr_y_yyz_yzz[i] * pa_x[i];

        tr_y_xyyz_zzz[i] = tr_y_yyz_zzz[i] * pa_x[i];
    }

    // Set up 230-240 components of targeted buffer : GF

    auto tr_y_xyzz_xxx = pbuffer.data(idx_dip_gf + 230);

    auto tr_y_xyzz_xxy = pbuffer.data(idx_dip_gf + 231);

    auto tr_y_xyzz_xxz = pbuffer.data(idx_dip_gf + 232);

    auto tr_y_xyzz_xyy = pbuffer.data(idx_dip_gf + 233);

    auto tr_y_xyzz_xyz = pbuffer.data(idx_dip_gf + 234);

    auto tr_y_xyzz_xzz = pbuffer.data(idx_dip_gf + 235);

    auto tr_y_xyzz_yyy = pbuffer.data(idx_dip_gf + 236);

    auto tr_y_xyzz_yyz = pbuffer.data(idx_dip_gf + 237);

    auto tr_y_xyzz_yzz = pbuffer.data(idx_dip_gf + 238);

    auto tr_y_xyzz_zzz = pbuffer.data(idx_dip_gf + 239);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xyzz_xxx, \
                             tr_y_xyzz_xxy, \
                             tr_y_xyzz_xxz, \
                             tr_y_xyzz_xyy, \
                             tr_y_xyzz_xyz, \
                             tr_y_xyzz_xzz, \
                             tr_y_xyzz_yyy, \
                             tr_y_xyzz_yyz, \
                             tr_y_xyzz_yzz, \
                             tr_y_xyzz_zzz, \
                             tr_y_yzz_xx,   \
                             tr_y_yzz_xxx,  \
                             tr_y_yzz_xxy,  \
                             tr_y_yzz_xxz,  \
                             tr_y_yzz_xy,   \
                             tr_y_yzz_xyy,  \
                             tr_y_yzz_xyz,  \
                             tr_y_yzz_xz,   \
                             tr_y_yzz_xzz,  \
                             tr_y_yzz_yy,   \
                             tr_y_yzz_yyy,  \
                             tr_y_yzz_yyz,  \
                             tr_y_yzz_yz,   \
                             tr_y_yzz_yzz,  \
                             tr_y_yzz_zz,   \
                             tr_y_yzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_xxx[i] = 3.0 * tr_y_yzz_xx[i] * fe_0 + tr_y_yzz_xxx[i] * pa_x[i];

        tr_y_xyzz_xxy[i] = 2.0 * tr_y_yzz_xy[i] * fe_0 + tr_y_yzz_xxy[i] * pa_x[i];

        tr_y_xyzz_xxz[i] = 2.0 * tr_y_yzz_xz[i] * fe_0 + tr_y_yzz_xxz[i] * pa_x[i];

        tr_y_xyzz_xyy[i] = tr_y_yzz_yy[i] * fe_0 + tr_y_yzz_xyy[i] * pa_x[i];

        tr_y_xyzz_xyz[i] = tr_y_yzz_yz[i] * fe_0 + tr_y_yzz_xyz[i] * pa_x[i];

        tr_y_xyzz_xzz[i] = tr_y_yzz_zz[i] * fe_0 + tr_y_yzz_xzz[i] * pa_x[i];

        tr_y_xyzz_yyy[i] = tr_y_yzz_yyy[i] * pa_x[i];

        tr_y_xyzz_yyz[i] = tr_y_yzz_yyz[i] * pa_x[i];

        tr_y_xyzz_yzz[i] = tr_y_yzz_yzz[i] * pa_x[i];

        tr_y_xyzz_zzz[i] = tr_y_yzz_zzz[i] * pa_x[i];
    }

    // Set up 240-250 components of targeted buffer : GF

    auto tr_y_xzzz_xxx = pbuffer.data(idx_dip_gf + 240);

    auto tr_y_xzzz_xxy = pbuffer.data(idx_dip_gf + 241);

    auto tr_y_xzzz_xxz = pbuffer.data(idx_dip_gf + 242);

    auto tr_y_xzzz_xyy = pbuffer.data(idx_dip_gf + 243);

    auto tr_y_xzzz_xyz = pbuffer.data(idx_dip_gf + 244);

    auto tr_y_xzzz_xzz = pbuffer.data(idx_dip_gf + 245);

    auto tr_y_xzzz_yyy = pbuffer.data(idx_dip_gf + 246);

    auto tr_y_xzzz_yyz = pbuffer.data(idx_dip_gf + 247);

    auto tr_y_xzzz_yzz = pbuffer.data(idx_dip_gf + 248);

    auto tr_y_xzzz_zzz = pbuffer.data(idx_dip_gf + 249);

#pragma omp simd aligned(pa_x,              \
                             tr_y_xzzz_xxx, \
                             tr_y_xzzz_xxy, \
                             tr_y_xzzz_xxz, \
                             tr_y_xzzz_xyy, \
                             tr_y_xzzz_xyz, \
                             tr_y_xzzz_xzz, \
                             tr_y_xzzz_yyy, \
                             tr_y_xzzz_yyz, \
                             tr_y_xzzz_yzz, \
                             tr_y_xzzz_zzz, \
                             tr_y_zzz_xx,   \
                             tr_y_zzz_xxx,  \
                             tr_y_zzz_xxy,  \
                             tr_y_zzz_xxz,  \
                             tr_y_zzz_xy,   \
                             tr_y_zzz_xyy,  \
                             tr_y_zzz_xyz,  \
                             tr_y_zzz_xz,   \
                             tr_y_zzz_xzz,  \
                             tr_y_zzz_yy,   \
                             tr_y_zzz_yyy,  \
                             tr_y_zzz_yyz,  \
                             tr_y_zzz_yz,   \
                             tr_y_zzz_yzz,  \
                             tr_y_zzz_zz,   \
                             tr_y_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_xxx[i] = 3.0 * tr_y_zzz_xx[i] * fe_0 + tr_y_zzz_xxx[i] * pa_x[i];

        tr_y_xzzz_xxy[i] = 2.0 * tr_y_zzz_xy[i] * fe_0 + tr_y_zzz_xxy[i] * pa_x[i];

        tr_y_xzzz_xxz[i] = 2.0 * tr_y_zzz_xz[i] * fe_0 + tr_y_zzz_xxz[i] * pa_x[i];

        tr_y_xzzz_xyy[i] = tr_y_zzz_yy[i] * fe_0 + tr_y_zzz_xyy[i] * pa_x[i];

        tr_y_xzzz_xyz[i] = tr_y_zzz_yz[i] * fe_0 + tr_y_zzz_xyz[i] * pa_x[i];

        tr_y_xzzz_xzz[i] = tr_y_zzz_zz[i] * fe_0 + tr_y_zzz_xzz[i] * pa_x[i];

        tr_y_xzzz_yyy[i] = tr_y_zzz_yyy[i] * pa_x[i];

        tr_y_xzzz_yyz[i] = tr_y_zzz_yyz[i] * pa_x[i];

        tr_y_xzzz_yzz[i] = tr_y_zzz_yzz[i] * pa_x[i];

        tr_y_xzzz_zzz[i] = tr_y_zzz_zzz[i] * pa_x[i];
    }

    // Set up 250-260 components of targeted buffer : GF

    auto tr_y_yyyy_xxx = pbuffer.data(idx_dip_gf + 250);

    auto tr_y_yyyy_xxy = pbuffer.data(idx_dip_gf + 251);

    auto tr_y_yyyy_xxz = pbuffer.data(idx_dip_gf + 252);

    auto tr_y_yyyy_xyy = pbuffer.data(idx_dip_gf + 253);

    auto tr_y_yyyy_xyz = pbuffer.data(idx_dip_gf + 254);

    auto tr_y_yyyy_xzz = pbuffer.data(idx_dip_gf + 255);

    auto tr_y_yyyy_yyy = pbuffer.data(idx_dip_gf + 256);

    auto tr_y_yyyy_yyz = pbuffer.data(idx_dip_gf + 257);

    auto tr_y_yyyy_yzz = pbuffer.data(idx_dip_gf + 258);

    auto tr_y_yyyy_zzz = pbuffer.data(idx_dip_gf + 259);

#pragma omp simd aligned(pa_y,              \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xxz,   \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyz,   \
                             tr_y_yy_xzz,   \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yzz,   \
                             tr_y_yy_zzz,   \
                             tr_y_yyy_xx,   \
                             tr_y_yyy_xxx,  \
                             tr_y_yyy_xxy,  \
                             tr_y_yyy_xxz,  \
                             tr_y_yyy_xy,   \
                             tr_y_yyy_xyy,  \
                             tr_y_yyy_xyz,  \
                             tr_y_yyy_xz,   \
                             tr_y_yyy_xzz,  \
                             tr_y_yyy_yy,   \
                             tr_y_yyy_yyy,  \
                             tr_y_yyy_yyz,  \
                             tr_y_yyy_yz,   \
                             tr_y_yyy_yzz,  \
                             tr_y_yyy_zz,   \
                             tr_y_yyy_zzz,  \
                             tr_y_yyyy_xxx, \
                             tr_y_yyyy_xxy, \
                             tr_y_yyyy_xxz, \
                             tr_y_yyyy_xyy, \
                             tr_y_yyyy_xyz, \
                             tr_y_yyyy_xzz, \
                             tr_y_yyyy_yyy, \
                             tr_y_yyyy_yyz, \
                             tr_y_yyyy_yzz, \
                             tr_y_yyyy_zzz, \
                             ts_yyy_xxx,    \
                             ts_yyy_xxy,    \
                             ts_yyy_xxz,    \
                             ts_yyy_xyy,    \
                             ts_yyy_xyz,    \
                             ts_yyy_xzz,    \
                             ts_yyy_yyy,    \
                             ts_yyy_yyz,    \
                             ts_yyy_yzz,    \
                             ts_yyy_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_xxx[i] = 3.0 * tr_y_yy_xxx[i] * fe_0 + ts_yyy_xxx[i] * fe_0 + tr_y_yyy_xxx[i] * pa_y[i];

        tr_y_yyyy_xxy[i] = 3.0 * tr_y_yy_xxy[i] * fe_0 + tr_y_yyy_xx[i] * fe_0 + ts_yyy_xxy[i] * fe_0 + tr_y_yyy_xxy[i] * pa_y[i];

        tr_y_yyyy_xxz[i] = 3.0 * tr_y_yy_xxz[i] * fe_0 + ts_yyy_xxz[i] * fe_0 + tr_y_yyy_xxz[i] * pa_y[i];

        tr_y_yyyy_xyy[i] = 3.0 * tr_y_yy_xyy[i] * fe_0 + 2.0 * tr_y_yyy_xy[i] * fe_0 + ts_yyy_xyy[i] * fe_0 + tr_y_yyy_xyy[i] * pa_y[i];

        tr_y_yyyy_xyz[i] = 3.0 * tr_y_yy_xyz[i] * fe_0 + tr_y_yyy_xz[i] * fe_0 + ts_yyy_xyz[i] * fe_0 + tr_y_yyy_xyz[i] * pa_y[i];

        tr_y_yyyy_xzz[i] = 3.0 * tr_y_yy_xzz[i] * fe_0 + ts_yyy_xzz[i] * fe_0 + tr_y_yyy_xzz[i] * pa_y[i];

        tr_y_yyyy_yyy[i] = 3.0 * tr_y_yy_yyy[i] * fe_0 + 3.0 * tr_y_yyy_yy[i] * fe_0 + ts_yyy_yyy[i] * fe_0 + tr_y_yyy_yyy[i] * pa_y[i];

        tr_y_yyyy_yyz[i] = 3.0 * tr_y_yy_yyz[i] * fe_0 + 2.0 * tr_y_yyy_yz[i] * fe_0 + ts_yyy_yyz[i] * fe_0 + tr_y_yyy_yyz[i] * pa_y[i];

        tr_y_yyyy_yzz[i] = 3.0 * tr_y_yy_yzz[i] * fe_0 + tr_y_yyy_zz[i] * fe_0 + ts_yyy_yzz[i] * fe_0 + tr_y_yyy_yzz[i] * pa_y[i];

        tr_y_yyyy_zzz[i] = 3.0 * tr_y_yy_zzz[i] * fe_0 + ts_yyy_zzz[i] * fe_0 + tr_y_yyy_zzz[i] * pa_y[i];
    }

    // Set up 260-270 components of targeted buffer : GF

    auto tr_y_yyyz_xxx = pbuffer.data(idx_dip_gf + 260);

    auto tr_y_yyyz_xxy = pbuffer.data(idx_dip_gf + 261);

    auto tr_y_yyyz_xxz = pbuffer.data(idx_dip_gf + 262);

    auto tr_y_yyyz_xyy = pbuffer.data(idx_dip_gf + 263);

    auto tr_y_yyyz_xyz = pbuffer.data(idx_dip_gf + 264);

    auto tr_y_yyyz_xzz = pbuffer.data(idx_dip_gf + 265);

    auto tr_y_yyyz_yyy = pbuffer.data(idx_dip_gf + 266);

    auto tr_y_yyyz_yyz = pbuffer.data(idx_dip_gf + 267);

    auto tr_y_yyyz_yzz = pbuffer.data(idx_dip_gf + 268);

    auto tr_y_yyyz_zzz = pbuffer.data(idx_dip_gf + 269);

#pragma omp simd aligned(pa_z,              \
                             tr_y_yyy_xx,   \
                             tr_y_yyy_xxx,  \
                             tr_y_yyy_xxy,  \
                             tr_y_yyy_xxz,  \
                             tr_y_yyy_xy,   \
                             tr_y_yyy_xyy,  \
                             tr_y_yyy_xyz,  \
                             tr_y_yyy_xz,   \
                             tr_y_yyy_xzz,  \
                             tr_y_yyy_yy,   \
                             tr_y_yyy_yyy,  \
                             tr_y_yyy_yyz,  \
                             tr_y_yyy_yz,   \
                             tr_y_yyy_yzz,  \
                             tr_y_yyy_zz,   \
                             tr_y_yyy_zzz,  \
                             tr_y_yyyz_xxx, \
                             tr_y_yyyz_xxy, \
                             tr_y_yyyz_xxz, \
                             tr_y_yyyz_xyy, \
                             tr_y_yyyz_xyz, \
                             tr_y_yyyz_xzz, \
                             tr_y_yyyz_yyy, \
                             tr_y_yyyz_yyz, \
                             tr_y_yyyz_yzz, \
                             tr_y_yyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_xxx[i] = tr_y_yyy_xxx[i] * pa_z[i];

        tr_y_yyyz_xxy[i] = tr_y_yyy_xxy[i] * pa_z[i];

        tr_y_yyyz_xxz[i] = tr_y_yyy_xx[i] * fe_0 + tr_y_yyy_xxz[i] * pa_z[i];

        tr_y_yyyz_xyy[i] = tr_y_yyy_xyy[i] * pa_z[i];

        tr_y_yyyz_xyz[i] = tr_y_yyy_xy[i] * fe_0 + tr_y_yyy_xyz[i] * pa_z[i];

        tr_y_yyyz_xzz[i] = 2.0 * tr_y_yyy_xz[i] * fe_0 + tr_y_yyy_xzz[i] * pa_z[i];

        tr_y_yyyz_yyy[i] = tr_y_yyy_yyy[i] * pa_z[i];

        tr_y_yyyz_yyz[i] = tr_y_yyy_yy[i] * fe_0 + tr_y_yyy_yyz[i] * pa_z[i];

        tr_y_yyyz_yzz[i] = 2.0 * tr_y_yyy_yz[i] * fe_0 + tr_y_yyy_yzz[i] * pa_z[i];

        tr_y_yyyz_zzz[i] = 3.0 * tr_y_yyy_zz[i] * fe_0 + tr_y_yyy_zzz[i] * pa_z[i];
    }

    // Set up 270-280 components of targeted buffer : GF

    auto tr_y_yyzz_xxx = pbuffer.data(idx_dip_gf + 270);

    auto tr_y_yyzz_xxy = pbuffer.data(idx_dip_gf + 271);

    auto tr_y_yyzz_xxz = pbuffer.data(idx_dip_gf + 272);

    auto tr_y_yyzz_xyy = pbuffer.data(idx_dip_gf + 273);

    auto tr_y_yyzz_xyz = pbuffer.data(idx_dip_gf + 274);

    auto tr_y_yyzz_xzz = pbuffer.data(idx_dip_gf + 275);

    auto tr_y_yyzz_yyy = pbuffer.data(idx_dip_gf + 276);

    auto tr_y_yyzz_yyz = pbuffer.data(idx_dip_gf + 277);

    auto tr_y_yyzz_yzz = pbuffer.data(idx_dip_gf + 278);

    auto tr_y_yyzz_zzz = pbuffer.data(idx_dip_gf + 279);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yy_xxx,   \
                             tr_y_yy_xxy,   \
                             tr_y_yy_xyy,   \
                             tr_y_yy_xyz,   \
                             tr_y_yy_yyy,   \
                             tr_y_yy_yyz,   \
                             tr_y_yy_yzz,   \
                             tr_y_yyz_xxx,  \
                             tr_y_yyz_xxy,  \
                             tr_y_yyz_xy,   \
                             tr_y_yyz_xyy,  \
                             tr_y_yyz_xyz,  \
                             tr_y_yyz_yy,   \
                             tr_y_yyz_yyy,  \
                             tr_y_yyz_yyz,  \
                             tr_y_yyz_yz,   \
                             tr_y_yyz_yzz,  \
                             tr_y_yyzz_xxx, \
                             tr_y_yyzz_xxy, \
                             tr_y_yyzz_xxz, \
                             tr_y_yyzz_xyy, \
                             tr_y_yyzz_xyz, \
                             tr_y_yyzz_xzz, \
                             tr_y_yyzz_yyy, \
                             tr_y_yyzz_yyz, \
                             tr_y_yyzz_yzz, \
                             tr_y_yyzz_zzz, \
                             tr_y_yzz_xxz,  \
                             tr_y_yzz_xzz,  \
                             tr_y_yzz_zzz,  \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xzz,   \
                             tr_y_zz_zzz,   \
                             ts_yzz_xxz,    \
                             ts_yzz_xzz,    \
                             ts_yzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_xxx[i] = tr_y_yy_xxx[i] * fe_0 + tr_y_yyz_xxx[i] * pa_z[i];

        tr_y_yyzz_xxy[i] = tr_y_yy_xxy[i] * fe_0 + tr_y_yyz_xxy[i] * pa_z[i];

        tr_y_yyzz_xxz[i] = tr_y_zz_xxz[i] * fe_0 + ts_yzz_xxz[i] * fe_0 + tr_y_yzz_xxz[i] * pa_y[i];

        tr_y_yyzz_xyy[i] = tr_y_yy_xyy[i] * fe_0 + tr_y_yyz_xyy[i] * pa_z[i];

        tr_y_yyzz_xyz[i] = tr_y_yy_xyz[i] * fe_0 + tr_y_yyz_xy[i] * fe_0 + tr_y_yyz_xyz[i] * pa_z[i];

        tr_y_yyzz_xzz[i] = tr_y_zz_xzz[i] * fe_0 + ts_yzz_xzz[i] * fe_0 + tr_y_yzz_xzz[i] * pa_y[i];

        tr_y_yyzz_yyy[i] = tr_y_yy_yyy[i] * fe_0 + tr_y_yyz_yyy[i] * pa_z[i];

        tr_y_yyzz_yyz[i] = tr_y_yy_yyz[i] * fe_0 + tr_y_yyz_yy[i] * fe_0 + tr_y_yyz_yyz[i] * pa_z[i];

        tr_y_yyzz_yzz[i] = tr_y_yy_yzz[i] * fe_0 + 2.0 * tr_y_yyz_yz[i] * fe_0 + tr_y_yyz_yzz[i] * pa_z[i];

        tr_y_yyzz_zzz[i] = tr_y_zz_zzz[i] * fe_0 + ts_yzz_zzz[i] * fe_0 + tr_y_yzz_zzz[i] * pa_y[i];
    }

    // Set up 280-290 components of targeted buffer : GF

    auto tr_y_yzzz_xxx = pbuffer.data(idx_dip_gf + 280);

    auto tr_y_yzzz_xxy = pbuffer.data(idx_dip_gf + 281);

    auto tr_y_yzzz_xxz = pbuffer.data(idx_dip_gf + 282);

    auto tr_y_yzzz_xyy = pbuffer.data(idx_dip_gf + 283);

    auto tr_y_yzzz_xyz = pbuffer.data(idx_dip_gf + 284);

    auto tr_y_yzzz_xzz = pbuffer.data(idx_dip_gf + 285);

    auto tr_y_yzzz_yyy = pbuffer.data(idx_dip_gf + 286);

    auto tr_y_yzzz_yyz = pbuffer.data(idx_dip_gf + 287);

    auto tr_y_yzzz_yzz = pbuffer.data(idx_dip_gf + 288);

    auto tr_y_yzzz_zzz = pbuffer.data(idx_dip_gf + 289);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_y_yz_xxy,   \
                             tr_y_yz_xyy,   \
                             tr_y_yz_yyy,   \
                             tr_y_yzz_xxy,  \
                             tr_y_yzz_xyy,  \
                             tr_y_yzz_yyy,  \
                             tr_y_yzzz_xxx, \
                             tr_y_yzzz_xxy, \
                             tr_y_yzzz_xxz, \
                             tr_y_yzzz_xyy, \
                             tr_y_yzzz_xyz, \
                             tr_y_yzzz_xzz, \
                             tr_y_yzzz_yyy, \
                             tr_y_yzzz_yyz, \
                             tr_y_yzzz_yzz, \
                             tr_y_yzzz_zzz, \
                             tr_y_zzz_xxx,  \
                             tr_y_zzz_xxz,  \
                             tr_y_zzz_xyz,  \
                             tr_y_zzz_xz,   \
                             tr_y_zzz_xzz,  \
                             tr_y_zzz_yyz,  \
                             tr_y_zzz_yz,   \
                             tr_y_zzz_yzz,  \
                             tr_y_zzz_zz,   \
                             tr_y_zzz_zzz,  \
                             ts_zzz_xxx,    \
                             ts_zzz_xxz,    \
                             ts_zzz_xyz,    \
                             ts_zzz_xzz,    \
                             ts_zzz_yyz,    \
                             ts_zzz_yzz,    \
                             ts_zzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_xxx[i] = ts_zzz_xxx[i] * fe_0 + tr_y_zzz_xxx[i] * pa_y[i];

        tr_y_yzzz_xxy[i] = 2.0 * tr_y_yz_xxy[i] * fe_0 + tr_y_yzz_xxy[i] * pa_z[i];

        tr_y_yzzz_xxz[i] = ts_zzz_xxz[i] * fe_0 + tr_y_zzz_xxz[i] * pa_y[i];

        tr_y_yzzz_xyy[i] = 2.0 * tr_y_yz_xyy[i] * fe_0 + tr_y_yzz_xyy[i] * pa_z[i];

        tr_y_yzzz_xyz[i] = tr_y_zzz_xz[i] * fe_0 + ts_zzz_xyz[i] * fe_0 + tr_y_zzz_xyz[i] * pa_y[i];

        tr_y_yzzz_xzz[i] = ts_zzz_xzz[i] * fe_0 + tr_y_zzz_xzz[i] * pa_y[i];

        tr_y_yzzz_yyy[i] = 2.0 * tr_y_yz_yyy[i] * fe_0 + tr_y_yzz_yyy[i] * pa_z[i];

        tr_y_yzzz_yyz[i] = 2.0 * tr_y_zzz_yz[i] * fe_0 + ts_zzz_yyz[i] * fe_0 + tr_y_zzz_yyz[i] * pa_y[i];

        tr_y_yzzz_yzz[i] = tr_y_zzz_zz[i] * fe_0 + ts_zzz_yzz[i] * fe_0 + tr_y_zzz_yzz[i] * pa_y[i];

        tr_y_yzzz_zzz[i] = ts_zzz_zzz[i] * fe_0 + tr_y_zzz_zzz[i] * pa_y[i];
    }

    // Set up 290-300 components of targeted buffer : GF

    auto tr_y_zzzz_xxx = pbuffer.data(idx_dip_gf + 290);

    auto tr_y_zzzz_xxy = pbuffer.data(idx_dip_gf + 291);

    auto tr_y_zzzz_xxz = pbuffer.data(idx_dip_gf + 292);

    auto tr_y_zzzz_xyy = pbuffer.data(idx_dip_gf + 293);

    auto tr_y_zzzz_xyz = pbuffer.data(idx_dip_gf + 294);

    auto tr_y_zzzz_xzz = pbuffer.data(idx_dip_gf + 295);

    auto tr_y_zzzz_yyy = pbuffer.data(idx_dip_gf + 296);

    auto tr_y_zzzz_yyz = pbuffer.data(idx_dip_gf + 297);

    auto tr_y_zzzz_yzz = pbuffer.data(idx_dip_gf + 298);

    auto tr_y_zzzz_zzz = pbuffer.data(idx_dip_gf + 299);

#pragma omp simd aligned(pa_z,              \
                             tr_y_zz_xxx,   \
                             tr_y_zz_xxy,   \
                             tr_y_zz_xxz,   \
                             tr_y_zz_xyy,   \
                             tr_y_zz_xyz,   \
                             tr_y_zz_xzz,   \
                             tr_y_zz_yyy,   \
                             tr_y_zz_yyz,   \
                             tr_y_zz_yzz,   \
                             tr_y_zz_zzz,   \
                             tr_y_zzz_xx,   \
                             tr_y_zzz_xxx,  \
                             tr_y_zzz_xxy,  \
                             tr_y_zzz_xxz,  \
                             tr_y_zzz_xy,   \
                             tr_y_zzz_xyy,  \
                             tr_y_zzz_xyz,  \
                             tr_y_zzz_xz,   \
                             tr_y_zzz_xzz,  \
                             tr_y_zzz_yy,   \
                             tr_y_zzz_yyy,  \
                             tr_y_zzz_yyz,  \
                             tr_y_zzz_yz,   \
                             tr_y_zzz_yzz,  \
                             tr_y_zzz_zz,   \
                             tr_y_zzz_zzz,  \
                             tr_y_zzzz_xxx, \
                             tr_y_zzzz_xxy, \
                             tr_y_zzzz_xxz, \
                             tr_y_zzzz_xyy, \
                             tr_y_zzzz_xyz, \
                             tr_y_zzzz_xzz, \
                             tr_y_zzzz_yyy, \
                             tr_y_zzzz_yyz, \
                             tr_y_zzzz_yzz, \
                             tr_y_zzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_xxx[i] = 3.0 * tr_y_zz_xxx[i] * fe_0 + tr_y_zzz_xxx[i] * pa_z[i];

        tr_y_zzzz_xxy[i] = 3.0 * tr_y_zz_xxy[i] * fe_0 + tr_y_zzz_xxy[i] * pa_z[i];

        tr_y_zzzz_xxz[i] = 3.0 * tr_y_zz_xxz[i] * fe_0 + tr_y_zzz_xx[i] * fe_0 + tr_y_zzz_xxz[i] * pa_z[i];

        tr_y_zzzz_xyy[i] = 3.0 * tr_y_zz_xyy[i] * fe_0 + tr_y_zzz_xyy[i] * pa_z[i];

        tr_y_zzzz_xyz[i] = 3.0 * tr_y_zz_xyz[i] * fe_0 + tr_y_zzz_xy[i] * fe_0 + tr_y_zzz_xyz[i] * pa_z[i];

        tr_y_zzzz_xzz[i] = 3.0 * tr_y_zz_xzz[i] * fe_0 + 2.0 * tr_y_zzz_xz[i] * fe_0 + tr_y_zzz_xzz[i] * pa_z[i];

        tr_y_zzzz_yyy[i] = 3.0 * tr_y_zz_yyy[i] * fe_0 + tr_y_zzz_yyy[i] * pa_z[i];

        tr_y_zzzz_yyz[i] = 3.0 * tr_y_zz_yyz[i] * fe_0 + tr_y_zzz_yy[i] * fe_0 + tr_y_zzz_yyz[i] * pa_z[i];

        tr_y_zzzz_yzz[i] = 3.0 * tr_y_zz_yzz[i] * fe_0 + 2.0 * tr_y_zzz_yz[i] * fe_0 + tr_y_zzz_yzz[i] * pa_z[i];

        tr_y_zzzz_zzz[i] = 3.0 * tr_y_zz_zzz[i] * fe_0 + 3.0 * tr_y_zzz_zz[i] * fe_0 + tr_y_zzz_zzz[i] * pa_z[i];
    }

    // Set up 300-310 components of targeted buffer : GF

    auto tr_z_xxxx_xxx = pbuffer.data(idx_dip_gf + 300);

    auto tr_z_xxxx_xxy = pbuffer.data(idx_dip_gf + 301);

    auto tr_z_xxxx_xxz = pbuffer.data(idx_dip_gf + 302);

    auto tr_z_xxxx_xyy = pbuffer.data(idx_dip_gf + 303);

    auto tr_z_xxxx_xyz = pbuffer.data(idx_dip_gf + 304);

    auto tr_z_xxxx_xzz = pbuffer.data(idx_dip_gf + 305);

    auto tr_z_xxxx_yyy = pbuffer.data(idx_dip_gf + 306);

    auto tr_z_xxxx_yyz = pbuffer.data(idx_dip_gf + 307);

    auto tr_z_xxxx_yzz = pbuffer.data(idx_dip_gf + 308);

    auto tr_z_xxxx_zzz = pbuffer.data(idx_dip_gf + 309);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xx_xxx,   \
                             tr_z_xx_xxy,   \
                             tr_z_xx_xxz,   \
                             tr_z_xx_xyy,   \
                             tr_z_xx_xyz,   \
                             tr_z_xx_xzz,   \
                             tr_z_xx_yyy,   \
                             tr_z_xx_yyz,   \
                             tr_z_xx_yzz,   \
                             tr_z_xx_zzz,   \
                             tr_z_xxx_xx,   \
                             tr_z_xxx_xxx,  \
                             tr_z_xxx_xxy,  \
                             tr_z_xxx_xxz,  \
                             tr_z_xxx_xy,   \
                             tr_z_xxx_xyy,  \
                             tr_z_xxx_xyz,  \
                             tr_z_xxx_xz,   \
                             tr_z_xxx_xzz,  \
                             tr_z_xxx_yy,   \
                             tr_z_xxx_yyy,  \
                             tr_z_xxx_yyz,  \
                             tr_z_xxx_yz,   \
                             tr_z_xxx_yzz,  \
                             tr_z_xxx_zz,   \
                             tr_z_xxx_zzz,  \
                             tr_z_xxxx_xxx, \
                             tr_z_xxxx_xxy, \
                             tr_z_xxxx_xxz, \
                             tr_z_xxxx_xyy, \
                             tr_z_xxxx_xyz, \
                             tr_z_xxxx_xzz, \
                             tr_z_xxxx_yyy, \
                             tr_z_xxxx_yyz, \
                             tr_z_xxxx_yzz, \
                             tr_z_xxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_xxx[i] = 3.0 * tr_z_xx_xxx[i] * fe_0 + 3.0 * tr_z_xxx_xx[i] * fe_0 + tr_z_xxx_xxx[i] * pa_x[i];

        tr_z_xxxx_xxy[i] = 3.0 * tr_z_xx_xxy[i] * fe_0 + 2.0 * tr_z_xxx_xy[i] * fe_0 + tr_z_xxx_xxy[i] * pa_x[i];

        tr_z_xxxx_xxz[i] = 3.0 * tr_z_xx_xxz[i] * fe_0 + 2.0 * tr_z_xxx_xz[i] * fe_0 + tr_z_xxx_xxz[i] * pa_x[i];

        tr_z_xxxx_xyy[i] = 3.0 * tr_z_xx_xyy[i] * fe_0 + tr_z_xxx_yy[i] * fe_0 + tr_z_xxx_xyy[i] * pa_x[i];

        tr_z_xxxx_xyz[i] = 3.0 * tr_z_xx_xyz[i] * fe_0 + tr_z_xxx_yz[i] * fe_0 + tr_z_xxx_xyz[i] * pa_x[i];

        tr_z_xxxx_xzz[i] = 3.0 * tr_z_xx_xzz[i] * fe_0 + tr_z_xxx_zz[i] * fe_0 + tr_z_xxx_xzz[i] * pa_x[i];

        tr_z_xxxx_yyy[i] = 3.0 * tr_z_xx_yyy[i] * fe_0 + tr_z_xxx_yyy[i] * pa_x[i];

        tr_z_xxxx_yyz[i] = 3.0 * tr_z_xx_yyz[i] * fe_0 + tr_z_xxx_yyz[i] * pa_x[i];

        tr_z_xxxx_yzz[i] = 3.0 * tr_z_xx_yzz[i] * fe_0 + tr_z_xxx_yzz[i] * pa_x[i];

        tr_z_xxxx_zzz[i] = 3.0 * tr_z_xx_zzz[i] * fe_0 + tr_z_xxx_zzz[i] * pa_x[i];
    }

    // Set up 310-320 components of targeted buffer : GF

    auto tr_z_xxxy_xxx = pbuffer.data(idx_dip_gf + 310);

    auto tr_z_xxxy_xxy = pbuffer.data(idx_dip_gf + 311);

    auto tr_z_xxxy_xxz = pbuffer.data(idx_dip_gf + 312);

    auto tr_z_xxxy_xyy = pbuffer.data(idx_dip_gf + 313);

    auto tr_z_xxxy_xyz = pbuffer.data(idx_dip_gf + 314);

    auto tr_z_xxxy_xzz = pbuffer.data(idx_dip_gf + 315);

    auto tr_z_xxxy_yyy = pbuffer.data(idx_dip_gf + 316);

    auto tr_z_xxxy_yyz = pbuffer.data(idx_dip_gf + 317);

    auto tr_z_xxxy_yzz = pbuffer.data(idx_dip_gf + 318);

    auto tr_z_xxxy_zzz = pbuffer.data(idx_dip_gf + 319);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxx_xx,   \
                             tr_z_xxx_xxx,  \
                             tr_z_xxx_xxy,  \
                             tr_z_xxx_xxz,  \
                             tr_z_xxx_xy,   \
                             tr_z_xxx_xyy,  \
                             tr_z_xxx_xyz,  \
                             tr_z_xxx_xz,   \
                             tr_z_xxx_xzz,  \
                             tr_z_xxx_zzz,  \
                             tr_z_xxxy_xxx, \
                             tr_z_xxxy_xxy, \
                             tr_z_xxxy_xxz, \
                             tr_z_xxxy_xyy, \
                             tr_z_xxxy_xyz, \
                             tr_z_xxxy_xzz, \
                             tr_z_xxxy_yyy, \
                             tr_z_xxxy_yyz, \
                             tr_z_xxxy_yzz, \
                             tr_z_xxxy_zzz, \
                             tr_z_xxy_yyy,  \
                             tr_z_xxy_yyz,  \
                             tr_z_xxy_yzz,  \
                             tr_z_xy_yyy,   \
                             tr_z_xy_yyz,   \
                             tr_z_xy_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_xxx[i] = tr_z_xxx_xxx[i] * pa_y[i];

        tr_z_xxxy_xxy[i] = tr_z_xxx_xx[i] * fe_0 + tr_z_xxx_xxy[i] * pa_y[i];

        tr_z_xxxy_xxz[i] = tr_z_xxx_xxz[i] * pa_y[i];

        tr_z_xxxy_xyy[i] = 2.0 * tr_z_xxx_xy[i] * fe_0 + tr_z_xxx_xyy[i] * pa_y[i];

        tr_z_xxxy_xyz[i] = tr_z_xxx_xz[i] * fe_0 + tr_z_xxx_xyz[i] * pa_y[i];

        tr_z_xxxy_xzz[i] = tr_z_xxx_xzz[i] * pa_y[i];

        tr_z_xxxy_yyy[i] = 2.0 * tr_z_xy_yyy[i] * fe_0 + tr_z_xxy_yyy[i] * pa_x[i];

        tr_z_xxxy_yyz[i] = 2.0 * tr_z_xy_yyz[i] * fe_0 + tr_z_xxy_yyz[i] * pa_x[i];

        tr_z_xxxy_yzz[i] = 2.0 * tr_z_xy_yzz[i] * fe_0 + tr_z_xxy_yzz[i] * pa_x[i];

        tr_z_xxxy_zzz[i] = tr_z_xxx_zzz[i] * pa_y[i];
    }

    // Set up 320-330 components of targeted buffer : GF

    auto tr_z_xxxz_xxx = pbuffer.data(idx_dip_gf + 320);

    auto tr_z_xxxz_xxy = pbuffer.data(idx_dip_gf + 321);

    auto tr_z_xxxz_xxz = pbuffer.data(idx_dip_gf + 322);

    auto tr_z_xxxz_xyy = pbuffer.data(idx_dip_gf + 323);

    auto tr_z_xxxz_xyz = pbuffer.data(idx_dip_gf + 324);

    auto tr_z_xxxz_xzz = pbuffer.data(idx_dip_gf + 325);

    auto tr_z_xxxz_yyy = pbuffer.data(idx_dip_gf + 326);

    auto tr_z_xxxz_yyz = pbuffer.data(idx_dip_gf + 327);

    auto tr_z_xxxz_yzz = pbuffer.data(idx_dip_gf + 328);

    auto tr_z_xxxz_zzz = pbuffer.data(idx_dip_gf + 329);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             tr_z_xxx_xxx,  \
                             tr_z_xxx_xxy,  \
                             tr_z_xxx_xyy,  \
                             tr_z_xxxz_xxx, \
                             tr_z_xxxz_xxy, \
                             tr_z_xxxz_xxz, \
                             tr_z_xxxz_xyy, \
                             tr_z_xxxz_xyz, \
                             tr_z_xxxz_xzz, \
                             tr_z_xxxz_yyy, \
                             tr_z_xxxz_yyz, \
                             tr_z_xxxz_yzz, \
                             tr_z_xxxz_zzz, \
                             tr_z_xxz_xxz,  \
                             tr_z_xxz_xyz,  \
                             tr_z_xxz_xz,   \
                             tr_z_xxz_xzz,  \
                             tr_z_xxz_yyy,  \
                             tr_z_xxz_yyz,  \
                             tr_z_xxz_yz,   \
                             tr_z_xxz_yzz,  \
                             tr_z_xxz_zz,   \
                             tr_z_xxz_zzz,  \
                             tr_z_xz_xxz,   \
                             tr_z_xz_xyz,   \
                             tr_z_xz_xzz,   \
                             tr_z_xz_yyy,   \
                             tr_z_xz_yyz,   \
                             tr_z_xz_yzz,   \
                             tr_z_xz_zzz,   \
                             ts_xxx_xxx,    \
                             ts_xxx_xxy,    \
                             ts_xxx_xyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_xxx[i] = ts_xxx_xxx[i] * fe_0 + tr_z_xxx_xxx[i] * pa_z[i];

        tr_z_xxxz_xxy[i] = ts_xxx_xxy[i] * fe_0 + tr_z_xxx_xxy[i] * pa_z[i];

        tr_z_xxxz_xxz[i] = 2.0 * tr_z_xz_xxz[i] * fe_0 + 2.0 * tr_z_xxz_xz[i] * fe_0 + tr_z_xxz_xxz[i] * pa_x[i];

        tr_z_xxxz_xyy[i] = ts_xxx_xyy[i] * fe_0 + tr_z_xxx_xyy[i] * pa_z[i];

        tr_z_xxxz_xyz[i] = 2.0 * tr_z_xz_xyz[i] * fe_0 + tr_z_xxz_yz[i] * fe_0 + tr_z_xxz_xyz[i] * pa_x[i];

        tr_z_xxxz_xzz[i] = 2.0 * tr_z_xz_xzz[i] * fe_0 + tr_z_xxz_zz[i] * fe_0 + tr_z_xxz_xzz[i] * pa_x[i];

        tr_z_xxxz_yyy[i] = 2.0 * tr_z_xz_yyy[i] * fe_0 + tr_z_xxz_yyy[i] * pa_x[i];

        tr_z_xxxz_yyz[i] = 2.0 * tr_z_xz_yyz[i] * fe_0 + tr_z_xxz_yyz[i] * pa_x[i];

        tr_z_xxxz_yzz[i] = 2.0 * tr_z_xz_yzz[i] * fe_0 + tr_z_xxz_yzz[i] * pa_x[i];

        tr_z_xxxz_zzz[i] = 2.0 * tr_z_xz_zzz[i] * fe_0 + tr_z_xxz_zzz[i] * pa_x[i];
    }

    // Set up 330-340 components of targeted buffer : GF

    auto tr_z_xxyy_xxx = pbuffer.data(idx_dip_gf + 330);

    auto tr_z_xxyy_xxy = pbuffer.data(idx_dip_gf + 331);

    auto tr_z_xxyy_xxz = pbuffer.data(idx_dip_gf + 332);

    auto tr_z_xxyy_xyy = pbuffer.data(idx_dip_gf + 333);

    auto tr_z_xxyy_xyz = pbuffer.data(idx_dip_gf + 334);

    auto tr_z_xxyy_xzz = pbuffer.data(idx_dip_gf + 335);

    auto tr_z_xxyy_yyy = pbuffer.data(idx_dip_gf + 336);

    auto tr_z_xxyy_yyz = pbuffer.data(idx_dip_gf + 337);

    auto tr_z_xxyy_yzz = pbuffer.data(idx_dip_gf + 338);

    auto tr_z_xxyy_zzz = pbuffer.data(idx_dip_gf + 339);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xx_xxx,   \
                             tr_z_xx_xxz,   \
                             tr_z_xx_xzz,   \
                             tr_z_xxy_xxx,  \
                             tr_z_xxy_xxz,  \
                             tr_z_xxy_xzz,  \
                             tr_z_xxyy_xxx, \
                             tr_z_xxyy_xxy, \
                             tr_z_xxyy_xxz, \
                             tr_z_xxyy_xyy, \
                             tr_z_xxyy_xyz, \
                             tr_z_xxyy_xzz, \
                             tr_z_xxyy_yyy, \
                             tr_z_xxyy_yyz, \
                             tr_z_xxyy_yzz, \
                             tr_z_xxyy_zzz, \
                             tr_z_xyy_xxy,  \
                             tr_z_xyy_xy,   \
                             tr_z_xyy_xyy,  \
                             tr_z_xyy_xyz,  \
                             tr_z_xyy_yy,   \
                             tr_z_xyy_yyy,  \
                             tr_z_xyy_yyz,  \
                             tr_z_xyy_yz,   \
                             tr_z_xyy_yzz,  \
                             tr_z_xyy_zzz,  \
                             tr_z_yy_xxy,   \
                             tr_z_yy_xyy,   \
                             tr_z_yy_xyz,   \
                             tr_z_yy_yyy,   \
                             tr_z_yy_yyz,   \
                             tr_z_yy_yzz,   \
                             tr_z_yy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_xxx[i] = tr_z_xx_xxx[i] * fe_0 + tr_z_xxy_xxx[i] * pa_y[i];

        tr_z_xxyy_xxy[i] = tr_z_yy_xxy[i] * fe_0 + 2.0 * tr_z_xyy_xy[i] * fe_0 + tr_z_xyy_xxy[i] * pa_x[i];

        tr_z_xxyy_xxz[i] = tr_z_xx_xxz[i] * fe_0 + tr_z_xxy_xxz[i] * pa_y[i];

        tr_z_xxyy_xyy[i] = tr_z_yy_xyy[i] * fe_0 + tr_z_xyy_yy[i] * fe_0 + tr_z_xyy_xyy[i] * pa_x[i];

        tr_z_xxyy_xyz[i] = tr_z_yy_xyz[i] * fe_0 + tr_z_xyy_yz[i] * fe_0 + tr_z_xyy_xyz[i] * pa_x[i];

        tr_z_xxyy_xzz[i] = tr_z_xx_xzz[i] * fe_0 + tr_z_xxy_xzz[i] * pa_y[i];

        tr_z_xxyy_yyy[i] = tr_z_yy_yyy[i] * fe_0 + tr_z_xyy_yyy[i] * pa_x[i];

        tr_z_xxyy_yyz[i] = tr_z_yy_yyz[i] * fe_0 + tr_z_xyy_yyz[i] * pa_x[i];

        tr_z_xxyy_yzz[i] = tr_z_yy_yzz[i] * fe_0 + tr_z_xyy_yzz[i] * pa_x[i];

        tr_z_xxyy_zzz[i] = tr_z_yy_zzz[i] * fe_0 + tr_z_xyy_zzz[i] * pa_x[i];
    }

    // Set up 340-350 components of targeted buffer : GF

    auto tr_z_xxyz_xxx = pbuffer.data(idx_dip_gf + 340);

    auto tr_z_xxyz_xxy = pbuffer.data(idx_dip_gf + 341);

    auto tr_z_xxyz_xxz = pbuffer.data(idx_dip_gf + 342);

    auto tr_z_xxyz_xyy = pbuffer.data(idx_dip_gf + 343);

    auto tr_z_xxyz_xyz = pbuffer.data(idx_dip_gf + 344);

    auto tr_z_xxyz_xzz = pbuffer.data(idx_dip_gf + 345);

    auto tr_z_xxyz_yyy = pbuffer.data(idx_dip_gf + 346);

    auto tr_z_xxyz_yyz = pbuffer.data(idx_dip_gf + 347);

    auto tr_z_xxyz_yzz = pbuffer.data(idx_dip_gf + 348);

    auto tr_z_xxyz_zzz = pbuffer.data(idx_dip_gf + 349);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xxyz_xxx, \
                             tr_z_xxyz_xxy, \
                             tr_z_xxyz_xxz, \
                             tr_z_xxyz_xyy, \
                             tr_z_xxyz_xyz, \
                             tr_z_xxyz_xzz, \
                             tr_z_xxyz_yyy, \
                             tr_z_xxyz_yyz, \
                             tr_z_xxyz_yzz, \
                             tr_z_xxyz_zzz, \
                             tr_z_xxz_xx,   \
                             tr_z_xxz_xxx,  \
                             tr_z_xxz_xxy,  \
                             tr_z_xxz_xxz,  \
                             tr_z_xxz_xy,   \
                             tr_z_xxz_xyy,  \
                             tr_z_xxz_xyz,  \
                             tr_z_xxz_xz,   \
                             tr_z_xxz_xzz,  \
                             tr_z_xxz_zzz,  \
                             tr_z_xyz_yyy,  \
                             tr_z_xyz_yyz,  \
                             tr_z_xyz_yzz,  \
                             tr_z_yz_yyy,   \
                             tr_z_yz_yyz,   \
                             tr_z_yz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_xxx[i] = tr_z_xxz_xxx[i] * pa_y[i];

        tr_z_xxyz_xxy[i] = tr_z_xxz_xx[i] * fe_0 + tr_z_xxz_xxy[i] * pa_y[i];

        tr_z_xxyz_xxz[i] = tr_z_xxz_xxz[i] * pa_y[i];

        tr_z_xxyz_xyy[i] = 2.0 * tr_z_xxz_xy[i] * fe_0 + tr_z_xxz_xyy[i] * pa_y[i];

        tr_z_xxyz_xyz[i] = tr_z_xxz_xz[i] * fe_0 + tr_z_xxz_xyz[i] * pa_y[i];

        tr_z_xxyz_xzz[i] = tr_z_xxz_xzz[i] * pa_y[i];

        tr_z_xxyz_yyy[i] = tr_z_yz_yyy[i] * fe_0 + tr_z_xyz_yyy[i] * pa_x[i];

        tr_z_xxyz_yyz[i] = tr_z_yz_yyz[i] * fe_0 + tr_z_xyz_yyz[i] * pa_x[i];

        tr_z_xxyz_yzz[i] = tr_z_yz_yzz[i] * fe_0 + tr_z_xyz_yzz[i] * pa_x[i];

        tr_z_xxyz_zzz[i] = tr_z_xxz_zzz[i] * pa_y[i];
    }

    // Set up 350-360 components of targeted buffer : GF

    auto tr_z_xxzz_xxx = pbuffer.data(idx_dip_gf + 350);

    auto tr_z_xxzz_xxy = pbuffer.data(idx_dip_gf + 351);

    auto tr_z_xxzz_xxz = pbuffer.data(idx_dip_gf + 352);

    auto tr_z_xxzz_xyy = pbuffer.data(idx_dip_gf + 353);

    auto tr_z_xxzz_xyz = pbuffer.data(idx_dip_gf + 354);

    auto tr_z_xxzz_xzz = pbuffer.data(idx_dip_gf + 355);

    auto tr_z_xxzz_yyy = pbuffer.data(idx_dip_gf + 356);

    auto tr_z_xxzz_yyz = pbuffer.data(idx_dip_gf + 357);

    auto tr_z_xxzz_yzz = pbuffer.data(idx_dip_gf + 358);

    auto tr_z_xxzz_zzz = pbuffer.data(idx_dip_gf + 359);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xxzz_xxx, \
                             tr_z_xxzz_xxy, \
                             tr_z_xxzz_xxz, \
                             tr_z_xxzz_xyy, \
                             tr_z_xxzz_xyz, \
                             tr_z_xxzz_xzz, \
                             tr_z_xxzz_yyy, \
                             tr_z_xxzz_yyz, \
                             tr_z_xxzz_yzz, \
                             tr_z_xxzz_zzz, \
                             tr_z_xzz_xx,   \
                             tr_z_xzz_xxx,  \
                             tr_z_xzz_xxy,  \
                             tr_z_xzz_xxz,  \
                             tr_z_xzz_xy,   \
                             tr_z_xzz_xyy,  \
                             tr_z_xzz_xyz,  \
                             tr_z_xzz_xz,   \
                             tr_z_xzz_xzz,  \
                             tr_z_xzz_yy,   \
                             tr_z_xzz_yyy,  \
                             tr_z_xzz_yyz,  \
                             tr_z_xzz_yz,   \
                             tr_z_xzz_yzz,  \
                             tr_z_xzz_zz,   \
                             tr_z_xzz_zzz,  \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xzz,   \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yzz,   \
                             tr_z_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_xxx[i] = tr_z_zz_xxx[i] * fe_0 + 3.0 * tr_z_xzz_xx[i] * fe_0 + tr_z_xzz_xxx[i] * pa_x[i];

        tr_z_xxzz_xxy[i] = tr_z_zz_xxy[i] * fe_0 + 2.0 * tr_z_xzz_xy[i] * fe_0 + tr_z_xzz_xxy[i] * pa_x[i];

        tr_z_xxzz_xxz[i] = tr_z_zz_xxz[i] * fe_0 + 2.0 * tr_z_xzz_xz[i] * fe_0 + tr_z_xzz_xxz[i] * pa_x[i];

        tr_z_xxzz_xyy[i] = tr_z_zz_xyy[i] * fe_0 + tr_z_xzz_yy[i] * fe_0 + tr_z_xzz_xyy[i] * pa_x[i];

        tr_z_xxzz_xyz[i] = tr_z_zz_xyz[i] * fe_0 + tr_z_xzz_yz[i] * fe_0 + tr_z_xzz_xyz[i] * pa_x[i];

        tr_z_xxzz_xzz[i] = tr_z_zz_xzz[i] * fe_0 + tr_z_xzz_zz[i] * fe_0 + tr_z_xzz_xzz[i] * pa_x[i];

        tr_z_xxzz_yyy[i] = tr_z_zz_yyy[i] * fe_0 + tr_z_xzz_yyy[i] * pa_x[i];

        tr_z_xxzz_yyz[i] = tr_z_zz_yyz[i] * fe_0 + tr_z_xzz_yyz[i] * pa_x[i];

        tr_z_xxzz_yzz[i] = tr_z_zz_yzz[i] * fe_0 + tr_z_xzz_yzz[i] * pa_x[i];

        tr_z_xxzz_zzz[i] = tr_z_zz_zzz[i] * fe_0 + tr_z_xzz_zzz[i] * pa_x[i];
    }

    // Set up 360-370 components of targeted buffer : GF

    auto tr_z_xyyy_xxx = pbuffer.data(idx_dip_gf + 360);

    auto tr_z_xyyy_xxy = pbuffer.data(idx_dip_gf + 361);

    auto tr_z_xyyy_xxz = pbuffer.data(idx_dip_gf + 362);

    auto tr_z_xyyy_xyy = pbuffer.data(idx_dip_gf + 363);

    auto tr_z_xyyy_xyz = pbuffer.data(idx_dip_gf + 364);

    auto tr_z_xyyy_xzz = pbuffer.data(idx_dip_gf + 365);

    auto tr_z_xyyy_yyy = pbuffer.data(idx_dip_gf + 366);

    auto tr_z_xyyy_yyz = pbuffer.data(idx_dip_gf + 367);

    auto tr_z_xyyy_yzz = pbuffer.data(idx_dip_gf + 368);

    auto tr_z_xyyy_zzz = pbuffer.data(idx_dip_gf + 369);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xyyy_xxx, \
                             tr_z_xyyy_xxy, \
                             tr_z_xyyy_xxz, \
                             tr_z_xyyy_xyy, \
                             tr_z_xyyy_xyz, \
                             tr_z_xyyy_xzz, \
                             tr_z_xyyy_yyy, \
                             tr_z_xyyy_yyz, \
                             tr_z_xyyy_yzz, \
                             tr_z_xyyy_zzz, \
                             tr_z_yyy_xx,   \
                             tr_z_yyy_xxx,  \
                             tr_z_yyy_xxy,  \
                             tr_z_yyy_xxz,  \
                             tr_z_yyy_xy,   \
                             tr_z_yyy_xyy,  \
                             tr_z_yyy_xyz,  \
                             tr_z_yyy_xz,   \
                             tr_z_yyy_xzz,  \
                             tr_z_yyy_yy,   \
                             tr_z_yyy_yyy,  \
                             tr_z_yyy_yyz,  \
                             tr_z_yyy_yz,   \
                             tr_z_yyy_yzz,  \
                             tr_z_yyy_zz,   \
                             tr_z_yyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_xxx[i] = 3.0 * tr_z_yyy_xx[i] * fe_0 + tr_z_yyy_xxx[i] * pa_x[i];

        tr_z_xyyy_xxy[i] = 2.0 * tr_z_yyy_xy[i] * fe_0 + tr_z_yyy_xxy[i] * pa_x[i];

        tr_z_xyyy_xxz[i] = 2.0 * tr_z_yyy_xz[i] * fe_0 + tr_z_yyy_xxz[i] * pa_x[i];

        tr_z_xyyy_xyy[i] = tr_z_yyy_yy[i] * fe_0 + tr_z_yyy_xyy[i] * pa_x[i];

        tr_z_xyyy_xyz[i] = tr_z_yyy_yz[i] * fe_0 + tr_z_yyy_xyz[i] * pa_x[i];

        tr_z_xyyy_xzz[i] = tr_z_yyy_zz[i] * fe_0 + tr_z_yyy_xzz[i] * pa_x[i];

        tr_z_xyyy_yyy[i] = tr_z_yyy_yyy[i] * pa_x[i];

        tr_z_xyyy_yyz[i] = tr_z_yyy_yyz[i] * pa_x[i];

        tr_z_xyyy_yzz[i] = tr_z_yyy_yzz[i] * pa_x[i];

        tr_z_xyyy_zzz[i] = tr_z_yyy_zzz[i] * pa_x[i];
    }

    // Set up 370-380 components of targeted buffer : GF

    auto tr_z_xyyz_xxx = pbuffer.data(idx_dip_gf + 370);

    auto tr_z_xyyz_xxy = pbuffer.data(idx_dip_gf + 371);

    auto tr_z_xyyz_xxz = pbuffer.data(idx_dip_gf + 372);

    auto tr_z_xyyz_xyy = pbuffer.data(idx_dip_gf + 373);

    auto tr_z_xyyz_xyz = pbuffer.data(idx_dip_gf + 374);

    auto tr_z_xyyz_xzz = pbuffer.data(idx_dip_gf + 375);

    auto tr_z_xyyz_yyy = pbuffer.data(idx_dip_gf + 376);

    auto tr_z_xyyz_yyz = pbuffer.data(idx_dip_gf + 377);

    auto tr_z_xyyz_yzz = pbuffer.data(idx_dip_gf + 378);

    auto tr_z_xyyz_zzz = pbuffer.data(idx_dip_gf + 379);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xyyz_xxx, \
                             tr_z_xyyz_xxy, \
                             tr_z_xyyz_xxz, \
                             tr_z_xyyz_xyy, \
                             tr_z_xyyz_xyz, \
                             tr_z_xyyz_xzz, \
                             tr_z_xyyz_yyy, \
                             tr_z_xyyz_yyz, \
                             tr_z_xyyz_yzz, \
                             tr_z_xyyz_zzz, \
                             tr_z_yyz_xx,   \
                             tr_z_yyz_xxx,  \
                             tr_z_yyz_xxy,  \
                             tr_z_yyz_xxz,  \
                             tr_z_yyz_xy,   \
                             tr_z_yyz_xyy,  \
                             tr_z_yyz_xyz,  \
                             tr_z_yyz_xz,   \
                             tr_z_yyz_xzz,  \
                             tr_z_yyz_yy,   \
                             tr_z_yyz_yyy,  \
                             tr_z_yyz_yyz,  \
                             tr_z_yyz_yz,   \
                             tr_z_yyz_yzz,  \
                             tr_z_yyz_zz,   \
                             tr_z_yyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_xxx[i] = 3.0 * tr_z_yyz_xx[i] * fe_0 + tr_z_yyz_xxx[i] * pa_x[i];

        tr_z_xyyz_xxy[i] = 2.0 * tr_z_yyz_xy[i] * fe_0 + tr_z_yyz_xxy[i] * pa_x[i];

        tr_z_xyyz_xxz[i] = 2.0 * tr_z_yyz_xz[i] * fe_0 + tr_z_yyz_xxz[i] * pa_x[i];

        tr_z_xyyz_xyy[i] = tr_z_yyz_yy[i] * fe_0 + tr_z_yyz_xyy[i] * pa_x[i];

        tr_z_xyyz_xyz[i] = tr_z_yyz_yz[i] * fe_0 + tr_z_yyz_xyz[i] * pa_x[i];

        tr_z_xyyz_xzz[i] = tr_z_yyz_zz[i] * fe_0 + tr_z_yyz_xzz[i] * pa_x[i];

        tr_z_xyyz_yyy[i] = tr_z_yyz_yyy[i] * pa_x[i];

        tr_z_xyyz_yyz[i] = tr_z_yyz_yyz[i] * pa_x[i];

        tr_z_xyyz_yzz[i] = tr_z_yyz_yzz[i] * pa_x[i];

        tr_z_xyyz_zzz[i] = tr_z_yyz_zzz[i] * pa_x[i];
    }

    // Set up 380-390 components of targeted buffer : GF

    auto tr_z_xyzz_xxx = pbuffer.data(idx_dip_gf + 380);

    auto tr_z_xyzz_xxy = pbuffer.data(idx_dip_gf + 381);

    auto tr_z_xyzz_xxz = pbuffer.data(idx_dip_gf + 382);

    auto tr_z_xyzz_xyy = pbuffer.data(idx_dip_gf + 383);

    auto tr_z_xyzz_xyz = pbuffer.data(idx_dip_gf + 384);

    auto tr_z_xyzz_xzz = pbuffer.data(idx_dip_gf + 385);

    auto tr_z_xyzz_yyy = pbuffer.data(idx_dip_gf + 386);

    auto tr_z_xyzz_yyz = pbuffer.data(idx_dip_gf + 387);

    auto tr_z_xyzz_yzz = pbuffer.data(idx_dip_gf + 388);

    auto tr_z_xyzz_zzz = pbuffer.data(idx_dip_gf + 389);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             tr_z_xyzz_xxx, \
                             tr_z_xyzz_xxy, \
                             tr_z_xyzz_xxz, \
                             tr_z_xyzz_xyy, \
                             tr_z_xyzz_xyz, \
                             tr_z_xyzz_xzz, \
                             tr_z_xyzz_yyy, \
                             tr_z_xyzz_yyz, \
                             tr_z_xyzz_yzz, \
                             tr_z_xyzz_zzz, \
                             tr_z_xzz_xxx,  \
                             tr_z_xzz_xxz,  \
                             tr_z_xzz_xzz,  \
                             tr_z_yzz_xxy,  \
                             tr_z_yzz_xy,   \
                             tr_z_yzz_xyy,  \
                             tr_z_yzz_xyz,  \
                             tr_z_yzz_yy,   \
                             tr_z_yzz_yyy,  \
                             tr_z_yzz_yyz,  \
                             tr_z_yzz_yz,   \
                             tr_z_yzz_yzz,  \
                             tr_z_yzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzz_xxx[i] = tr_z_xzz_xxx[i] * pa_y[i];

        tr_z_xyzz_xxy[i] = 2.0 * tr_z_yzz_xy[i] * fe_0 + tr_z_yzz_xxy[i] * pa_x[i];

        tr_z_xyzz_xxz[i] = tr_z_xzz_xxz[i] * pa_y[i];

        tr_z_xyzz_xyy[i] = tr_z_yzz_yy[i] * fe_0 + tr_z_yzz_xyy[i] * pa_x[i];

        tr_z_xyzz_xyz[i] = tr_z_yzz_yz[i] * fe_0 + tr_z_yzz_xyz[i] * pa_x[i];

        tr_z_xyzz_xzz[i] = tr_z_xzz_xzz[i] * pa_y[i];

        tr_z_xyzz_yyy[i] = tr_z_yzz_yyy[i] * pa_x[i];

        tr_z_xyzz_yyz[i] = tr_z_yzz_yyz[i] * pa_x[i];

        tr_z_xyzz_yzz[i] = tr_z_yzz_yzz[i] * pa_x[i];

        tr_z_xyzz_zzz[i] = tr_z_yzz_zzz[i] * pa_x[i];
    }

    // Set up 390-400 components of targeted buffer : GF

    auto tr_z_xzzz_xxx = pbuffer.data(idx_dip_gf + 390);

    auto tr_z_xzzz_xxy = pbuffer.data(idx_dip_gf + 391);

    auto tr_z_xzzz_xxz = pbuffer.data(idx_dip_gf + 392);

    auto tr_z_xzzz_xyy = pbuffer.data(idx_dip_gf + 393);

    auto tr_z_xzzz_xyz = pbuffer.data(idx_dip_gf + 394);

    auto tr_z_xzzz_xzz = pbuffer.data(idx_dip_gf + 395);

    auto tr_z_xzzz_yyy = pbuffer.data(idx_dip_gf + 396);

    auto tr_z_xzzz_yyz = pbuffer.data(idx_dip_gf + 397);

    auto tr_z_xzzz_yzz = pbuffer.data(idx_dip_gf + 398);

    auto tr_z_xzzz_zzz = pbuffer.data(idx_dip_gf + 399);

#pragma omp simd aligned(pa_x,              \
                             tr_z_xzzz_xxx, \
                             tr_z_xzzz_xxy, \
                             tr_z_xzzz_xxz, \
                             tr_z_xzzz_xyy, \
                             tr_z_xzzz_xyz, \
                             tr_z_xzzz_xzz, \
                             tr_z_xzzz_yyy, \
                             tr_z_xzzz_yyz, \
                             tr_z_xzzz_yzz, \
                             tr_z_xzzz_zzz, \
                             tr_z_zzz_xx,   \
                             tr_z_zzz_xxx,  \
                             tr_z_zzz_xxy,  \
                             tr_z_zzz_xxz,  \
                             tr_z_zzz_xy,   \
                             tr_z_zzz_xyy,  \
                             tr_z_zzz_xyz,  \
                             tr_z_zzz_xz,   \
                             tr_z_zzz_xzz,  \
                             tr_z_zzz_yy,   \
                             tr_z_zzz_yyy,  \
                             tr_z_zzz_yyz,  \
                             tr_z_zzz_yz,   \
                             tr_z_zzz_yzz,  \
                             tr_z_zzz_zz,   \
                             tr_z_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_xxx[i] = 3.0 * tr_z_zzz_xx[i] * fe_0 + tr_z_zzz_xxx[i] * pa_x[i];

        tr_z_xzzz_xxy[i] = 2.0 * tr_z_zzz_xy[i] * fe_0 + tr_z_zzz_xxy[i] * pa_x[i];

        tr_z_xzzz_xxz[i] = 2.0 * tr_z_zzz_xz[i] * fe_0 + tr_z_zzz_xxz[i] * pa_x[i];

        tr_z_xzzz_xyy[i] = tr_z_zzz_yy[i] * fe_0 + tr_z_zzz_xyy[i] * pa_x[i];

        tr_z_xzzz_xyz[i] = tr_z_zzz_yz[i] * fe_0 + tr_z_zzz_xyz[i] * pa_x[i];

        tr_z_xzzz_xzz[i] = tr_z_zzz_zz[i] * fe_0 + tr_z_zzz_xzz[i] * pa_x[i];

        tr_z_xzzz_yyy[i] = tr_z_zzz_yyy[i] * pa_x[i];

        tr_z_xzzz_yyz[i] = tr_z_zzz_yyz[i] * pa_x[i];

        tr_z_xzzz_yzz[i] = tr_z_zzz_yzz[i] * pa_x[i];

        tr_z_xzzz_zzz[i] = tr_z_zzz_zzz[i] * pa_x[i];
    }

    // Set up 400-410 components of targeted buffer : GF

    auto tr_z_yyyy_xxx = pbuffer.data(idx_dip_gf + 400);

    auto tr_z_yyyy_xxy = pbuffer.data(idx_dip_gf + 401);

    auto tr_z_yyyy_xxz = pbuffer.data(idx_dip_gf + 402);

    auto tr_z_yyyy_xyy = pbuffer.data(idx_dip_gf + 403);

    auto tr_z_yyyy_xyz = pbuffer.data(idx_dip_gf + 404);

    auto tr_z_yyyy_xzz = pbuffer.data(idx_dip_gf + 405);

    auto tr_z_yyyy_yyy = pbuffer.data(idx_dip_gf + 406);

    auto tr_z_yyyy_yyz = pbuffer.data(idx_dip_gf + 407);

    auto tr_z_yyyy_yzz = pbuffer.data(idx_dip_gf + 408);

    auto tr_z_yyyy_zzz = pbuffer.data(idx_dip_gf + 409);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yy_xxx,   \
                             tr_z_yy_xxy,   \
                             tr_z_yy_xxz,   \
                             tr_z_yy_xyy,   \
                             tr_z_yy_xyz,   \
                             tr_z_yy_xzz,   \
                             tr_z_yy_yyy,   \
                             tr_z_yy_yyz,   \
                             tr_z_yy_yzz,   \
                             tr_z_yy_zzz,   \
                             tr_z_yyy_xx,   \
                             tr_z_yyy_xxx,  \
                             tr_z_yyy_xxy,  \
                             tr_z_yyy_xxz,  \
                             tr_z_yyy_xy,   \
                             tr_z_yyy_xyy,  \
                             tr_z_yyy_xyz,  \
                             tr_z_yyy_xz,   \
                             tr_z_yyy_xzz,  \
                             tr_z_yyy_yy,   \
                             tr_z_yyy_yyy,  \
                             tr_z_yyy_yyz,  \
                             tr_z_yyy_yz,   \
                             tr_z_yyy_yzz,  \
                             tr_z_yyy_zz,   \
                             tr_z_yyy_zzz,  \
                             tr_z_yyyy_xxx, \
                             tr_z_yyyy_xxy, \
                             tr_z_yyyy_xxz, \
                             tr_z_yyyy_xyy, \
                             tr_z_yyyy_xyz, \
                             tr_z_yyyy_xzz, \
                             tr_z_yyyy_yyy, \
                             tr_z_yyyy_yyz, \
                             tr_z_yyyy_yzz, \
                             tr_z_yyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_xxx[i] = 3.0 * tr_z_yy_xxx[i] * fe_0 + tr_z_yyy_xxx[i] * pa_y[i];

        tr_z_yyyy_xxy[i] = 3.0 * tr_z_yy_xxy[i] * fe_0 + tr_z_yyy_xx[i] * fe_0 + tr_z_yyy_xxy[i] * pa_y[i];

        tr_z_yyyy_xxz[i] = 3.0 * tr_z_yy_xxz[i] * fe_0 + tr_z_yyy_xxz[i] * pa_y[i];

        tr_z_yyyy_xyy[i] = 3.0 * tr_z_yy_xyy[i] * fe_0 + 2.0 * tr_z_yyy_xy[i] * fe_0 + tr_z_yyy_xyy[i] * pa_y[i];

        tr_z_yyyy_xyz[i] = 3.0 * tr_z_yy_xyz[i] * fe_0 + tr_z_yyy_xz[i] * fe_0 + tr_z_yyy_xyz[i] * pa_y[i];

        tr_z_yyyy_xzz[i] = 3.0 * tr_z_yy_xzz[i] * fe_0 + tr_z_yyy_xzz[i] * pa_y[i];

        tr_z_yyyy_yyy[i] = 3.0 * tr_z_yy_yyy[i] * fe_0 + 3.0 * tr_z_yyy_yy[i] * fe_0 + tr_z_yyy_yyy[i] * pa_y[i];

        tr_z_yyyy_yyz[i] = 3.0 * tr_z_yy_yyz[i] * fe_0 + 2.0 * tr_z_yyy_yz[i] * fe_0 + tr_z_yyy_yyz[i] * pa_y[i];

        tr_z_yyyy_yzz[i] = 3.0 * tr_z_yy_yzz[i] * fe_0 + tr_z_yyy_zz[i] * fe_0 + tr_z_yyy_yzz[i] * pa_y[i];

        tr_z_yyyy_zzz[i] = 3.0 * tr_z_yy_zzz[i] * fe_0 + tr_z_yyy_zzz[i] * pa_y[i];
    }

    // Set up 410-420 components of targeted buffer : GF

    auto tr_z_yyyz_xxx = pbuffer.data(idx_dip_gf + 410);

    auto tr_z_yyyz_xxy = pbuffer.data(idx_dip_gf + 411);

    auto tr_z_yyyz_xxz = pbuffer.data(idx_dip_gf + 412);

    auto tr_z_yyyz_xyy = pbuffer.data(idx_dip_gf + 413);

    auto tr_z_yyyz_xyz = pbuffer.data(idx_dip_gf + 414);

    auto tr_z_yyyz_xzz = pbuffer.data(idx_dip_gf + 415);

    auto tr_z_yyyz_yyy = pbuffer.data(idx_dip_gf + 416);

    auto tr_z_yyyz_yyz = pbuffer.data(idx_dip_gf + 417);

    auto tr_z_yyyz_yzz = pbuffer.data(idx_dip_gf + 418);

    auto tr_z_yyyz_zzz = pbuffer.data(idx_dip_gf + 419);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             tr_z_yyy_xxy,  \
                             tr_z_yyy_xyy,  \
                             tr_z_yyy_yyy,  \
                             tr_z_yyyz_xxx, \
                             tr_z_yyyz_xxy, \
                             tr_z_yyyz_xxz, \
                             tr_z_yyyz_xyy, \
                             tr_z_yyyz_xyz, \
                             tr_z_yyyz_xzz, \
                             tr_z_yyyz_yyy, \
                             tr_z_yyyz_yyz, \
                             tr_z_yyyz_yzz, \
                             tr_z_yyyz_zzz, \
                             tr_z_yyz_xxx,  \
                             tr_z_yyz_xxz,  \
                             tr_z_yyz_xyz,  \
                             tr_z_yyz_xz,   \
                             tr_z_yyz_xzz,  \
                             tr_z_yyz_yyz,  \
                             tr_z_yyz_yz,   \
                             tr_z_yyz_yzz,  \
                             tr_z_yyz_zz,   \
                             tr_z_yyz_zzz,  \
                             tr_z_yz_xxx,   \
                             tr_z_yz_xxz,   \
                             tr_z_yz_xyz,   \
                             tr_z_yz_xzz,   \
                             tr_z_yz_yyz,   \
                             tr_z_yz_yzz,   \
                             tr_z_yz_zzz,   \
                             ts_yyy_xxy,    \
                             ts_yyy_xyy,    \
                             ts_yyy_yyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_xxx[i] = 2.0 * tr_z_yz_xxx[i] * fe_0 + tr_z_yyz_xxx[i] * pa_y[i];

        tr_z_yyyz_xxy[i] = ts_yyy_xxy[i] * fe_0 + tr_z_yyy_xxy[i] * pa_z[i];

        tr_z_yyyz_xxz[i] = 2.0 * tr_z_yz_xxz[i] * fe_0 + tr_z_yyz_xxz[i] * pa_y[i];

        tr_z_yyyz_xyy[i] = ts_yyy_xyy[i] * fe_0 + tr_z_yyy_xyy[i] * pa_z[i];

        tr_z_yyyz_xyz[i] = 2.0 * tr_z_yz_xyz[i] * fe_0 + tr_z_yyz_xz[i] * fe_0 + tr_z_yyz_xyz[i] * pa_y[i];

        tr_z_yyyz_xzz[i] = 2.0 * tr_z_yz_xzz[i] * fe_0 + tr_z_yyz_xzz[i] * pa_y[i];

        tr_z_yyyz_yyy[i] = ts_yyy_yyy[i] * fe_0 + tr_z_yyy_yyy[i] * pa_z[i];

        tr_z_yyyz_yyz[i] = 2.0 * tr_z_yz_yyz[i] * fe_0 + 2.0 * tr_z_yyz_yz[i] * fe_0 + tr_z_yyz_yyz[i] * pa_y[i];

        tr_z_yyyz_yzz[i] = 2.0 * tr_z_yz_yzz[i] * fe_0 + tr_z_yyz_zz[i] * fe_0 + tr_z_yyz_yzz[i] * pa_y[i];

        tr_z_yyyz_zzz[i] = 2.0 * tr_z_yz_zzz[i] * fe_0 + tr_z_yyz_zzz[i] * pa_y[i];
    }

    // Set up 420-430 components of targeted buffer : GF

    auto tr_z_yyzz_xxx = pbuffer.data(idx_dip_gf + 420);

    auto tr_z_yyzz_xxy = pbuffer.data(idx_dip_gf + 421);

    auto tr_z_yyzz_xxz = pbuffer.data(idx_dip_gf + 422);

    auto tr_z_yyzz_xyy = pbuffer.data(idx_dip_gf + 423);

    auto tr_z_yyzz_xyz = pbuffer.data(idx_dip_gf + 424);

    auto tr_z_yyzz_xzz = pbuffer.data(idx_dip_gf + 425);

    auto tr_z_yyzz_yyy = pbuffer.data(idx_dip_gf + 426);

    auto tr_z_yyzz_yyz = pbuffer.data(idx_dip_gf + 427);

    auto tr_z_yyzz_yzz = pbuffer.data(idx_dip_gf + 428);

    auto tr_z_yyzz_zzz = pbuffer.data(idx_dip_gf + 429);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yyzz_xxx, \
                             tr_z_yyzz_xxy, \
                             tr_z_yyzz_xxz, \
                             tr_z_yyzz_xyy, \
                             tr_z_yyzz_xyz, \
                             tr_z_yyzz_xzz, \
                             tr_z_yyzz_yyy, \
                             tr_z_yyzz_yyz, \
                             tr_z_yyzz_yzz, \
                             tr_z_yyzz_zzz, \
                             tr_z_yzz_xx,   \
                             tr_z_yzz_xxx,  \
                             tr_z_yzz_xxy,  \
                             tr_z_yzz_xxz,  \
                             tr_z_yzz_xy,   \
                             tr_z_yzz_xyy,  \
                             tr_z_yzz_xyz,  \
                             tr_z_yzz_xz,   \
                             tr_z_yzz_xzz,  \
                             tr_z_yzz_yy,   \
                             tr_z_yzz_yyy,  \
                             tr_z_yzz_yyz,  \
                             tr_z_yzz_yz,   \
                             tr_z_yzz_yzz,  \
                             tr_z_yzz_zz,   \
                             tr_z_yzz_zzz,  \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xzz,   \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yzz,   \
                             tr_z_zz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_xxx[i] = tr_z_zz_xxx[i] * fe_0 + tr_z_yzz_xxx[i] * pa_y[i];

        tr_z_yyzz_xxy[i] = tr_z_zz_xxy[i] * fe_0 + tr_z_yzz_xx[i] * fe_0 + tr_z_yzz_xxy[i] * pa_y[i];

        tr_z_yyzz_xxz[i] = tr_z_zz_xxz[i] * fe_0 + tr_z_yzz_xxz[i] * pa_y[i];

        tr_z_yyzz_xyy[i] = tr_z_zz_xyy[i] * fe_0 + 2.0 * tr_z_yzz_xy[i] * fe_0 + tr_z_yzz_xyy[i] * pa_y[i];

        tr_z_yyzz_xyz[i] = tr_z_zz_xyz[i] * fe_0 + tr_z_yzz_xz[i] * fe_0 + tr_z_yzz_xyz[i] * pa_y[i];

        tr_z_yyzz_xzz[i] = tr_z_zz_xzz[i] * fe_0 + tr_z_yzz_xzz[i] * pa_y[i];

        tr_z_yyzz_yyy[i] = tr_z_zz_yyy[i] * fe_0 + 3.0 * tr_z_yzz_yy[i] * fe_0 + tr_z_yzz_yyy[i] * pa_y[i];

        tr_z_yyzz_yyz[i] = tr_z_zz_yyz[i] * fe_0 + 2.0 * tr_z_yzz_yz[i] * fe_0 + tr_z_yzz_yyz[i] * pa_y[i];

        tr_z_yyzz_yzz[i] = tr_z_zz_yzz[i] * fe_0 + tr_z_yzz_zz[i] * fe_0 + tr_z_yzz_yzz[i] * pa_y[i];

        tr_z_yyzz_zzz[i] = tr_z_zz_zzz[i] * fe_0 + tr_z_yzz_zzz[i] * pa_y[i];
    }

    // Set up 430-440 components of targeted buffer : GF

    auto tr_z_yzzz_xxx = pbuffer.data(idx_dip_gf + 430);

    auto tr_z_yzzz_xxy = pbuffer.data(idx_dip_gf + 431);

    auto tr_z_yzzz_xxz = pbuffer.data(idx_dip_gf + 432);

    auto tr_z_yzzz_xyy = pbuffer.data(idx_dip_gf + 433);

    auto tr_z_yzzz_xyz = pbuffer.data(idx_dip_gf + 434);

    auto tr_z_yzzz_xzz = pbuffer.data(idx_dip_gf + 435);

    auto tr_z_yzzz_yyy = pbuffer.data(idx_dip_gf + 436);

    auto tr_z_yzzz_yyz = pbuffer.data(idx_dip_gf + 437);

    auto tr_z_yzzz_yzz = pbuffer.data(idx_dip_gf + 438);

    auto tr_z_yzzz_zzz = pbuffer.data(idx_dip_gf + 439);

#pragma omp simd aligned(pa_y,              \
                             tr_z_yzzz_xxx, \
                             tr_z_yzzz_xxy, \
                             tr_z_yzzz_xxz, \
                             tr_z_yzzz_xyy, \
                             tr_z_yzzz_xyz, \
                             tr_z_yzzz_xzz, \
                             tr_z_yzzz_yyy, \
                             tr_z_yzzz_yyz, \
                             tr_z_yzzz_yzz, \
                             tr_z_yzzz_zzz, \
                             tr_z_zzz_xx,   \
                             tr_z_zzz_xxx,  \
                             tr_z_zzz_xxy,  \
                             tr_z_zzz_xxz,  \
                             tr_z_zzz_xy,   \
                             tr_z_zzz_xyy,  \
                             tr_z_zzz_xyz,  \
                             tr_z_zzz_xz,   \
                             tr_z_zzz_xzz,  \
                             tr_z_zzz_yy,   \
                             tr_z_zzz_yyy,  \
                             tr_z_zzz_yyz,  \
                             tr_z_zzz_yz,   \
                             tr_z_zzz_yzz,  \
                             tr_z_zzz_zz,   \
                             tr_z_zzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_xxx[i] = tr_z_zzz_xxx[i] * pa_y[i];

        tr_z_yzzz_xxy[i] = tr_z_zzz_xx[i] * fe_0 + tr_z_zzz_xxy[i] * pa_y[i];

        tr_z_yzzz_xxz[i] = tr_z_zzz_xxz[i] * pa_y[i];

        tr_z_yzzz_xyy[i] = 2.0 * tr_z_zzz_xy[i] * fe_0 + tr_z_zzz_xyy[i] * pa_y[i];

        tr_z_yzzz_xyz[i] = tr_z_zzz_xz[i] * fe_0 + tr_z_zzz_xyz[i] * pa_y[i];

        tr_z_yzzz_xzz[i] = tr_z_zzz_xzz[i] * pa_y[i];

        tr_z_yzzz_yyy[i] = 3.0 * tr_z_zzz_yy[i] * fe_0 + tr_z_zzz_yyy[i] * pa_y[i];

        tr_z_yzzz_yyz[i] = 2.0 * tr_z_zzz_yz[i] * fe_0 + tr_z_zzz_yyz[i] * pa_y[i];

        tr_z_yzzz_yzz[i] = tr_z_zzz_zz[i] * fe_0 + tr_z_zzz_yzz[i] * pa_y[i];

        tr_z_yzzz_zzz[i] = tr_z_zzz_zzz[i] * pa_y[i];
    }

    // Set up 440-450 components of targeted buffer : GF

    auto tr_z_zzzz_xxx = pbuffer.data(idx_dip_gf + 440);

    auto tr_z_zzzz_xxy = pbuffer.data(idx_dip_gf + 441);

    auto tr_z_zzzz_xxz = pbuffer.data(idx_dip_gf + 442);

    auto tr_z_zzzz_xyy = pbuffer.data(idx_dip_gf + 443);

    auto tr_z_zzzz_xyz = pbuffer.data(idx_dip_gf + 444);

    auto tr_z_zzzz_xzz = pbuffer.data(idx_dip_gf + 445);

    auto tr_z_zzzz_yyy = pbuffer.data(idx_dip_gf + 446);

    auto tr_z_zzzz_yyz = pbuffer.data(idx_dip_gf + 447);

    auto tr_z_zzzz_yzz = pbuffer.data(idx_dip_gf + 448);

    auto tr_z_zzzz_zzz = pbuffer.data(idx_dip_gf + 449);

#pragma omp simd aligned(pa_z,              \
                             tr_z_zz_xxx,   \
                             tr_z_zz_xxy,   \
                             tr_z_zz_xxz,   \
                             tr_z_zz_xyy,   \
                             tr_z_zz_xyz,   \
                             tr_z_zz_xzz,   \
                             tr_z_zz_yyy,   \
                             tr_z_zz_yyz,   \
                             tr_z_zz_yzz,   \
                             tr_z_zz_zzz,   \
                             tr_z_zzz_xx,   \
                             tr_z_zzz_xxx,  \
                             tr_z_zzz_xxy,  \
                             tr_z_zzz_xxz,  \
                             tr_z_zzz_xy,   \
                             tr_z_zzz_xyy,  \
                             tr_z_zzz_xyz,  \
                             tr_z_zzz_xz,   \
                             tr_z_zzz_xzz,  \
                             tr_z_zzz_yy,   \
                             tr_z_zzz_yyy,  \
                             tr_z_zzz_yyz,  \
                             tr_z_zzz_yz,   \
                             tr_z_zzz_yzz,  \
                             tr_z_zzz_zz,   \
                             tr_z_zzz_zzz,  \
                             tr_z_zzzz_xxx, \
                             tr_z_zzzz_xxy, \
                             tr_z_zzzz_xxz, \
                             tr_z_zzzz_xyy, \
                             tr_z_zzzz_xyz, \
                             tr_z_zzzz_xzz, \
                             tr_z_zzzz_yyy, \
                             tr_z_zzzz_yyz, \
                             tr_z_zzzz_yzz, \
                             tr_z_zzzz_zzz, \
                             ts_zzz_xxx,    \
                             ts_zzz_xxy,    \
                             ts_zzz_xxz,    \
                             ts_zzz_xyy,    \
                             ts_zzz_xyz,    \
                             ts_zzz_xzz,    \
                             ts_zzz_yyy,    \
                             ts_zzz_yyz,    \
                             ts_zzz_yzz,    \
                             ts_zzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_xxx[i] = 3.0 * tr_z_zz_xxx[i] * fe_0 + ts_zzz_xxx[i] * fe_0 + tr_z_zzz_xxx[i] * pa_z[i];

        tr_z_zzzz_xxy[i] = 3.0 * tr_z_zz_xxy[i] * fe_0 + ts_zzz_xxy[i] * fe_0 + tr_z_zzz_xxy[i] * pa_z[i];

        tr_z_zzzz_xxz[i] = 3.0 * tr_z_zz_xxz[i] * fe_0 + tr_z_zzz_xx[i] * fe_0 + ts_zzz_xxz[i] * fe_0 + tr_z_zzz_xxz[i] * pa_z[i];

        tr_z_zzzz_xyy[i] = 3.0 * tr_z_zz_xyy[i] * fe_0 + ts_zzz_xyy[i] * fe_0 + tr_z_zzz_xyy[i] * pa_z[i];

        tr_z_zzzz_xyz[i] = 3.0 * tr_z_zz_xyz[i] * fe_0 + tr_z_zzz_xy[i] * fe_0 + ts_zzz_xyz[i] * fe_0 + tr_z_zzz_xyz[i] * pa_z[i];

        tr_z_zzzz_xzz[i] = 3.0 * tr_z_zz_xzz[i] * fe_0 + 2.0 * tr_z_zzz_xz[i] * fe_0 + ts_zzz_xzz[i] * fe_0 + tr_z_zzz_xzz[i] * pa_z[i];

        tr_z_zzzz_yyy[i] = 3.0 * tr_z_zz_yyy[i] * fe_0 + ts_zzz_yyy[i] * fe_0 + tr_z_zzz_yyy[i] * pa_z[i];

        tr_z_zzzz_yyz[i] = 3.0 * tr_z_zz_yyz[i] * fe_0 + tr_z_zzz_yy[i] * fe_0 + ts_zzz_yyz[i] * fe_0 + tr_z_zzz_yyz[i] * pa_z[i];

        tr_z_zzzz_yzz[i] = 3.0 * tr_z_zz_yzz[i] * fe_0 + 2.0 * tr_z_zzz_yz[i] * fe_0 + ts_zzz_yzz[i] * fe_0 + tr_z_zzz_yzz[i] * pa_z[i];

        tr_z_zzzz_zzz[i] = 3.0 * tr_z_zz_zzz[i] * fe_0 + 3.0 * tr_z_zzz_zz[i] * fe_0 + ts_zzz_zzz[i] * fe_0 + tr_z_zzz_zzz[i] * pa_z[i];
    }
}

}  // namespace diprec
