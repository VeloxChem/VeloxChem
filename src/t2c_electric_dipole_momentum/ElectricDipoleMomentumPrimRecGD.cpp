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

#include "ElectricDipoleMomentumPrimRecGD.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_gd(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_gd,
                                      const size_t              idx_dip_dd,
                                      const size_t              idx_dip_fp,
                                      const size_t              idx_ovl_fd,
                                      const size_t              idx_dip_fd,
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

    // Set up components of auxiliary buffer : DD

    auto tr_x_xx_xx = pbuffer.data(idx_dip_dd);

    auto tr_x_xx_xy = pbuffer.data(idx_dip_dd + 1);

    auto tr_x_xx_xz = pbuffer.data(idx_dip_dd + 2);

    auto tr_x_xx_yy = pbuffer.data(idx_dip_dd + 3);

    auto tr_x_xx_yz = pbuffer.data(idx_dip_dd + 4);

    auto tr_x_xx_zz = pbuffer.data(idx_dip_dd + 5);

    auto tr_x_xy_xx = pbuffer.data(idx_dip_dd + 6);

    auto tr_x_xy_xz = pbuffer.data(idx_dip_dd + 8);

    auto tr_x_xz_xx = pbuffer.data(idx_dip_dd + 12);

    auto tr_x_xz_xy = pbuffer.data(idx_dip_dd + 13);

    auto tr_x_xz_xz = pbuffer.data(idx_dip_dd + 14);

    auto tr_x_yy_xx = pbuffer.data(idx_dip_dd + 18);

    auto tr_x_yy_xy = pbuffer.data(idx_dip_dd + 19);

    auto tr_x_yy_xz = pbuffer.data(idx_dip_dd + 20);

    auto tr_x_yy_yy = pbuffer.data(idx_dip_dd + 21);

    auto tr_x_yy_yz = pbuffer.data(idx_dip_dd + 22);

    auto tr_x_yy_zz = pbuffer.data(idx_dip_dd + 23);

    auto tr_x_yz_xz = pbuffer.data(idx_dip_dd + 26);

    auto tr_x_yz_zz = pbuffer.data(idx_dip_dd + 29);

    auto tr_x_zz_xx = pbuffer.data(idx_dip_dd + 30);

    auto tr_x_zz_xy = pbuffer.data(idx_dip_dd + 31);

    auto tr_x_zz_xz = pbuffer.data(idx_dip_dd + 32);

    auto tr_x_zz_yy = pbuffer.data(idx_dip_dd + 33);

    auto tr_x_zz_yz = pbuffer.data(idx_dip_dd + 34);

    auto tr_x_zz_zz = pbuffer.data(idx_dip_dd + 35);

    auto tr_y_xx_xx = pbuffer.data(idx_dip_dd + 36);

    auto tr_y_xx_xy = pbuffer.data(idx_dip_dd + 37);

    auto tr_y_xx_xz = pbuffer.data(idx_dip_dd + 38);

    auto tr_y_xx_yy = pbuffer.data(idx_dip_dd + 39);

    auto tr_y_xx_yz = pbuffer.data(idx_dip_dd + 40);

    auto tr_y_xx_zz = pbuffer.data(idx_dip_dd + 41);

    auto tr_y_xy_xy = pbuffer.data(idx_dip_dd + 43);

    auto tr_y_xy_yy = pbuffer.data(idx_dip_dd + 45);

    auto tr_y_xy_yz = pbuffer.data(idx_dip_dd + 46);

    auto tr_y_xy_zz = pbuffer.data(idx_dip_dd + 47);

    auto tr_y_xz_yz = pbuffer.data(idx_dip_dd + 52);

    auto tr_y_xz_zz = pbuffer.data(idx_dip_dd + 53);

    auto tr_y_yy_xx = pbuffer.data(idx_dip_dd + 54);

    auto tr_y_yy_xy = pbuffer.data(idx_dip_dd + 55);

    auto tr_y_yy_xz = pbuffer.data(idx_dip_dd + 56);

    auto tr_y_yy_yy = pbuffer.data(idx_dip_dd + 57);

    auto tr_y_yy_yz = pbuffer.data(idx_dip_dd + 58);

    auto tr_y_yy_zz = pbuffer.data(idx_dip_dd + 59);

    auto tr_y_yz_xy = pbuffer.data(idx_dip_dd + 61);

    auto tr_y_yz_yy = pbuffer.data(idx_dip_dd + 63);

    auto tr_y_yz_yz = pbuffer.data(idx_dip_dd + 64);

    auto tr_y_yz_zz = pbuffer.data(idx_dip_dd + 65);

    auto tr_y_zz_xx = pbuffer.data(idx_dip_dd + 66);

    auto tr_y_zz_xy = pbuffer.data(idx_dip_dd + 67);

    auto tr_y_zz_xz = pbuffer.data(idx_dip_dd + 68);

    auto tr_y_zz_yy = pbuffer.data(idx_dip_dd + 69);

    auto tr_y_zz_yz = pbuffer.data(idx_dip_dd + 70);

    auto tr_y_zz_zz = pbuffer.data(idx_dip_dd + 71);

    auto tr_z_xx_xx = pbuffer.data(idx_dip_dd + 72);

    auto tr_z_xx_xy = pbuffer.data(idx_dip_dd + 73);

    auto tr_z_xx_xz = pbuffer.data(idx_dip_dd + 74);

    auto tr_z_xx_yy = pbuffer.data(idx_dip_dd + 75);

    auto tr_z_xx_yz = pbuffer.data(idx_dip_dd + 76);

    auto tr_z_xx_zz = pbuffer.data(idx_dip_dd + 77);

    auto tr_z_xy_yy = pbuffer.data(idx_dip_dd + 81);

    auto tr_z_xy_yz = pbuffer.data(idx_dip_dd + 82);

    auto tr_z_xz_xz = pbuffer.data(idx_dip_dd + 86);

    auto tr_z_xz_yy = pbuffer.data(idx_dip_dd + 87);

    auto tr_z_xz_yz = pbuffer.data(idx_dip_dd + 88);

    auto tr_z_xz_zz = pbuffer.data(idx_dip_dd + 89);

    auto tr_z_yy_xx = pbuffer.data(idx_dip_dd + 90);

    auto tr_z_yy_xy = pbuffer.data(idx_dip_dd + 91);

    auto tr_z_yy_xz = pbuffer.data(idx_dip_dd + 92);

    auto tr_z_yy_yy = pbuffer.data(idx_dip_dd + 93);

    auto tr_z_yy_yz = pbuffer.data(idx_dip_dd + 94);

    auto tr_z_yy_zz = pbuffer.data(idx_dip_dd + 95);

    auto tr_z_yz_xx = pbuffer.data(idx_dip_dd + 96);

    auto tr_z_yz_xz = pbuffer.data(idx_dip_dd + 98);

    auto tr_z_yz_yy = pbuffer.data(idx_dip_dd + 99);

    auto tr_z_yz_yz = pbuffer.data(idx_dip_dd + 100);

    auto tr_z_yz_zz = pbuffer.data(idx_dip_dd + 101);

    auto tr_z_zz_xx = pbuffer.data(idx_dip_dd + 102);

    auto tr_z_zz_xy = pbuffer.data(idx_dip_dd + 103);

    auto tr_z_zz_xz = pbuffer.data(idx_dip_dd + 104);

    auto tr_z_zz_yy = pbuffer.data(idx_dip_dd + 105);

    auto tr_z_zz_yz = pbuffer.data(idx_dip_dd + 106);

    auto tr_z_zz_zz = pbuffer.data(idx_dip_dd + 107);

    // Set up components of auxiliary buffer : FP

    auto tr_x_xxx_x = pbuffer.data(idx_dip_fp);

    auto tr_x_xxx_y = pbuffer.data(idx_dip_fp + 1);

    auto tr_x_xxx_z = pbuffer.data(idx_dip_fp + 2);

    auto tr_x_xxy_x = pbuffer.data(idx_dip_fp + 3);

    auto tr_x_xxz_x = pbuffer.data(idx_dip_fp + 6);

    auto tr_x_xxz_z = pbuffer.data(idx_dip_fp + 8);

    auto tr_x_xzz_x = pbuffer.data(idx_dip_fp + 15);

    auto tr_x_yyy_x = pbuffer.data(idx_dip_fp + 18);

    auto tr_x_yyy_y = pbuffer.data(idx_dip_fp + 19);

    auto tr_x_yyy_z = pbuffer.data(idx_dip_fp + 20);

    auto tr_x_yzz_z = pbuffer.data(idx_dip_fp + 26);

    auto tr_x_zzz_x = pbuffer.data(idx_dip_fp + 27);

    auto tr_x_zzz_y = pbuffer.data(idx_dip_fp + 28);

    auto tr_x_zzz_z = pbuffer.data(idx_dip_fp + 29);

    auto tr_y_xxx_x = pbuffer.data(idx_dip_fp + 30);

    auto tr_y_xxx_y = pbuffer.data(idx_dip_fp + 31);

    auto tr_y_xxx_z = pbuffer.data(idx_dip_fp + 32);

    auto tr_y_xxy_y = pbuffer.data(idx_dip_fp + 34);

    auto tr_y_xyy_x = pbuffer.data(idx_dip_fp + 39);

    auto tr_y_xyy_y = pbuffer.data(idx_dip_fp + 40);

    auto tr_y_xyy_z = pbuffer.data(idx_dip_fp + 41);

    auto tr_y_xzz_z = pbuffer.data(idx_dip_fp + 47);

    auto tr_y_yyy_x = pbuffer.data(idx_dip_fp + 48);

    auto tr_y_yyy_y = pbuffer.data(idx_dip_fp + 49);

    auto tr_y_yyy_z = pbuffer.data(idx_dip_fp + 50);

    auto tr_y_yyz_y = pbuffer.data(idx_dip_fp + 52);

    auto tr_y_yyz_z = pbuffer.data(idx_dip_fp + 53);

    auto tr_y_yzz_x = pbuffer.data(idx_dip_fp + 54);

    auto tr_y_yzz_y = pbuffer.data(idx_dip_fp + 55);

    auto tr_y_yzz_z = pbuffer.data(idx_dip_fp + 56);

    auto tr_y_zzz_x = pbuffer.data(idx_dip_fp + 57);

    auto tr_y_zzz_y = pbuffer.data(idx_dip_fp + 58);

    auto tr_y_zzz_z = pbuffer.data(idx_dip_fp + 59);

    auto tr_z_xxx_x = pbuffer.data(idx_dip_fp + 60);

    auto tr_z_xxx_y = pbuffer.data(idx_dip_fp + 61);

    auto tr_z_xxx_z = pbuffer.data(idx_dip_fp + 62);

    auto tr_z_xxz_x = pbuffer.data(idx_dip_fp + 66);

    auto tr_z_xxz_z = pbuffer.data(idx_dip_fp + 68);

    auto tr_z_xyy_y = pbuffer.data(idx_dip_fp + 70);

    auto tr_z_xzz_x = pbuffer.data(idx_dip_fp + 75);

    auto tr_z_xzz_y = pbuffer.data(idx_dip_fp + 76);

    auto tr_z_xzz_z = pbuffer.data(idx_dip_fp + 77);

    auto tr_z_yyy_x = pbuffer.data(idx_dip_fp + 78);

    auto tr_z_yyy_y = pbuffer.data(idx_dip_fp + 79);

    auto tr_z_yyy_z = pbuffer.data(idx_dip_fp + 80);

    auto tr_z_yyz_x = pbuffer.data(idx_dip_fp + 81);

    auto tr_z_yyz_y = pbuffer.data(idx_dip_fp + 82);

    auto tr_z_yyz_z = pbuffer.data(idx_dip_fp + 83);

    auto tr_z_yzz_x = pbuffer.data(idx_dip_fp + 84);

    auto tr_z_yzz_y = pbuffer.data(idx_dip_fp + 85);

    auto tr_z_yzz_z = pbuffer.data(idx_dip_fp + 86);

    auto tr_z_zzz_x = pbuffer.data(idx_dip_fp + 87);

    auto tr_z_zzz_y = pbuffer.data(idx_dip_fp + 88);

    auto tr_z_zzz_z = pbuffer.data(idx_dip_fp + 89);

    // Set up components of auxiliary buffer : FD

    auto ts_xxx_xx = pbuffer.data(idx_ovl_fd);

    auto ts_xxx_xy = pbuffer.data(idx_ovl_fd + 1);

    auto ts_xxx_xz = pbuffer.data(idx_ovl_fd + 2);

    auto ts_xxx_yy = pbuffer.data(idx_ovl_fd + 3);

    auto ts_xxx_yz = pbuffer.data(idx_ovl_fd + 4);

    auto ts_xxx_zz = pbuffer.data(idx_ovl_fd + 5);

    auto ts_xxz_xz = pbuffer.data(idx_ovl_fd + 14);

    auto ts_xyy_yy = pbuffer.data(idx_ovl_fd + 21);

    auto ts_xyy_yz = pbuffer.data(idx_ovl_fd + 22);

    auto ts_xzz_yz = pbuffer.data(idx_ovl_fd + 34);

    auto ts_xzz_zz = pbuffer.data(idx_ovl_fd + 35);

    auto ts_yyy_xx = pbuffer.data(idx_ovl_fd + 36);

    auto ts_yyy_xy = pbuffer.data(idx_ovl_fd + 37);

    auto ts_yyy_xz = pbuffer.data(idx_ovl_fd + 38);

    auto ts_yyy_yy = pbuffer.data(idx_ovl_fd + 39);

    auto ts_yyy_yz = pbuffer.data(idx_ovl_fd + 40);

    auto ts_yyy_zz = pbuffer.data(idx_ovl_fd + 41);

    auto ts_yyz_yz = pbuffer.data(idx_ovl_fd + 46);

    auto ts_yyz_zz = pbuffer.data(idx_ovl_fd + 47);

    auto ts_yzz_xz = pbuffer.data(idx_ovl_fd + 50);

    auto ts_yzz_yy = pbuffer.data(idx_ovl_fd + 51);

    auto ts_yzz_yz = pbuffer.data(idx_ovl_fd + 52);

    auto ts_yzz_zz = pbuffer.data(idx_ovl_fd + 53);

    auto ts_zzz_xx = pbuffer.data(idx_ovl_fd + 54);

    auto ts_zzz_xy = pbuffer.data(idx_ovl_fd + 55);

    auto ts_zzz_xz = pbuffer.data(idx_ovl_fd + 56);

    auto ts_zzz_yy = pbuffer.data(idx_ovl_fd + 57);

    auto ts_zzz_yz = pbuffer.data(idx_ovl_fd + 58);

    auto ts_zzz_zz = pbuffer.data(idx_ovl_fd + 59);

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

    auto tr_x_xxy_yy = pbuffer.data(idx_dip_fd + 9);

    auto tr_x_xxy_zz = pbuffer.data(idx_dip_fd + 11);

    auto tr_x_xxz_xx = pbuffer.data(idx_dip_fd + 12);

    auto tr_x_xxz_xy = pbuffer.data(idx_dip_fd + 13);

    auto tr_x_xxz_xz = pbuffer.data(idx_dip_fd + 14);

    auto tr_x_xxz_yy = pbuffer.data(idx_dip_fd + 15);

    auto tr_x_xxz_yz = pbuffer.data(idx_dip_fd + 16);

    auto tr_x_xxz_zz = pbuffer.data(idx_dip_fd + 17);

    auto tr_x_xyy_xx = pbuffer.data(idx_dip_fd + 18);

    auto tr_x_xyy_xy = pbuffer.data(idx_dip_fd + 19);

    auto tr_x_xyy_xz = pbuffer.data(idx_dip_fd + 20);

    auto tr_x_xyy_yy = pbuffer.data(idx_dip_fd + 21);

    auto tr_x_xyy_yz = pbuffer.data(idx_dip_fd + 22);

    auto tr_x_xyz_xz = pbuffer.data(idx_dip_fd + 26);

    auto tr_x_xzz_xx = pbuffer.data(idx_dip_fd + 30);

    auto tr_x_xzz_xy = pbuffer.data(idx_dip_fd + 31);

    auto tr_x_xzz_xz = pbuffer.data(idx_dip_fd + 32);

    auto tr_x_xzz_yz = pbuffer.data(idx_dip_fd + 34);

    auto tr_x_xzz_zz = pbuffer.data(idx_dip_fd + 35);

    auto tr_x_yyy_xx = pbuffer.data(idx_dip_fd + 36);

    auto tr_x_yyy_xy = pbuffer.data(idx_dip_fd + 37);

    auto tr_x_yyy_xz = pbuffer.data(idx_dip_fd + 38);

    auto tr_x_yyy_yy = pbuffer.data(idx_dip_fd + 39);

    auto tr_x_yyy_yz = pbuffer.data(idx_dip_fd + 40);

    auto tr_x_yyy_zz = pbuffer.data(idx_dip_fd + 41);

    auto tr_x_yyz_xy = pbuffer.data(idx_dip_fd + 43);

    auto tr_x_yyz_xz = pbuffer.data(idx_dip_fd + 44);

    auto tr_x_yyz_yy = pbuffer.data(idx_dip_fd + 45);

    auto tr_x_yyz_yz = pbuffer.data(idx_dip_fd + 46);

    auto tr_x_yyz_zz = pbuffer.data(idx_dip_fd + 47);

    auto tr_x_yzz_xx = pbuffer.data(idx_dip_fd + 48);

    auto tr_x_yzz_xz = pbuffer.data(idx_dip_fd + 50);

    auto tr_x_yzz_yy = pbuffer.data(idx_dip_fd + 51);

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

    auto tr_y_xxy_xx = pbuffer.data(idx_dip_fd + 66);

    auto tr_y_xxy_xy = pbuffer.data(idx_dip_fd + 67);

    auto tr_y_xxy_yy = pbuffer.data(idx_dip_fd + 69);

    auto tr_y_xxy_yz = pbuffer.data(idx_dip_fd + 70);

    auto tr_y_xxy_zz = pbuffer.data(idx_dip_fd + 71);

    auto tr_y_xxz_xx = pbuffer.data(idx_dip_fd + 72);

    auto tr_y_xxz_xy = pbuffer.data(idx_dip_fd + 73);

    auto tr_y_xxz_xz = pbuffer.data(idx_dip_fd + 74);

    auto tr_y_xxz_yz = pbuffer.data(idx_dip_fd + 76);

    auto tr_y_xxz_zz = pbuffer.data(idx_dip_fd + 77);

    auto tr_y_xyy_xx = pbuffer.data(idx_dip_fd + 78);

    auto tr_y_xyy_xy = pbuffer.data(idx_dip_fd + 79);

    auto tr_y_xyy_xz = pbuffer.data(idx_dip_fd + 80);

    auto tr_y_xyy_yy = pbuffer.data(idx_dip_fd + 81);

    auto tr_y_xyy_yz = pbuffer.data(idx_dip_fd + 82);

    auto tr_y_xyy_zz = pbuffer.data(idx_dip_fd + 83);

    auto tr_y_xyz_yz = pbuffer.data(idx_dip_fd + 88);

    auto tr_y_xyz_zz = pbuffer.data(idx_dip_fd + 89);

    auto tr_y_xzz_xz = pbuffer.data(idx_dip_fd + 92);

    auto tr_y_xzz_yy = pbuffer.data(idx_dip_fd + 93);

    auto tr_y_xzz_yz = pbuffer.data(idx_dip_fd + 94);

    auto tr_y_xzz_zz = pbuffer.data(idx_dip_fd + 95);

    auto tr_y_yyy_xx = pbuffer.data(idx_dip_fd + 96);

    auto tr_y_yyy_xy = pbuffer.data(idx_dip_fd + 97);

    auto tr_y_yyy_xz = pbuffer.data(idx_dip_fd + 98);

    auto tr_y_yyy_yy = pbuffer.data(idx_dip_fd + 99);

    auto tr_y_yyy_yz = pbuffer.data(idx_dip_fd + 100);

    auto tr_y_yyy_zz = pbuffer.data(idx_dip_fd + 101);

    auto tr_y_yyz_xx = pbuffer.data(idx_dip_fd + 102);

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

    auto tr_z_xxy_xx = pbuffer.data(idx_dip_fd + 126);

    auto tr_z_xxy_xz = pbuffer.data(idx_dip_fd + 128);

    auto tr_z_xxy_yy = pbuffer.data(idx_dip_fd + 129);

    auto tr_z_xxy_yz = pbuffer.data(idx_dip_fd + 130);

    auto tr_z_xxz_xx = pbuffer.data(idx_dip_fd + 132);

    auto tr_z_xxz_xy = pbuffer.data(idx_dip_fd + 133);

    auto tr_z_xxz_xz = pbuffer.data(idx_dip_fd + 134);

    auto tr_z_xxz_yy = pbuffer.data(idx_dip_fd + 135);

    auto tr_z_xxz_yz = pbuffer.data(idx_dip_fd + 136);

    auto tr_z_xxz_zz = pbuffer.data(idx_dip_fd + 137);

    auto tr_z_xyy_xy = pbuffer.data(idx_dip_fd + 139);

    auto tr_z_xyy_yy = pbuffer.data(idx_dip_fd + 141);

    auto tr_z_xyy_yz = pbuffer.data(idx_dip_fd + 142);

    auto tr_z_xyy_zz = pbuffer.data(idx_dip_fd + 143);

    auto tr_z_xyz_yy = pbuffer.data(idx_dip_fd + 147);

    auto tr_z_xyz_yz = pbuffer.data(idx_dip_fd + 148);

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

    // Set up 0-6 components of targeted buffer : GD

    auto tr_x_xxxx_xx = pbuffer.data(idx_dip_gd);

    auto tr_x_xxxx_xy = pbuffer.data(idx_dip_gd + 1);

    auto tr_x_xxxx_xz = pbuffer.data(idx_dip_gd + 2);

    auto tr_x_xxxx_yy = pbuffer.data(idx_dip_gd + 3);

    auto tr_x_xxxx_yz = pbuffer.data(idx_dip_gd + 4);

    auto tr_x_xxxx_zz = pbuffer.data(idx_dip_gd + 5);

#pragma omp simd aligned(pa_x,             \
                             tr_x_xx_xx,   \
                             tr_x_xx_xy,   \
                             tr_x_xx_xz,   \
                             tr_x_xx_yy,   \
                             tr_x_xx_yz,   \
                             tr_x_xx_zz,   \
                             tr_x_xxx_x,   \
                             tr_x_xxx_xx,  \
                             tr_x_xxx_xy,  \
                             tr_x_xxx_xz,  \
                             tr_x_xxx_y,   \
                             tr_x_xxx_yy,  \
                             tr_x_xxx_yz,  \
                             tr_x_xxx_z,   \
                             tr_x_xxx_zz,  \
                             tr_x_xxxx_xx, \
                             tr_x_xxxx_xy, \
                             tr_x_xxxx_xz, \
                             tr_x_xxxx_yy, \
                             tr_x_xxxx_yz, \
                             tr_x_xxxx_zz, \
                             ts_xxx_xx,    \
                             ts_xxx_xy,    \
                             ts_xxx_xz,    \
                             ts_xxx_yy,    \
                             ts_xxx_yz,    \
                             ts_xxx_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_xx[i] = 3.0 * tr_x_xx_xx[i] * fe_0 + 2.0 * tr_x_xxx_x[i] * fe_0 + ts_xxx_xx[i] * fe_0 + tr_x_xxx_xx[i] * pa_x[i];

        tr_x_xxxx_xy[i] = 3.0 * tr_x_xx_xy[i] * fe_0 + tr_x_xxx_y[i] * fe_0 + ts_xxx_xy[i] * fe_0 + tr_x_xxx_xy[i] * pa_x[i];

        tr_x_xxxx_xz[i] = 3.0 * tr_x_xx_xz[i] * fe_0 + tr_x_xxx_z[i] * fe_0 + ts_xxx_xz[i] * fe_0 + tr_x_xxx_xz[i] * pa_x[i];

        tr_x_xxxx_yy[i] = 3.0 * tr_x_xx_yy[i] * fe_0 + ts_xxx_yy[i] * fe_0 + tr_x_xxx_yy[i] * pa_x[i];

        tr_x_xxxx_yz[i] = 3.0 * tr_x_xx_yz[i] * fe_0 + ts_xxx_yz[i] * fe_0 + tr_x_xxx_yz[i] * pa_x[i];

        tr_x_xxxx_zz[i] = 3.0 * tr_x_xx_zz[i] * fe_0 + ts_xxx_zz[i] * fe_0 + tr_x_xxx_zz[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : GD

    auto tr_x_xxxy_xx = pbuffer.data(idx_dip_gd + 6);

    auto tr_x_xxxy_xy = pbuffer.data(idx_dip_gd + 7);

    auto tr_x_xxxy_xz = pbuffer.data(idx_dip_gd + 8);

    auto tr_x_xxxy_yy = pbuffer.data(idx_dip_gd + 9);

    auto tr_x_xxxy_yz = pbuffer.data(idx_dip_gd + 10);

    auto tr_x_xxxy_zz = pbuffer.data(idx_dip_gd + 11);

#pragma omp simd aligned(pa_y,             \
                             tr_x_xxx_x,   \
                             tr_x_xxx_xx,  \
                             tr_x_xxx_xy,  \
                             tr_x_xxx_xz,  \
                             tr_x_xxx_y,   \
                             tr_x_xxx_yy,  \
                             tr_x_xxx_yz,  \
                             tr_x_xxx_z,   \
                             tr_x_xxx_zz,  \
                             tr_x_xxxy_xx, \
                             tr_x_xxxy_xy, \
                             tr_x_xxxy_xz, \
                             tr_x_xxxy_yy, \
                             tr_x_xxxy_yz, \
                             tr_x_xxxy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_xx[i] = tr_x_xxx_xx[i] * pa_y[i];

        tr_x_xxxy_xy[i] = tr_x_xxx_x[i] * fe_0 + tr_x_xxx_xy[i] * pa_y[i];

        tr_x_xxxy_xz[i] = tr_x_xxx_xz[i] * pa_y[i];

        tr_x_xxxy_yy[i] = 2.0 * tr_x_xxx_y[i] * fe_0 + tr_x_xxx_yy[i] * pa_y[i];

        tr_x_xxxy_yz[i] = tr_x_xxx_z[i] * fe_0 + tr_x_xxx_yz[i] * pa_y[i];

        tr_x_xxxy_zz[i] = tr_x_xxx_zz[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : GD

    auto tr_x_xxxz_xx = pbuffer.data(idx_dip_gd + 12);

    auto tr_x_xxxz_xy = pbuffer.data(idx_dip_gd + 13);

    auto tr_x_xxxz_xz = pbuffer.data(idx_dip_gd + 14);

    auto tr_x_xxxz_yy = pbuffer.data(idx_dip_gd + 15);

    auto tr_x_xxxz_yz = pbuffer.data(idx_dip_gd + 16);

    auto tr_x_xxxz_zz = pbuffer.data(idx_dip_gd + 17);

#pragma omp simd aligned(pa_z,             \
                             tr_x_xxx_x,   \
                             tr_x_xxx_xx,  \
                             tr_x_xxx_xy,  \
                             tr_x_xxx_xz,  \
                             tr_x_xxx_y,   \
                             tr_x_xxx_yy,  \
                             tr_x_xxx_yz,  \
                             tr_x_xxx_z,   \
                             tr_x_xxx_zz,  \
                             tr_x_xxxz_xx, \
                             tr_x_xxxz_xy, \
                             tr_x_xxxz_xz, \
                             tr_x_xxxz_yy, \
                             tr_x_xxxz_yz, \
                             tr_x_xxxz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_xx[i] = tr_x_xxx_xx[i] * pa_z[i];

        tr_x_xxxz_xy[i] = tr_x_xxx_xy[i] * pa_z[i];

        tr_x_xxxz_xz[i] = tr_x_xxx_x[i] * fe_0 + tr_x_xxx_xz[i] * pa_z[i];

        tr_x_xxxz_yy[i] = tr_x_xxx_yy[i] * pa_z[i];

        tr_x_xxxz_yz[i] = tr_x_xxx_y[i] * fe_0 + tr_x_xxx_yz[i] * pa_z[i];

        tr_x_xxxz_zz[i] = 2.0 * tr_x_xxx_z[i] * fe_0 + tr_x_xxx_zz[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : GD

    auto tr_x_xxyy_xx = pbuffer.data(idx_dip_gd + 18);

    auto tr_x_xxyy_xy = pbuffer.data(idx_dip_gd + 19);

    auto tr_x_xxyy_xz = pbuffer.data(idx_dip_gd + 20);

    auto tr_x_xxyy_yy = pbuffer.data(idx_dip_gd + 21);

    auto tr_x_xxyy_yz = pbuffer.data(idx_dip_gd + 22);

    auto tr_x_xxyy_zz = pbuffer.data(idx_dip_gd + 23);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_xx_xx,   \
                             tr_x_xx_xy,   \
                             tr_x_xx_xz,   \
                             tr_x_xx_zz,   \
                             tr_x_xxy_x,   \
                             tr_x_xxy_xx,  \
                             tr_x_xxy_xy,  \
                             tr_x_xxy_xz,  \
                             tr_x_xxy_zz,  \
                             tr_x_xxyy_xx, \
                             tr_x_xxyy_xy, \
                             tr_x_xxyy_xz, \
                             tr_x_xxyy_yy, \
                             tr_x_xxyy_yz, \
                             tr_x_xxyy_zz, \
                             tr_x_xyy_yy,  \
                             tr_x_xyy_yz,  \
                             tr_x_yy_yy,   \
                             tr_x_yy_yz,   \
                             ts_xyy_yy,    \
                             ts_xyy_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_xx[i] = tr_x_xx_xx[i] * fe_0 + tr_x_xxy_xx[i] * pa_y[i];

        tr_x_xxyy_xy[i] = tr_x_xx_xy[i] * fe_0 + tr_x_xxy_x[i] * fe_0 + tr_x_xxy_xy[i] * pa_y[i];

        tr_x_xxyy_xz[i] = tr_x_xx_xz[i] * fe_0 + tr_x_xxy_xz[i] * pa_y[i];

        tr_x_xxyy_yy[i] = tr_x_yy_yy[i] * fe_0 + ts_xyy_yy[i] * fe_0 + tr_x_xyy_yy[i] * pa_x[i];

        tr_x_xxyy_yz[i] = tr_x_yy_yz[i] * fe_0 + ts_xyy_yz[i] * fe_0 + tr_x_xyy_yz[i] * pa_x[i];

        tr_x_xxyy_zz[i] = tr_x_xx_zz[i] * fe_0 + tr_x_xxy_zz[i] * pa_y[i];
    }

    // Set up 24-30 components of targeted buffer : GD

    auto tr_x_xxyz_xx = pbuffer.data(idx_dip_gd + 24);

    auto tr_x_xxyz_xy = pbuffer.data(idx_dip_gd + 25);

    auto tr_x_xxyz_xz = pbuffer.data(idx_dip_gd + 26);

    auto tr_x_xxyz_yy = pbuffer.data(idx_dip_gd + 27);

    auto tr_x_xxyz_yz = pbuffer.data(idx_dip_gd + 28);

    auto tr_x_xxyz_zz = pbuffer.data(idx_dip_gd + 29);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_xxy_xy,  \
                             tr_x_xxy_yy,  \
                             tr_x_xxyz_xx, \
                             tr_x_xxyz_xy, \
                             tr_x_xxyz_xz, \
                             tr_x_xxyz_yy, \
                             tr_x_xxyz_yz, \
                             tr_x_xxyz_zz, \
                             tr_x_xxz_xx,  \
                             tr_x_xxz_xz,  \
                             tr_x_xxz_yz,  \
                             tr_x_xxz_z,   \
                             tr_x_xxz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyz_xx[i] = tr_x_xxz_xx[i] * pa_y[i];

        tr_x_xxyz_xy[i] = tr_x_xxy_xy[i] * pa_z[i];

        tr_x_xxyz_xz[i] = tr_x_xxz_xz[i] * pa_y[i];

        tr_x_xxyz_yy[i] = tr_x_xxy_yy[i] * pa_z[i];

        tr_x_xxyz_yz[i] = tr_x_xxz_z[i] * fe_0 + tr_x_xxz_yz[i] * pa_y[i];

        tr_x_xxyz_zz[i] = tr_x_xxz_zz[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : GD

    auto tr_x_xxzz_xx = pbuffer.data(idx_dip_gd + 30);

    auto tr_x_xxzz_xy = pbuffer.data(idx_dip_gd + 31);

    auto tr_x_xxzz_xz = pbuffer.data(idx_dip_gd + 32);

    auto tr_x_xxzz_yy = pbuffer.data(idx_dip_gd + 33);

    auto tr_x_xxzz_yz = pbuffer.data(idx_dip_gd + 34);

    auto tr_x_xxzz_zz = pbuffer.data(idx_dip_gd + 35);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_x_xx_xx,   \
                             tr_x_xx_xy,   \
                             tr_x_xx_xz,   \
                             tr_x_xx_yy,   \
                             tr_x_xxz_x,   \
                             tr_x_xxz_xx,  \
                             tr_x_xxz_xy,  \
                             tr_x_xxz_xz,  \
                             tr_x_xxz_yy,  \
                             tr_x_xxzz_xx, \
                             tr_x_xxzz_xy, \
                             tr_x_xxzz_xz, \
                             tr_x_xxzz_yy, \
                             tr_x_xxzz_yz, \
                             tr_x_xxzz_zz, \
                             tr_x_xzz_yz,  \
                             tr_x_xzz_zz,  \
                             tr_x_zz_yz,   \
                             tr_x_zz_zz,   \
                             ts_xzz_yz,    \
                             ts_xzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_xx[i] = tr_x_xx_xx[i] * fe_0 + tr_x_xxz_xx[i] * pa_z[i];

        tr_x_xxzz_xy[i] = tr_x_xx_xy[i] * fe_0 + tr_x_xxz_xy[i] * pa_z[i];

        tr_x_xxzz_xz[i] = tr_x_xx_xz[i] * fe_0 + tr_x_xxz_x[i] * fe_0 + tr_x_xxz_xz[i] * pa_z[i];

        tr_x_xxzz_yy[i] = tr_x_xx_yy[i] * fe_0 + tr_x_xxz_yy[i] * pa_z[i];

        tr_x_xxzz_yz[i] = tr_x_zz_yz[i] * fe_0 + ts_xzz_yz[i] * fe_0 + tr_x_xzz_yz[i] * pa_x[i];

        tr_x_xxzz_zz[i] = tr_x_zz_zz[i] * fe_0 + ts_xzz_zz[i] * fe_0 + tr_x_xzz_zz[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : GD

    auto tr_x_xyyy_xx = pbuffer.data(idx_dip_gd + 36);

    auto tr_x_xyyy_xy = pbuffer.data(idx_dip_gd + 37);

    auto tr_x_xyyy_xz = pbuffer.data(idx_dip_gd + 38);

    auto tr_x_xyyy_yy = pbuffer.data(idx_dip_gd + 39);

    auto tr_x_xyyy_yz = pbuffer.data(idx_dip_gd + 40);

    auto tr_x_xyyy_zz = pbuffer.data(idx_dip_gd + 41);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_xy_xx,   \
                             tr_x_xy_xz,   \
                             tr_x_xyy_xx,  \
                             tr_x_xyy_xz,  \
                             tr_x_xyyy_xx, \
                             tr_x_xyyy_xy, \
                             tr_x_xyyy_xz, \
                             tr_x_xyyy_yy, \
                             tr_x_xyyy_yz, \
                             tr_x_xyyy_zz, \
                             tr_x_yyy_xy,  \
                             tr_x_yyy_y,   \
                             tr_x_yyy_yy,  \
                             tr_x_yyy_yz,  \
                             tr_x_yyy_zz,  \
                             ts_yyy_xy,    \
                             ts_yyy_yy,    \
                             ts_yyy_yz,    \
                             ts_yyy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_xx[i] = 2.0 * tr_x_xy_xx[i] * fe_0 + tr_x_xyy_xx[i] * pa_y[i];

        tr_x_xyyy_xy[i] = tr_x_yyy_y[i] * fe_0 + ts_yyy_xy[i] * fe_0 + tr_x_yyy_xy[i] * pa_x[i];

        tr_x_xyyy_xz[i] = 2.0 * tr_x_xy_xz[i] * fe_0 + tr_x_xyy_xz[i] * pa_y[i];

        tr_x_xyyy_yy[i] = ts_yyy_yy[i] * fe_0 + tr_x_yyy_yy[i] * pa_x[i];

        tr_x_xyyy_yz[i] = ts_yyy_yz[i] * fe_0 + tr_x_yyy_yz[i] * pa_x[i];

        tr_x_xyyy_zz[i] = ts_yyy_zz[i] * fe_0 + tr_x_yyy_zz[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : GD

    auto tr_x_xyyz_xx = pbuffer.data(idx_dip_gd + 42);

    auto tr_x_xyyz_xy = pbuffer.data(idx_dip_gd + 43);

    auto tr_x_xyyz_xz = pbuffer.data(idx_dip_gd + 44);

    auto tr_x_xyyz_yy = pbuffer.data(idx_dip_gd + 45);

    auto tr_x_xyyz_yz = pbuffer.data(idx_dip_gd + 46);

    auto tr_x_xyyz_zz = pbuffer.data(idx_dip_gd + 47);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             tr_x_xyy_xx,  \
                             tr_x_xyy_xy,  \
                             tr_x_xyy_yy,  \
                             tr_x_xyyz_xx, \
                             tr_x_xyyz_xy, \
                             tr_x_xyyz_xz, \
                             tr_x_xyyz_yy, \
                             tr_x_xyyz_yz, \
                             tr_x_xyyz_zz, \
                             tr_x_xyz_xz,  \
                             tr_x_xz_xz,   \
                             tr_x_yyz_yz,  \
                             tr_x_yyz_zz,  \
                             ts_yyz_yz,    \
                             ts_yyz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_xx[i] = tr_x_xyy_xx[i] * pa_z[i];

        tr_x_xyyz_xy[i] = tr_x_xyy_xy[i] * pa_z[i];

        tr_x_xyyz_xz[i] = tr_x_xz_xz[i] * fe_0 + tr_x_xyz_xz[i] * pa_y[i];

        tr_x_xyyz_yy[i] = tr_x_xyy_yy[i] * pa_z[i];

        tr_x_xyyz_yz[i] = ts_yyz_yz[i] * fe_0 + tr_x_yyz_yz[i] * pa_x[i];

        tr_x_xyyz_zz[i] = ts_yyz_zz[i] * fe_0 + tr_x_yyz_zz[i] * pa_x[i];
    }

    // Set up 48-54 components of targeted buffer : GD

    auto tr_x_xyzz_xx = pbuffer.data(idx_dip_gd + 48);

    auto tr_x_xyzz_xy = pbuffer.data(idx_dip_gd + 49);

    auto tr_x_xyzz_xz = pbuffer.data(idx_dip_gd + 50);

    auto tr_x_xyzz_yy = pbuffer.data(idx_dip_gd + 51);

    auto tr_x_xyzz_yz = pbuffer.data(idx_dip_gd + 52);

    auto tr_x_xyzz_zz = pbuffer.data(idx_dip_gd + 53);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_x_xyzz_xx, \
                             tr_x_xyzz_xy, \
                             tr_x_xyzz_xz, \
                             tr_x_xyzz_yy, \
                             tr_x_xyzz_yz, \
                             tr_x_xyzz_zz, \
                             tr_x_xzz_x,   \
                             tr_x_xzz_xx,  \
                             tr_x_xzz_xy,  \
                             tr_x_xzz_xz,  \
                             tr_x_xzz_zz,  \
                             tr_x_yzz_yy,  \
                             tr_x_yzz_yz,  \
                             ts_yzz_yy,    \
                             ts_yzz_yz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_xx[i] = tr_x_xzz_xx[i] * pa_y[i];

        tr_x_xyzz_xy[i] = tr_x_xzz_x[i] * fe_0 + tr_x_xzz_xy[i] * pa_y[i];

        tr_x_xyzz_xz[i] = tr_x_xzz_xz[i] * pa_y[i];

        tr_x_xyzz_yy[i] = ts_yzz_yy[i] * fe_0 + tr_x_yzz_yy[i] * pa_x[i];

        tr_x_xyzz_yz[i] = ts_yzz_yz[i] * fe_0 + tr_x_yzz_yz[i] * pa_x[i];

        tr_x_xyzz_zz[i] = tr_x_xzz_zz[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : GD

    auto tr_x_xzzz_xx = pbuffer.data(idx_dip_gd + 54);

    auto tr_x_xzzz_xy = pbuffer.data(idx_dip_gd + 55);

    auto tr_x_xzzz_xz = pbuffer.data(idx_dip_gd + 56);

    auto tr_x_xzzz_yy = pbuffer.data(idx_dip_gd + 57);

    auto tr_x_xzzz_yz = pbuffer.data(idx_dip_gd + 58);

    auto tr_x_xzzz_zz = pbuffer.data(idx_dip_gd + 59);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_x_xz_xx,   \
                             tr_x_xz_xy,   \
                             tr_x_xzz_xx,  \
                             tr_x_xzz_xy,  \
                             tr_x_xzzz_xx, \
                             tr_x_xzzz_xy, \
                             tr_x_xzzz_xz, \
                             tr_x_xzzz_yy, \
                             tr_x_xzzz_yz, \
                             tr_x_xzzz_zz, \
                             tr_x_zzz_xz,  \
                             tr_x_zzz_yy,  \
                             tr_x_zzz_yz,  \
                             tr_x_zzz_z,   \
                             tr_x_zzz_zz,  \
                             ts_zzz_xz,    \
                             ts_zzz_yy,    \
                             ts_zzz_yz,    \
                             ts_zzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_xx[i] = 2.0 * tr_x_xz_xx[i] * fe_0 + tr_x_xzz_xx[i] * pa_z[i];

        tr_x_xzzz_xy[i] = 2.0 * tr_x_xz_xy[i] * fe_0 + tr_x_xzz_xy[i] * pa_z[i];

        tr_x_xzzz_xz[i] = tr_x_zzz_z[i] * fe_0 + ts_zzz_xz[i] * fe_0 + tr_x_zzz_xz[i] * pa_x[i];

        tr_x_xzzz_yy[i] = ts_zzz_yy[i] * fe_0 + tr_x_zzz_yy[i] * pa_x[i];

        tr_x_xzzz_yz[i] = ts_zzz_yz[i] * fe_0 + tr_x_zzz_yz[i] * pa_x[i];

        tr_x_xzzz_zz[i] = ts_zzz_zz[i] * fe_0 + tr_x_zzz_zz[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : GD

    auto tr_x_yyyy_xx = pbuffer.data(idx_dip_gd + 60);

    auto tr_x_yyyy_xy = pbuffer.data(idx_dip_gd + 61);

    auto tr_x_yyyy_xz = pbuffer.data(idx_dip_gd + 62);

    auto tr_x_yyyy_yy = pbuffer.data(idx_dip_gd + 63);

    auto tr_x_yyyy_yz = pbuffer.data(idx_dip_gd + 64);

    auto tr_x_yyyy_zz = pbuffer.data(idx_dip_gd + 65);

#pragma omp simd aligned(pa_y,             \
                             tr_x_yy_xx,   \
                             tr_x_yy_xy,   \
                             tr_x_yy_xz,   \
                             tr_x_yy_yy,   \
                             tr_x_yy_yz,   \
                             tr_x_yy_zz,   \
                             tr_x_yyy_x,   \
                             tr_x_yyy_xx,  \
                             tr_x_yyy_xy,  \
                             tr_x_yyy_xz,  \
                             tr_x_yyy_y,   \
                             tr_x_yyy_yy,  \
                             tr_x_yyy_yz,  \
                             tr_x_yyy_z,   \
                             tr_x_yyy_zz,  \
                             tr_x_yyyy_xx, \
                             tr_x_yyyy_xy, \
                             tr_x_yyyy_xz, \
                             tr_x_yyyy_yy, \
                             tr_x_yyyy_yz, \
                             tr_x_yyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_xx[i] = 3.0 * tr_x_yy_xx[i] * fe_0 + tr_x_yyy_xx[i] * pa_y[i];

        tr_x_yyyy_xy[i] = 3.0 * tr_x_yy_xy[i] * fe_0 + tr_x_yyy_x[i] * fe_0 + tr_x_yyy_xy[i] * pa_y[i];

        tr_x_yyyy_xz[i] = 3.0 * tr_x_yy_xz[i] * fe_0 + tr_x_yyy_xz[i] * pa_y[i];

        tr_x_yyyy_yy[i] = 3.0 * tr_x_yy_yy[i] * fe_0 + 2.0 * tr_x_yyy_y[i] * fe_0 + tr_x_yyy_yy[i] * pa_y[i];

        tr_x_yyyy_yz[i] = 3.0 * tr_x_yy_yz[i] * fe_0 + tr_x_yyy_z[i] * fe_0 + tr_x_yyy_yz[i] * pa_y[i];

        tr_x_yyyy_zz[i] = 3.0 * tr_x_yy_zz[i] * fe_0 + tr_x_yyy_zz[i] * pa_y[i];
    }

    // Set up 66-72 components of targeted buffer : GD

    auto tr_x_yyyz_xx = pbuffer.data(idx_dip_gd + 66);

    auto tr_x_yyyz_xy = pbuffer.data(idx_dip_gd + 67);

    auto tr_x_yyyz_xz = pbuffer.data(idx_dip_gd + 68);

    auto tr_x_yyyz_yy = pbuffer.data(idx_dip_gd + 69);

    auto tr_x_yyyz_yz = pbuffer.data(idx_dip_gd + 70);

    auto tr_x_yyyz_zz = pbuffer.data(idx_dip_gd + 71);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_yyy_xx,  \
                             tr_x_yyy_xy,  \
                             tr_x_yyy_y,   \
                             tr_x_yyy_yy,  \
                             tr_x_yyy_yz,  \
                             tr_x_yyyz_xx, \
                             tr_x_yyyz_xy, \
                             tr_x_yyyz_xz, \
                             tr_x_yyyz_yy, \
                             tr_x_yyyz_yz, \
                             tr_x_yyyz_zz, \
                             tr_x_yyz_xz,  \
                             tr_x_yyz_zz,  \
                             tr_x_yz_xz,   \
                             tr_x_yz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_xx[i] = tr_x_yyy_xx[i] * pa_z[i];

        tr_x_yyyz_xy[i] = tr_x_yyy_xy[i] * pa_z[i];

        tr_x_yyyz_xz[i] = 2.0 * tr_x_yz_xz[i] * fe_0 + tr_x_yyz_xz[i] * pa_y[i];

        tr_x_yyyz_yy[i] = tr_x_yyy_yy[i] * pa_z[i];

        tr_x_yyyz_yz[i] = tr_x_yyy_y[i] * fe_0 + tr_x_yyy_yz[i] * pa_z[i];

        tr_x_yyyz_zz[i] = 2.0 * tr_x_yz_zz[i] * fe_0 + tr_x_yyz_zz[i] * pa_y[i];
    }

    // Set up 72-78 components of targeted buffer : GD

    auto tr_x_yyzz_xx = pbuffer.data(idx_dip_gd + 72);

    auto tr_x_yyzz_xy = pbuffer.data(idx_dip_gd + 73);

    auto tr_x_yyzz_xz = pbuffer.data(idx_dip_gd + 74);

    auto tr_x_yyzz_yy = pbuffer.data(idx_dip_gd + 75);

    auto tr_x_yyzz_yz = pbuffer.data(idx_dip_gd + 76);

    auto tr_x_yyzz_zz = pbuffer.data(idx_dip_gd + 77);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_x_yy_xy,   \
                             tr_x_yy_yy,   \
                             tr_x_yyz_xy,  \
                             tr_x_yyz_yy,  \
                             tr_x_yyzz_xx, \
                             tr_x_yyzz_xy, \
                             tr_x_yyzz_xz, \
                             tr_x_yyzz_yy, \
                             tr_x_yyzz_yz, \
                             tr_x_yyzz_zz, \
                             tr_x_yzz_xx,  \
                             tr_x_yzz_xz,  \
                             tr_x_yzz_yz,  \
                             tr_x_yzz_z,   \
                             tr_x_yzz_zz,  \
                             tr_x_zz_xx,   \
                             tr_x_zz_xz,   \
                             tr_x_zz_yz,   \
                             tr_x_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_xx[i] = tr_x_zz_xx[i] * fe_0 + tr_x_yzz_xx[i] * pa_y[i];

        tr_x_yyzz_xy[i] = tr_x_yy_xy[i] * fe_0 + tr_x_yyz_xy[i] * pa_z[i];

        tr_x_yyzz_xz[i] = tr_x_zz_xz[i] * fe_0 + tr_x_yzz_xz[i] * pa_y[i];

        tr_x_yyzz_yy[i] = tr_x_yy_yy[i] * fe_0 + tr_x_yyz_yy[i] * pa_z[i];

        tr_x_yyzz_yz[i] = tr_x_zz_yz[i] * fe_0 + tr_x_yzz_z[i] * fe_0 + tr_x_yzz_yz[i] * pa_y[i];

        tr_x_yyzz_zz[i] = tr_x_zz_zz[i] * fe_0 + tr_x_yzz_zz[i] * pa_y[i];
    }

    // Set up 78-84 components of targeted buffer : GD

    auto tr_x_yzzz_xx = pbuffer.data(idx_dip_gd + 78);

    auto tr_x_yzzz_xy = pbuffer.data(idx_dip_gd + 79);

    auto tr_x_yzzz_xz = pbuffer.data(idx_dip_gd + 80);

    auto tr_x_yzzz_yy = pbuffer.data(idx_dip_gd + 81);

    auto tr_x_yzzz_yz = pbuffer.data(idx_dip_gd + 82);

    auto tr_x_yzzz_zz = pbuffer.data(idx_dip_gd + 83);

#pragma omp simd aligned(pa_y,             \
                             tr_x_yzzz_xx, \
                             tr_x_yzzz_xy, \
                             tr_x_yzzz_xz, \
                             tr_x_yzzz_yy, \
                             tr_x_yzzz_yz, \
                             tr_x_yzzz_zz, \
                             tr_x_zzz_x,   \
                             tr_x_zzz_xx,  \
                             tr_x_zzz_xy,  \
                             tr_x_zzz_xz,  \
                             tr_x_zzz_y,   \
                             tr_x_zzz_yy,  \
                             tr_x_zzz_yz,  \
                             tr_x_zzz_z,   \
                             tr_x_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_xx[i] = tr_x_zzz_xx[i] * pa_y[i];

        tr_x_yzzz_xy[i] = tr_x_zzz_x[i] * fe_0 + tr_x_zzz_xy[i] * pa_y[i];

        tr_x_yzzz_xz[i] = tr_x_zzz_xz[i] * pa_y[i];

        tr_x_yzzz_yy[i] = 2.0 * tr_x_zzz_y[i] * fe_0 + tr_x_zzz_yy[i] * pa_y[i];

        tr_x_yzzz_yz[i] = tr_x_zzz_z[i] * fe_0 + tr_x_zzz_yz[i] * pa_y[i];

        tr_x_yzzz_zz[i] = tr_x_zzz_zz[i] * pa_y[i];
    }

    // Set up 84-90 components of targeted buffer : GD

    auto tr_x_zzzz_xx = pbuffer.data(idx_dip_gd + 84);

    auto tr_x_zzzz_xy = pbuffer.data(idx_dip_gd + 85);

    auto tr_x_zzzz_xz = pbuffer.data(idx_dip_gd + 86);

    auto tr_x_zzzz_yy = pbuffer.data(idx_dip_gd + 87);

    auto tr_x_zzzz_yz = pbuffer.data(idx_dip_gd + 88);

    auto tr_x_zzzz_zz = pbuffer.data(idx_dip_gd + 89);

#pragma omp simd aligned(pa_z,             \
                             tr_x_zz_xx,   \
                             tr_x_zz_xy,   \
                             tr_x_zz_xz,   \
                             tr_x_zz_yy,   \
                             tr_x_zz_yz,   \
                             tr_x_zz_zz,   \
                             tr_x_zzz_x,   \
                             tr_x_zzz_xx,  \
                             tr_x_zzz_xy,  \
                             tr_x_zzz_xz,  \
                             tr_x_zzz_y,   \
                             tr_x_zzz_yy,  \
                             tr_x_zzz_yz,  \
                             tr_x_zzz_z,   \
                             tr_x_zzz_zz,  \
                             tr_x_zzzz_xx, \
                             tr_x_zzzz_xy, \
                             tr_x_zzzz_xz, \
                             tr_x_zzzz_yy, \
                             tr_x_zzzz_yz, \
                             tr_x_zzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_xx[i] = 3.0 * tr_x_zz_xx[i] * fe_0 + tr_x_zzz_xx[i] * pa_z[i];

        tr_x_zzzz_xy[i] = 3.0 * tr_x_zz_xy[i] * fe_0 + tr_x_zzz_xy[i] * pa_z[i];

        tr_x_zzzz_xz[i] = 3.0 * tr_x_zz_xz[i] * fe_0 + tr_x_zzz_x[i] * fe_0 + tr_x_zzz_xz[i] * pa_z[i];

        tr_x_zzzz_yy[i] = 3.0 * tr_x_zz_yy[i] * fe_0 + tr_x_zzz_yy[i] * pa_z[i];

        tr_x_zzzz_yz[i] = 3.0 * tr_x_zz_yz[i] * fe_0 + tr_x_zzz_y[i] * fe_0 + tr_x_zzz_yz[i] * pa_z[i];

        tr_x_zzzz_zz[i] = 3.0 * tr_x_zz_zz[i] * fe_0 + 2.0 * tr_x_zzz_z[i] * fe_0 + tr_x_zzz_zz[i] * pa_z[i];
    }

    // Set up 90-96 components of targeted buffer : GD

    auto tr_y_xxxx_xx = pbuffer.data(idx_dip_gd + 90);

    auto tr_y_xxxx_xy = pbuffer.data(idx_dip_gd + 91);

    auto tr_y_xxxx_xz = pbuffer.data(idx_dip_gd + 92);

    auto tr_y_xxxx_yy = pbuffer.data(idx_dip_gd + 93);

    auto tr_y_xxxx_yz = pbuffer.data(idx_dip_gd + 94);

    auto tr_y_xxxx_zz = pbuffer.data(idx_dip_gd + 95);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xx_xx,   \
                             tr_y_xx_xy,   \
                             tr_y_xx_xz,   \
                             tr_y_xx_yy,   \
                             tr_y_xx_yz,   \
                             tr_y_xx_zz,   \
                             tr_y_xxx_x,   \
                             tr_y_xxx_xx,  \
                             tr_y_xxx_xy,  \
                             tr_y_xxx_xz,  \
                             tr_y_xxx_y,   \
                             tr_y_xxx_yy,  \
                             tr_y_xxx_yz,  \
                             tr_y_xxx_z,   \
                             tr_y_xxx_zz,  \
                             tr_y_xxxx_xx, \
                             tr_y_xxxx_xy, \
                             tr_y_xxxx_xz, \
                             tr_y_xxxx_yy, \
                             tr_y_xxxx_yz, \
                             tr_y_xxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_xx[i] = 3.0 * tr_y_xx_xx[i] * fe_0 + 2.0 * tr_y_xxx_x[i] * fe_0 + tr_y_xxx_xx[i] * pa_x[i];

        tr_y_xxxx_xy[i] = 3.0 * tr_y_xx_xy[i] * fe_0 + tr_y_xxx_y[i] * fe_0 + tr_y_xxx_xy[i] * pa_x[i];

        tr_y_xxxx_xz[i] = 3.0 * tr_y_xx_xz[i] * fe_0 + tr_y_xxx_z[i] * fe_0 + tr_y_xxx_xz[i] * pa_x[i];

        tr_y_xxxx_yy[i] = 3.0 * tr_y_xx_yy[i] * fe_0 + tr_y_xxx_yy[i] * pa_x[i];

        tr_y_xxxx_yz[i] = 3.0 * tr_y_xx_yz[i] * fe_0 + tr_y_xxx_yz[i] * pa_x[i];

        tr_y_xxxx_zz[i] = 3.0 * tr_y_xx_zz[i] * fe_0 + tr_y_xxx_zz[i] * pa_x[i];
    }

    // Set up 96-102 components of targeted buffer : GD

    auto tr_y_xxxy_xx = pbuffer.data(idx_dip_gd + 96);

    auto tr_y_xxxy_xy = pbuffer.data(idx_dip_gd + 97);

    auto tr_y_xxxy_xz = pbuffer.data(idx_dip_gd + 98);

    auto tr_y_xxxy_yy = pbuffer.data(idx_dip_gd + 99);

    auto tr_y_xxxy_yz = pbuffer.data(idx_dip_gd + 100);

    auto tr_y_xxxy_zz = pbuffer.data(idx_dip_gd + 101);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_y_xxx_xx,  \
                             tr_y_xxx_xz,  \
                             tr_y_xxxy_xx, \
                             tr_y_xxxy_xy, \
                             tr_y_xxxy_xz, \
                             tr_y_xxxy_yy, \
                             tr_y_xxxy_yz, \
                             tr_y_xxxy_zz, \
                             tr_y_xxy_xy,  \
                             tr_y_xxy_y,   \
                             tr_y_xxy_yy,  \
                             tr_y_xxy_yz,  \
                             tr_y_xxy_zz,  \
                             tr_y_xy_xy,   \
                             tr_y_xy_yy,   \
                             tr_y_xy_yz,   \
                             tr_y_xy_zz,   \
                             ts_xxx_xx,    \
                             ts_xxx_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_xx[i] = ts_xxx_xx[i] * fe_0 + tr_y_xxx_xx[i] * pa_y[i];

        tr_y_xxxy_xy[i] = 2.0 * tr_y_xy_xy[i] * fe_0 + tr_y_xxy_y[i] * fe_0 + tr_y_xxy_xy[i] * pa_x[i];

        tr_y_xxxy_xz[i] = ts_xxx_xz[i] * fe_0 + tr_y_xxx_xz[i] * pa_y[i];

        tr_y_xxxy_yy[i] = 2.0 * tr_y_xy_yy[i] * fe_0 + tr_y_xxy_yy[i] * pa_x[i];

        tr_y_xxxy_yz[i] = 2.0 * tr_y_xy_yz[i] * fe_0 + tr_y_xxy_yz[i] * pa_x[i];

        tr_y_xxxy_zz[i] = 2.0 * tr_y_xy_zz[i] * fe_0 + tr_y_xxy_zz[i] * pa_x[i];
    }

    // Set up 102-108 components of targeted buffer : GD

    auto tr_y_xxxz_xx = pbuffer.data(idx_dip_gd + 102);

    auto tr_y_xxxz_xy = pbuffer.data(idx_dip_gd + 103);

    auto tr_y_xxxz_xz = pbuffer.data(idx_dip_gd + 104);

    auto tr_y_xxxz_yy = pbuffer.data(idx_dip_gd + 105);

    auto tr_y_xxxz_yz = pbuffer.data(idx_dip_gd + 106);

    auto tr_y_xxxz_zz = pbuffer.data(idx_dip_gd + 107);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_xxx_x,   \
                             tr_y_xxx_xx,  \
                             tr_y_xxx_xy,  \
                             tr_y_xxx_xz,  \
                             tr_y_xxx_yy,  \
                             tr_y_xxxz_xx, \
                             tr_y_xxxz_xy, \
                             tr_y_xxxz_xz, \
                             tr_y_xxxz_yy, \
                             tr_y_xxxz_yz, \
                             tr_y_xxxz_zz, \
                             tr_y_xxz_yz,  \
                             tr_y_xxz_zz,  \
                             tr_y_xz_yz,   \
                             tr_y_xz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_xx[i] = tr_y_xxx_xx[i] * pa_z[i];

        tr_y_xxxz_xy[i] = tr_y_xxx_xy[i] * pa_z[i];

        tr_y_xxxz_xz[i] = tr_y_xxx_x[i] * fe_0 + tr_y_xxx_xz[i] * pa_z[i];

        tr_y_xxxz_yy[i] = tr_y_xxx_yy[i] * pa_z[i];

        tr_y_xxxz_yz[i] = 2.0 * tr_y_xz_yz[i] * fe_0 + tr_y_xxz_yz[i] * pa_x[i];

        tr_y_xxxz_zz[i] = 2.0 * tr_y_xz_zz[i] * fe_0 + tr_y_xxz_zz[i] * pa_x[i];
    }

    // Set up 108-114 components of targeted buffer : GD

    auto tr_y_xxyy_xx = pbuffer.data(idx_dip_gd + 108);

    auto tr_y_xxyy_xy = pbuffer.data(idx_dip_gd + 109);

    auto tr_y_xxyy_xz = pbuffer.data(idx_dip_gd + 110);

    auto tr_y_xxyy_yy = pbuffer.data(idx_dip_gd + 111);

    auto tr_y_xxyy_yz = pbuffer.data(idx_dip_gd + 112);

    auto tr_y_xxyy_zz = pbuffer.data(idx_dip_gd + 113);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xxyy_xx, \
                             tr_y_xxyy_xy, \
                             tr_y_xxyy_xz, \
                             tr_y_xxyy_yy, \
                             tr_y_xxyy_yz, \
                             tr_y_xxyy_zz, \
                             tr_y_xyy_x,   \
                             tr_y_xyy_xx,  \
                             tr_y_xyy_xy,  \
                             tr_y_xyy_xz,  \
                             tr_y_xyy_y,   \
                             tr_y_xyy_yy,  \
                             tr_y_xyy_yz,  \
                             tr_y_xyy_z,   \
                             tr_y_xyy_zz,  \
                             tr_y_yy_xx,   \
                             tr_y_yy_xy,   \
                             tr_y_yy_xz,   \
                             tr_y_yy_yy,   \
                             tr_y_yy_yz,   \
                             tr_y_yy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_xx[i] = tr_y_yy_xx[i] * fe_0 + 2.0 * tr_y_xyy_x[i] * fe_0 + tr_y_xyy_xx[i] * pa_x[i];

        tr_y_xxyy_xy[i] = tr_y_yy_xy[i] * fe_0 + tr_y_xyy_y[i] * fe_0 + tr_y_xyy_xy[i] * pa_x[i];

        tr_y_xxyy_xz[i] = tr_y_yy_xz[i] * fe_0 + tr_y_xyy_z[i] * fe_0 + tr_y_xyy_xz[i] * pa_x[i];

        tr_y_xxyy_yy[i] = tr_y_yy_yy[i] * fe_0 + tr_y_xyy_yy[i] * pa_x[i];

        tr_y_xxyy_yz[i] = tr_y_yy_yz[i] * fe_0 + tr_y_xyy_yz[i] * pa_x[i];

        tr_y_xxyy_zz[i] = tr_y_yy_zz[i] * fe_0 + tr_y_xyy_zz[i] * pa_x[i];
    }

    // Set up 114-120 components of targeted buffer : GD

    auto tr_y_xxyz_xx = pbuffer.data(idx_dip_gd + 114);

    auto tr_y_xxyz_xy = pbuffer.data(idx_dip_gd + 115);

    auto tr_y_xxyz_xz = pbuffer.data(idx_dip_gd + 116);

    auto tr_y_xxyz_yy = pbuffer.data(idx_dip_gd + 117);

    auto tr_y_xxyz_yz = pbuffer.data(idx_dip_gd + 118);

    auto tr_y_xxyz_zz = pbuffer.data(idx_dip_gd + 119);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             pa_z,         \
                             tr_y_xxy_xx,  \
                             tr_y_xxy_xy,  \
                             tr_y_xxy_yy,  \
                             tr_y_xxyz_xx, \
                             tr_y_xxyz_xy, \
                             tr_y_xxyz_xz, \
                             tr_y_xxyz_yy, \
                             tr_y_xxyz_yz, \
                             tr_y_xxyz_zz, \
                             tr_y_xxz_xz,  \
                             tr_y_xyz_yz,  \
                             tr_y_xyz_zz,  \
                             tr_y_yz_yz,   \
                             tr_y_yz_zz,   \
                             ts_xxz_xz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_xx[i] = tr_y_xxy_xx[i] * pa_z[i];

        tr_y_xxyz_xy[i] = tr_y_xxy_xy[i] * pa_z[i];

        tr_y_xxyz_xz[i] = ts_xxz_xz[i] * fe_0 + tr_y_xxz_xz[i] * pa_y[i];

        tr_y_xxyz_yy[i] = tr_y_xxy_yy[i] * pa_z[i];

        tr_y_xxyz_yz[i] = tr_y_yz_yz[i] * fe_0 + tr_y_xyz_yz[i] * pa_x[i];

        tr_y_xxyz_zz[i] = tr_y_yz_zz[i] * fe_0 + tr_y_xyz_zz[i] * pa_x[i];
    }

    // Set up 120-126 components of targeted buffer : GD

    auto tr_y_xxzz_xx = pbuffer.data(idx_dip_gd + 120);

    auto tr_y_xxzz_xy = pbuffer.data(idx_dip_gd + 121);

    auto tr_y_xxzz_xz = pbuffer.data(idx_dip_gd + 122);

    auto tr_y_xxzz_yy = pbuffer.data(idx_dip_gd + 123);

    auto tr_y_xxzz_yz = pbuffer.data(idx_dip_gd + 124);

    auto tr_y_xxzz_zz = pbuffer.data(idx_dip_gd + 125);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_xx_xx,   \
                             tr_y_xx_xy,   \
                             tr_y_xxz_xx,  \
                             tr_y_xxz_xy,  \
                             tr_y_xxzz_xx, \
                             tr_y_xxzz_xy, \
                             tr_y_xxzz_xz, \
                             tr_y_xxzz_yy, \
                             tr_y_xxzz_yz, \
                             tr_y_xxzz_zz, \
                             tr_y_xzz_xz,  \
                             tr_y_xzz_yy,  \
                             tr_y_xzz_yz,  \
                             tr_y_xzz_z,   \
                             tr_y_xzz_zz,  \
                             tr_y_zz_xz,   \
                             tr_y_zz_yy,   \
                             tr_y_zz_yz,   \
                             tr_y_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_xx[i] = tr_y_xx_xx[i] * fe_0 + tr_y_xxz_xx[i] * pa_z[i];

        tr_y_xxzz_xy[i] = tr_y_xx_xy[i] * fe_0 + tr_y_xxz_xy[i] * pa_z[i];

        tr_y_xxzz_xz[i] = tr_y_zz_xz[i] * fe_0 + tr_y_xzz_z[i] * fe_0 + tr_y_xzz_xz[i] * pa_x[i];

        tr_y_xxzz_yy[i] = tr_y_zz_yy[i] * fe_0 + tr_y_xzz_yy[i] * pa_x[i];

        tr_y_xxzz_yz[i] = tr_y_zz_yz[i] * fe_0 + tr_y_xzz_yz[i] * pa_x[i];

        tr_y_xxzz_zz[i] = tr_y_zz_zz[i] * fe_0 + tr_y_xzz_zz[i] * pa_x[i];
    }

    // Set up 126-132 components of targeted buffer : GD

    auto tr_y_xyyy_xx = pbuffer.data(idx_dip_gd + 126);

    auto tr_y_xyyy_xy = pbuffer.data(idx_dip_gd + 127);

    auto tr_y_xyyy_xz = pbuffer.data(idx_dip_gd + 128);

    auto tr_y_xyyy_yy = pbuffer.data(idx_dip_gd + 129);

    auto tr_y_xyyy_yz = pbuffer.data(idx_dip_gd + 130);

    auto tr_y_xyyy_zz = pbuffer.data(idx_dip_gd + 131);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xyyy_xx, \
                             tr_y_xyyy_xy, \
                             tr_y_xyyy_xz, \
                             tr_y_xyyy_yy, \
                             tr_y_xyyy_yz, \
                             tr_y_xyyy_zz, \
                             tr_y_yyy_x,   \
                             tr_y_yyy_xx,  \
                             tr_y_yyy_xy,  \
                             tr_y_yyy_xz,  \
                             tr_y_yyy_y,   \
                             tr_y_yyy_yy,  \
                             tr_y_yyy_yz,  \
                             tr_y_yyy_z,   \
                             tr_y_yyy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_xx[i] = 2.0 * tr_y_yyy_x[i] * fe_0 + tr_y_yyy_xx[i] * pa_x[i];

        tr_y_xyyy_xy[i] = tr_y_yyy_y[i] * fe_0 + tr_y_yyy_xy[i] * pa_x[i];

        tr_y_xyyy_xz[i] = tr_y_yyy_z[i] * fe_0 + tr_y_yyy_xz[i] * pa_x[i];

        tr_y_xyyy_yy[i] = tr_y_yyy_yy[i] * pa_x[i];

        tr_y_xyyy_yz[i] = tr_y_yyy_yz[i] * pa_x[i];

        tr_y_xyyy_zz[i] = tr_y_yyy_zz[i] * pa_x[i];
    }

    // Set up 132-138 components of targeted buffer : GD

    auto tr_y_xyyz_xx = pbuffer.data(idx_dip_gd + 132);

    auto tr_y_xyyz_xy = pbuffer.data(idx_dip_gd + 133);

    auto tr_y_xyyz_xz = pbuffer.data(idx_dip_gd + 134);

    auto tr_y_xyyz_yy = pbuffer.data(idx_dip_gd + 135);

    auto tr_y_xyyz_yz = pbuffer.data(idx_dip_gd + 136);

    auto tr_y_xyyz_zz = pbuffer.data(idx_dip_gd + 137);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_y_xyy_xx,  \
                             tr_y_xyy_xy,  \
                             tr_y_xyyz_xx, \
                             tr_y_xyyz_xy, \
                             tr_y_xyyz_xz, \
                             tr_y_xyyz_yy, \
                             tr_y_xyyz_yz, \
                             tr_y_xyyz_zz, \
                             tr_y_yyz_xz,  \
                             tr_y_yyz_yy,  \
                             tr_y_yyz_yz,  \
                             tr_y_yyz_z,   \
                             tr_y_yyz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyz_xx[i] = tr_y_xyy_xx[i] * pa_z[i];

        tr_y_xyyz_xy[i] = tr_y_xyy_xy[i] * pa_z[i];

        tr_y_xyyz_xz[i] = tr_y_yyz_z[i] * fe_0 + tr_y_yyz_xz[i] * pa_x[i];

        tr_y_xyyz_yy[i] = tr_y_yyz_yy[i] * pa_x[i];

        tr_y_xyyz_yz[i] = tr_y_yyz_yz[i] * pa_x[i];

        tr_y_xyyz_zz[i] = tr_y_yyz_zz[i] * pa_x[i];
    }

    // Set up 138-144 components of targeted buffer : GD

    auto tr_y_xyzz_xx = pbuffer.data(idx_dip_gd + 138);

    auto tr_y_xyzz_xy = pbuffer.data(idx_dip_gd + 139);

    auto tr_y_xyzz_xz = pbuffer.data(idx_dip_gd + 140);

    auto tr_y_xyzz_yy = pbuffer.data(idx_dip_gd + 141);

    auto tr_y_xyzz_yz = pbuffer.data(idx_dip_gd + 142);

    auto tr_y_xyzz_zz = pbuffer.data(idx_dip_gd + 143);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xyzz_xx, \
                             tr_y_xyzz_xy, \
                             tr_y_xyzz_xz, \
                             tr_y_xyzz_yy, \
                             tr_y_xyzz_yz, \
                             tr_y_xyzz_zz, \
                             tr_y_yzz_x,   \
                             tr_y_yzz_xx,  \
                             tr_y_yzz_xy,  \
                             tr_y_yzz_xz,  \
                             tr_y_yzz_y,   \
                             tr_y_yzz_yy,  \
                             tr_y_yzz_yz,  \
                             tr_y_yzz_z,   \
                             tr_y_yzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_xx[i] = 2.0 * tr_y_yzz_x[i] * fe_0 + tr_y_yzz_xx[i] * pa_x[i];

        tr_y_xyzz_xy[i] = tr_y_yzz_y[i] * fe_0 + tr_y_yzz_xy[i] * pa_x[i];

        tr_y_xyzz_xz[i] = tr_y_yzz_z[i] * fe_0 + tr_y_yzz_xz[i] * pa_x[i];

        tr_y_xyzz_yy[i] = tr_y_yzz_yy[i] * pa_x[i];

        tr_y_xyzz_yz[i] = tr_y_yzz_yz[i] * pa_x[i];

        tr_y_xyzz_zz[i] = tr_y_yzz_zz[i] * pa_x[i];
    }

    // Set up 144-150 components of targeted buffer : GD

    auto tr_y_xzzz_xx = pbuffer.data(idx_dip_gd + 144);

    auto tr_y_xzzz_xy = pbuffer.data(idx_dip_gd + 145);

    auto tr_y_xzzz_xz = pbuffer.data(idx_dip_gd + 146);

    auto tr_y_xzzz_yy = pbuffer.data(idx_dip_gd + 147);

    auto tr_y_xzzz_yz = pbuffer.data(idx_dip_gd + 148);

    auto tr_y_xzzz_zz = pbuffer.data(idx_dip_gd + 149);

#pragma omp simd aligned(pa_x,             \
                             tr_y_xzzz_xx, \
                             tr_y_xzzz_xy, \
                             tr_y_xzzz_xz, \
                             tr_y_xzzz_yy, \
                             tr_y_xzzz_yz, \
                             tr_y_xzzz_zz, \
                             tr_y_zzz_x,   \
                             tr_y_zzz_xx,  \
                             tr_y_zzz_xy,  \
                             tr_y_zzz_xz,  \
                             tr_y_zzz_y,   \
                             tr_y_zzz_yy,  \
                             tr_y_zzz_yz,  \
                             tr_y_zzz_z,   \
                             tr_y_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_xx[i] = 2.0 * tr_y_zzz_x[i] * fe_0 + tr_y_zzz_xx[i] * pa_x[i];

        tr_y_xzzz_xy[i] = tr_y_zzz_y[i] * fe_0 + tr_y_zzz_xy[i] * pa_x[i];

        tr_y_xzzz_xz[i] = tr_y_zzz_z[i] * fe_0 + tr_y_zzz_xz[i] * pa_x[i];

        tr_y_xzzz_yy[i] = tr_y_zzz_yy[i] * pa_x[i];

        tr_y_xzzz_yz[i] = tr_y_zzz_yz[i] * pa_x[i];

        tr_y_xzzz_zz[i] = tr_y_zzz_zz[i] * pa_x[i];
    }

    // Set up 150-156 components of targeted buffer : GD

    auto tr_y_yyyy_xx = pbuffer.data(idx_dip_gd + 150);

    auto tr_y_yyyy_xy = pbuffer.data(idx_dip_gd + 151);

    auto tr_y_yyyy_xz = pbuffer.data(idx_dip_gd + 152);

    auto tr_y_yyyy_yy = pbuffer.data(idx_dip_gd + 153);

    auto tr_y_yyyy_yz = pbuffer.data(idx_dip_gd + 154);

    auto tr_y_yyyy_zz = pbuffer.data(idx_dip_gd + 155);

#pragma omp simd aligned(pa_y,             \
                             tr_y_yy_xx,   \
                             tr_y_yy_xy,   \
                             tr_y_yy_xz,   \
                             tr_y_yy_yy,   \
                             tr_y_yy_yz,   \
                             tr_y_yy_zz,   \
                             tr_y_yyy_x,   \
                             tr_y_yyy_xx,  \
                             tr_y_yyy_xy,  \
                             tr_y_yyy_xz,  \
                             tr_y_yyy_y,   \
                             tr_y_yyy_yy,  \
                             tr_y_yyy_yz,  \
                             tr_y_yyy_z,   \
                             tr_y_yyy_zz,  \
                             tr_y_yyyy_xx, \
                             tr_y_yyyy_xy, \
                             tr_y_yyyy_xz, \
                             tr_y_yyyy_yy, \
                             tr_y_yyyy_yz, \
                             tr_y_yyyy_zz, \
                             ts_yyy_xx,    \
                             ts_yyy_xy,    \
                             ts_yyy_xz,    \
                             ts_yyy_yy,    \
                             ts_yyy_yz,    \
                             ts_yyy_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_xx[i] = 3.0 * tr_y_yy_xx[i] * fe_0 + ts_yyy_xx[i] * fe_0 + tr_y_yyy_xx[i] * pa_y[i];

        tr_y_yyyy_xy[i] = 3.0 * tr_y_yy_xy[i] * fe_0 + tr_y_yyy_x[i] * fe_0 + ts_yyy_xy[i] * fe_0 + tr_y_yyy_xy[i] * pa_y[i];

        tr_y_yyyy_xz[i] = 3.0 * tr_y_yy_xz[i] * fe_0 + ts_yyy_xz[i] * fe_0 + tr_y_yyy_xz[i] * pa_y[i];

        tr_y_yyyy_yy[i] = 3.0 * tr_y_yy_yy[i] * fe_0 + 2.0 * tr_y_yyy_y[i] * fe_0 + ts_yyy_yy[i] * fe_0 + tr_y_yyy_yy[i] * pa_y[i];

        tr_y_yyyy_yz[i] = 3.0 * tr_y_yy_yz[i] * fe_0 + tr_y_yyy_z[i] * fe_0 + ts_yyy_yz[i] * fe_0 + tr_y_yyy_yz[i] * pa_y[i];

        tr_y_yyyy_zz[i] = 3.0 * tr_y_yy_zz[i] * fe_0 + ts_yyy_zz[i] * fe_0 + tr_y_yyy_zz[i] * pa_y[i];
    }

    // Set up 156-162 components of targeted buffer : GD

    auto tr_y_yyyz_xx = pbuffer.data(idx_dip_gd + 156);

    auto tr_y_yyyz_xy = pbuffer.data(idx_dip_gd + 157);

    auto tr_y_yyyz_xz = pbuffer.data(idx_dip_gd + 158);

    auto tr_y_yyyz_yy = pbuffer.data(idx_dip_gd + 159);

    auto tr_y_yyyz_yz = pbuffer.data(idx_dip_gd + 160);

    auto tr_y_yyyz_zz = pbuffer.data(idx_dip_gd + 161);

#pragma omp simd aligned(pa_z,             \
                             tr_y_yyy_x,   \
                             tr_y_yyy_xx,  \
                             tr_y_yyy_xy,  \
                             tr_y_yyy_xz,  \
                             tr_y_yyy_y,   \
                             tr_y_yyy_yy,  \
                             tr_y_yyy_yz,  \
                             tr_y_yyy_z,   \
                             tr_y_yyy_zz,  \
                             tr_y_yyyz_xx, \
                             tr_y_yyyz_xy, \
                             tr_y_yyyz_xz, \
                             tr_y_yyyz_yy, \
                             tr_y_yyyz_yz, \
                             tr_y_yyyz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_xx[i] = tr_y_yyy_xx[i] * pa_z[i];

        tr_y_yyyz_xy[i] = tr_y_yyy_xy[i] * pa_z[i];

        tr_y_yyyz_xz[i] = tr_y_yyy_x[i] * fe_0 + tr_y_yyy_xz[i] * pa_z[i];

        tr_y_yyyz_yy[i] = tr_y_yyy_yy[i] * pa_z[i];

        tr_y_yyyz_yz[i] = tr_y_yyy_y[i] * fe_0 + tr_y_yyy_yz[i] * pa_z[i];

        tr_y_yyyz_zz[i] = 2.0 * tr_y_yyy_z[i] * fe_0 + tr_y_yyy_zz[i] * pa_z[i];
    }

    // Set up 162-168 components of targeted buffer : GD

    auto tr_y_yyzz_xx = pbuffer.data(idx_dip_gd + 162);

    auto tr_y_yyzz_xy = pbuffer.data(idx_dip_gd + 163);

    auto tr_y_yyzz_xz = pbuffer.data(idx_dip_gd + 164);

    auto tr_y_yyzz_yy = pbuffer.data(idx_dip_gd + 165);

    auto tr_y_yyzz_yz = pbuffer.data(idx_dip_gd + 166);

    auto tr_y_yyzz_zz = pbuffer.data(idx_dip_gd + 167);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_y_yy_xx,   \
                             tr_y_yy_xy,   \
                             tr_y_yy_yy,   \
                             tr_y_yy_yz,   \
                             tr_y_yyz_xx,  \
                             tr_y_yyz_xy,  \
                             tr_y_yyz_y,   \
                             tr_y_yyz_yy,  \
                             tr_y_yyz_yz,  \
                             tr_y_yyzz_xx, \
                             tr_y_yyzz_xy, \
                             tr_y_yyzz_xz, \
                             tr_y_yyzz_yy, \
                             tr_y_yyzz_yz, \
                             tr_y_yyzz_zz, \
                             tr_y_yzz_xz,  \
                             tr_y_yzz_zz,  \
                             tr_y_zz_xz,   \
                             tr_y_zz_zz,   \
                             ts_yzz_xz,    \
                             ts_yzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_xx[i] = tr_y_yy_xx[i] * fe_0 + tr_y_yyz_xx[i] * pa_z[i];

        tr_y_yyzz_xy[i] = tr_y_yy_xy[i] * fe_0 + tr_y_yyz_xy[i] * pa_z[i];

        tr_y_yyzz_xz[i] = tr_y_zz_xz[i] * fe_0 + ts_yzz_xz[i] * fe_0 + tr_y_yzz_xz[i] * pa_y[i];

        tr_y_yyzz_yy[i] = tr_y_yy_yy[i] * fe_0 + tr_y_yyz_yy[i] * pa_z[i];

        tr_y_yyzz_yz[i] = tr_y_yy_yz[i] * fe_0 + tr_y_yyz_y[i] * fe_0 + tr_y_yyz_yz[i] * pa_z[i];

        tr_y_yyzz_zz[i] = tr_y_zz_zz[i] * fe_0 + ts_yzz_zz[i] * fe_0 + tr_y_yzz_zz[i] * pa_y[i];
    }

    // Set up 168-174 components of targeted buffer : GD

    auto tr_y_yzzz_xx = pbuffer.data(idx_dip_gd + 168);

    auto tr_y_yzzz_xy = pbuffer.data(idx_dip_gd + 169);

    auto tr_y_yzzz_xz = pbuffer.data(idx_dip_gd + 170);

    auto tr_y_yzzz_yy = pbuffer.data(idx_dip_gd + 171);

    auto tr_y_yzzz_yz = pbuffer.data(idx_dip_gd + 172);

    auto tr_y_yzzz_zz = pbuffer.data(idx_dip_gd + 173);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_y_yz_xy,   \
                             tr_y_yz_yy,   \
                             tr_y_yzz_xy,  \
                             tr_y_yzz_yy,  \
                             tr_y_yzzz_xx, \
                             tr_y_yzzz_xy, \
                             tr_y_yzzz_xz, \
                             tr_y_yzzz_yy, \
                             tr_y_yzzz_yz, \
                             tr_y_yzzz_zz, \
                             tr_y_zzz_xx,  \
                             tr_y_zzz_xz,  \
                             tr_y_zzz_yz,  \
                             tr_y_zzz_z,   \
                             tr_y_zzz_zz,  \
                             ts_zzz_xx,    \
                             ts_zzz_xz,    \
                             ts_zzz_yz,    \
                             ts_zzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_xx[i] = ts_zzz_xx[i] * fe_0 + tr_y_zzz_xx[i] * pa_y[i];

        tr_y_yzzz_xy[i] = 2.0 * tr_y_yz_xy[i] * fe_0 + tr_y_yzz_xy[i] * pa_z[i];

        tr_y_yzzz_xz[i] = ts_zzz_xz[i] * fe_0 + tr_y_zzz_xz[i] * pa_y[i];

        tr_y_yzzz_yy[i] = 2.0 * tr_y_yz_yy[i] * fe_0 + tr_y_yzz_yy[i] * pa_z[i];

        tr_y_yzzz_yz[i] = tr_y_zzz_z[i] * fe_0 + ts_zzz_yz[i] * fe_0 + tr_y_zzz_yz[i] * pa_y[i];

        tr_y_yzzz_zz[i] = ts_zzz_zz[i] * fe_0 + tr_y_zzz_zz[i] * pa_y[i];
    }

    // Set up 174-180 components of targeted buffer : GD

    auto tr_y_zzzz_xx = pbuffer.data(idx_dip_gd + 174);

    auto tr_y_zzzz_xy = pbuffer.data(idx_dip_gd + 175);

    auto tr_y_zzzz_xz = pbuffer.data(idx_dip_gd + 176);

    auto tr_y_zzzz_yy = pbuffer.data(idx_dip_gd + 177);

    auto tr_y_zzzz_yz = pbuffer.data(idx_dip_gd + 178);

    auto tr_y_zzzz_zz = pbuffer.data(idx_dip_gd + 179);

#pragma omp simd aligned(pa_z,             \
                             tr_y_zz_xx,   \
                             tr_y_zz_xy,   \
                             tr_y_zz_xz,   \
                             tr_y_zz_yy,   \
                             tr_y_zz_yz,   \
                             tr_y_zz_zz,   \
                             tr_y_zzz_x,   \
                             tr_y_zzz_xx,  \
                             tr_y_zzz_xy,  \
                             tr_y_zzz_xz,  \
                             tr_y_zzz_y,   \
                             tr_y_zzz_yy,  \
                             tr_y_zzz_yz,  \
                             tr_y_zzz_z,   \
                             tr_y_zzz_zz,  \
                             tr_y_zzzz_xx, \
                             tr_y_zzzz_xy, \
                             tr_y_zzzz_xz, \
                             tr_y_zzzz_yy, \
                             tr_y_zzzz_yz, \
                             tr_y_zzzz_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_xx[i] = 3.0 * tr_y_zz_xx[i] * fe_0 + tr_y_zzz_xx[i] * pa_z[i];

        tr_y_zzzz_xy[i] = 3.0 * tr_y_zz_xy[i] * fe_0 + tr_y_zzz_xy[i] * pa_z[i];

        tr_y_zzzz_xz[i] = 3.0 * tr_y_zz_xz[i] * fe_0 + tr_y_zzz_x[i] * fe_0 + tr_y_zzz_xz[i] * pa_z[i];

        tr_y_zzzz_yy[i] = 3.0 * tr_y_zz_yy[i] * fe_0 + tr_y_zzz_yy[i] * pa_z[i];

        tr_y_zzzz_yz[i] = 3.0 * tr_y_zz_yz[i] * fe_0 + tr_y_zzz_y[i] * fe_0 + tr_y_zzz_yz[i] * pa_z[i];

        tr_y_zzzz_zz[i] = 3.0 * tr_y_zz_zz[i] * fe_0 + 2.0 * tr_y_zzz_z[i] * fe_0 + tr_y_zzz_zz[i] * pa_z[i];
    }

    // Set up 180-186 components of targeted buffer : GD

    auto tr_z_xxxx_xx = pbuffer.data(idx_dip_gd + 180);

    auto tr_z_xxxx_xy = pbuffer.data(idx_dip_gd + 181);

    auto tr_z_xxxx_xz = pbuffer.data(idx_dip_gd + 182);

    auto tr_z_xxxx_yy = pbuffer.data(idx_dip_gd + 183);

    auto tr_z_xxxx_yz = pbuffer.data(idx_dip_gd + 184);

    auto tr_z_xxxx_zz = pbuffer.data(idx_dip_gd + 185);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xx_xx,   \
                             tr_z_xx_xy,   \
                             tr_z_xx_xz,   \
                             tr_z_xx_yy,   \
                             tr_z_xx_yz,   \
                             tr_z_xx_zz,   \
                             tr_z_xxx_x,   \
                             tr_z_xxx_xx,  \
                             tr_z_xxx_xy,  \
                             tr_z_xxx_xz,  \
                             tr_z_xxx_y,   \
                             tr_z_xxx_yy,  \
                             tr_z_xxx_yz,  \
                             tr_z_xxx_z,   \
                             tr_z_xxx_zz,  \
                             tr_z_xxxx_xx, \
                             tr_z_xxxx_xy, \
                             tr_z_xxxx_xz, \
                             tr_z_xxxx_yy, \
                             tr_z_xxxx_yz, \
                             tr_z_xxxx_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_xx[i] = 3.0 * tr_z_xx_xx[i] * fe_0 + 2.0 * tr_z_xxx_x[i] * fe_0 + tr_z_xxx_xx[i] * pa_x[i];

        tr_z_xxxx_xy[i] = 3.0 * tr_z_xx_xy[i] * fe_0 + tr_z_xxx_y[i] * fe_0 + tr_z_xxx_xy[i] * pa_x[i];

        tr_z_xxxx_xz[i] = 3.0 * tr_z_xx_xz[i] * fe_0 + tr_z_xxx_z[i] * fe_0 + tr_z_xxx_xz[i] * pa_x[i];

        tr_z_xxxx_yy[i] = 3.0 * tr_z_xx_yy[i] * fe_0 + tr_z_xxx_yy[i] * pa_x[i];

        tr_z_xxxx_yz[i] = 3.0 * tr_z_xx_yz[i] * fe_0 + tr_z_xxx_yz[i] * pa_x[i];

        tr_z_xxxx_zz[i] = 3.0 * tr_z_xx_zz[i] * fe_0 + tr_z_xxx_zz[i] * pa_x[i];
    }

    // Set up 186-192 components of targeted buffer : GD

    auto tr_z_xxxy_xx = pbuffer.data(idx_dip_gd + 186);

    auto tr_z_xxxy_xy = pbuffer.data(idx_dip_gd + 187);

    auto tr_z_xxxy_xz = pbuffer.data(idx_dip_gd + 188);

    auto tr_z_xxxy_yy = pbuffer.data(idx_dip_gd + 189);

    auto tr_z_xxxy_yz = pbuffer.data(idx_dip_gd + 190);

    auto tr_z_xxxy_zz = pbuffer.data(idx_dip_gd + 191);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xxx_x,   \
                             tr_z_xxx_xx,  \
                             tr_z_xxx_xy,  \
                             tr_z_xxx_xz,  \
                             tr_z_xxx_zz,  \
                             tr_z_xxxy_xx, \
                             tr_z_xxxy_xy, \
                             tr_z_xxxy_xz, \
                             tr_z_xxxy_yy, \
                             tr_z_xxxy_yz, \
                             tr_z_xxxy_zz, \
                             tr_z_xxy_yy,  \
                             tr_z_xxy_yz,  \
                             tr_z_xy_yy,   \
                             tr_z_xy_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_xx[i] = tr_z_xxx_xx[i] * pa_y[i];

        tr_z_xxxy_xy[i] = tr_z_xxx_x[i] * fe_0 + tr_z_xxx_xy[i] * pa_y[i];

        tr_z_xxxy_xz[i] = tr_z_xxx_xz[i] * pa_y[i];

        tr_z_xxxy_yy[i] = 2.0 * tr_z_xy_yy[i] * fe_0 + tr_z_xxy_yy[i] * pa_x[i];

        tr_z_xxxy_yz[i] = 2.0 * tr_z_xy_yz[i] * fe_0 + tr_z_xxy_yz[i] * pa_x[i];

        tr_z_xxxy_zz[i] = tr_z_xxx_zz[i] * pa_y[i];
    }

    // Set up 192-198 components of targeted buffer : GD

    auto tr_z_xxxz_xx = pbuffer.data(idx_dip_gd + 192);

    auto tr_z_xxxz_xy = pbuffer.data(idx_dip_gd + 193);

    auto tr_z_xxxz_xz = pbuffer.data(idx_dip_gd + 194);

    auto tr_z_xxxz_yy = pbuffer.data(idx_dip_gd + 195);

    auto tr_z_xxxz_yz = pbuffer.data(idx_dip_gd + 196);

    auto tr_z_xxxz_zz = pbuffer.data(idx_dip_gd + 197);

#pragma omp simd aligned(pa_x,             \
                             pa_z,         \
                             tr_z_xxx_xx,  \
                             tr_z_xxx_xy,  \
                             tr_z_xxxz_xx, \
                             tr_z_xxxz_xy, \
                             tr_z_xxxz_xz, \
                             tr_z_xxxz_yy, \
                             tr_z_xxxz_yz, \
                             tr_z_xxxz_zz, \
                             tr_z_xxz_xz,  \
                             tr_z_xxz_yy,  \
                             tr_z_xxz_yz,  \
                             tr_z_xxz_z,   \
                             tr_z_xxz_zz,  \
                             tr_z_xz_xz,   \
                             tr_z_xz_yy,   \
                             tr_z_xz_yz,   \
                             tr_z_xz_zz,   \
                             ts_xxx_xx,    \
                             ts_xxx_xy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_xx[i] = ts_xxx_xx[i] * fe_0 + tr_z_xxx_xx[i] * pa_z[i];

        tr_z_xxxz_xy[i] = ts_xxx_xy[i] * fe_0 + tr_z_xxx_xy[i] * pa_z[i];

        tr_z_xxxz_xz[i] = 2.0 * tr_z_xz_xz[i] * fe_0 + tr_z_xxz_z[i] * fe_0 + tr_z_xxz_xz[i] * pa_x[i];

        tr_z_xxxz_yy[i] = 2.0 * tr_z_xz_yy[i] * fe_0 + tr_z_xxz_yy[i] * pa_x[i];

        tr_z_xxxz_yz[i] = 2.0 * tr_z_xz_yz[i] * fe_0 + tr_z_xxz_yz[i] * pa_x[i];

        tr_z_xxxz_zz[i] = 2.0 * tr_z_xz_zz[i] * fe_0 + tr_z_xxz_zz[i] * pa_x[i];
    }

    // Set up 198-204 components of targeted buffer : GD

    auto tr_z_xxyy_xx = pbuffer.data(idx_dip_gd + 198);

    auto tr_z_xxyy_xy = pbuffer.data(idx_dip_gd + 199);

    auto tr_z_xxyy_xz = pbuffer.data(idx_dip_gd + 200);

    auto tr_z_xxyy_yy = pbuffer.data(idx_dip_gd + 201);

    auto tr_z_xxyy_yz = pbuffer.data(idx_dip_gd + 202);

    auto tr_z_xxyy_zz = pbuffer.data(idx_dip_gd + 203);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xx_xx,   \
                             tr_z_xx_xz,   \
                             tr_z_xxy_xx,  \
                             tr_z_xxy_xz,  \
                             tr_z_xxyy_xx, \
                             tr_z_xxyy_xy, \
                             tr_z_xxyy_xz, \
                             tr_z_xxyy_yy, \
                             tr_z_xxyy_yz, \
                             tr_z_xxyy_zz, \
                             tr_z_xyy_xy,  \
                             tr_z_xyy_y,   \
                             tr_z_xyy_yy,  \
                             tr_z_xyy_yz,  \
                             tr_z_xyy_zz,  \
                             tr_z_yy_xy,   \
                             tr_z_yy_yy,   \
                             tr_z_yy_yz,   \
                             tr_z_yy_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_xx[i] = tr_z_xx_xx[i] * fe_0 + tr_z_xxy_xx[i] * pa_y[i];

        tr_z_xxyy_xy[i] = tr_z_yy_xy[i] * fe_0 + tr_z_xyy_y[i] * fe_0 + tr_z_xyy_xy[i] * pa_x[i];

        tr_z_xxyy_xz[i] = tr_z_xx_xz[i] * fe_0 + tr_z_xxy_xz[i] * pa_y[i];

        tr_z_xxyy_yy[i] = tr_z_yy_yy[i] * fe_0 + tr_z_xyy_yy[i] * pa_x[i];

        tr_z_xxyy_yz[i] = tr_z_yy_yz[i] * fe_0 + tr_z_xyy_yz[i] * pa_x[i];

        tr_z_xxyy_zz[i] = tr_z_yy_zz[i] * fe_0 + tr_z_xyy_zz[i] * pa_x[i];
    }

    // Set up 204-210 components of targeted buffer : GD

    auto tr_z_xxyz_xx = pbuffer.data(idx_dip_gd + 204);

    auto tr_z_xxyz_xy = pbuffer.data(idx_dip_gd + 205);

    auto tr_z_xxyz_xz = pbuffer.data(idx_dip_gd + 206);

    auto tr_z_xxyz_yy = pbuffer.data(idx_dip_gd + 207);

    auto tr_z_xxyz_yz = pbuffer.data(idx_dip_gd + 208);

    auto tr_z_xxyz_zz = pbuffer.data(idx_dip_gd + 209);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xxyz_xx, \
                             tr_z_xxyz_xy, \
                             tr_z_xxyz_xz, \
                             tr_z_xxyz_yy, \
                             tr_z_xxyz_yz, \
                             tr_z_xxyz_zz, \
                             tr_z_xxz_x,   \
                             tr_z_xxz_xx,  \
                             tr_z_xxz_xy,  \
                             tr_z_xxz_xz,  \
                             tr_z_xxz_zz,  \
                             tr_z_xyz_yy,  \
                             tr_z_xyz_yz,  \
                             tr_z_yz_yy,   \
                             tr_z_yz_yz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_xx[i] = tr_z_xxz_xx[i] * pa_y[i];

        tr_z_xxyz_xy[i] = tr_z_xxz_x[i] * fe_0 + tr_z_xxz_xy[i] * pa_y[i];

        tr_z_xxyz_xz[i] = tr_z_xxz_xz[i] * pa_y[i];

        tr_z_xxyz_yy[i] = tr_z_yz_yy[i] * fe_0 + tr_z_xyz_yy[i] * pa_x[i];

        tr_z_xxyz_yz[i] = tr_z_yz_yz[i] * fe_0 + tr_z_xyz_yz[i] * pa_x[i];

        tr_z_xxyz_zz[i] = tr_z_xxz_zz[i] * pa_y[i];
    }

    // Set up 210-216 components of targeted buffer : GD

    auto tr_z_xxzz_xx = pbuffer.data(idx_dip_gd + 210);

    auto tr_z_xxzz_xy = pbuffer.data(idx_dip_gd + 211);

    auto tr_z_xxzz_xz = pbuffer.data(idx_dip_gd + 212);

    auto tr_z_xxzz_yy = pbuffer.data(idx_dip_gd + 213);

    auto tr_z_xxzz_yz = pbuffer.data(idx_dip_gd + 214);

    auto tr_z_xxzz_zz = pbuffer.data(idx_dip_gd + 215);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xxzz_xx, \
                             tr_z_xxzz_xy, \
                             tr_z_xxzz_xz, \
                             tr_z_xxzz_yy, \
                             tr_z_xxzz_yz, \
                             tr_z_xxzz_zz, \
                             tr_z_xzz_x,   \
                             tr_z_xzz_xx,  \
                             tr_z_xzz_xy,  \
                             tr_z_xzz_xz,  \
                             tr_z_xzz_y,   \
                             tr_z_xzz_yy,  \
                             tr_z_xzz_yz,  \
                             tr_z_xzz_z,   \
                             tr_z_xzz_zz,  \
                             tr_z_zz_xx,   \
                             tr_z_zz_xy,   \
                             tr_z_zz_xz,   \
                             tr_z_zz_yy,   \
                             tr_z_zz_yz,   \
                             tr_z_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_xx[i] = tr_z_zz_xx[i] * fe_0 + 2.0 * tr_z_xzz_x[i] * fe_0 + tr_z_xzz_xx[i] * pa_x[i];

        tr_z_xxzz_xy[i] = tr_z_zz_xy[i] * fe_0 + tr_z_xzz_y[i] * fe_0 + tr_z_xzz_xy[i] * pa_x[i];

        tr_z_xxzz_xz[i] = tr_z_zz_xz[i] * fe_0 + tr_z_xzz_z[i] * fe_0 + tr_z_xzz_xz[i] * pa_x[i];

        tr_z_xxzz_yy[i] = tr_z_zz_yy[i] * fe_0 + tr_z_xzz_yy[i] * pa_x[i];

        tr_z_xxzz_yz[i] = tr_z_zz_yz[i] * fe_0 + tr_z_xzz_yz[i] * pa_x[i];

        tr_z_xxzz_zz[i] = tr_z_zz_zz[i] * fe_0 + tr_z_xzz_zz[i] * pa_x[i];
    }

    // Set up 216-222 components of targeted buffer : GD

    auto tr_z_xyyy_xx = pbuffer.data(idx_dip_gd + 216);

    auto tr_z_xyyy_xy = pbuffer.data(idx_dip_gd + 217);

    auto tr_z_xyyy_xz = pbuffer.data(idx_dip_gd + 218);

    auto tr_z_xyyy_yy = pbuffer.data(idx_dip_gd + 219);

    auto tr_z_xyyy_yz = pbuffer.data(idx_dip_gd + 220);

    auto tr_z_xyyy_zz = pbuffer.data(idx_dip_gd + 221);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xyyy_xx, \
                             tr_z_xyyy_xy, \
                             tr_z_xyyy_xz, \
                             tr_z_xyyy_yy, \
                             tr_z_xyyy_yz, \
                             tr_z_xyyy_zz, \
                             tr_z_yyy_x,   \
                             tr_z_yyy_xx,  \
                             tr_z_yyy_xy,  \
                             tr_z_yyy_xz,  \
                             tr_z_yyy_y,   \
                             tr_z_yyy_yy,  \
                             tr_z_yyy_yz,  \
                             tr_z_yyy_z,   \
                             tr_z_yyy_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_xx[i] = 2.0 * tr_z_yyy_x[i] * fe_0 + tr_z_yyy_xx[i] * pa_x[i];

        tr_z_xyyy_xy[i] = tr_z_yyy_y[i] * fe_0 + tr_z_yyy_xy[i] * pa_x[i];

        tr_z_xyyy_xz[i] = tr_z_yyy_z[i] * fe_0 + tr_z_yyy_xz[i] * pa_x[i];

        tr_z_xyyy_yy[i] = tr_z_yyy_yy[i] * pa_x[i];

        tr_z_xyyy_yz[i] = tr_z_yyy_yz[i] * pa_x[i];

        tr_z_xyyy_zz[i] = tr_z_yyy_zz[i] * pa_x[i];
    }

    // Set up 222-228 components of targeted buffer : GD

    auto tr_z_xyyz_xx = pbuffer.data(idx_dip_gd + 222);

    auto tr_z_xyyz_xy = pbuffer.data(idx_dip_gd + 223);

    auto tr_z_xyyz_xz = pbuffer.data(idx_dip_gd + 224);

    auto tr_z_xyyz_yy = pbuffer.data(idx_dip_gd + 225);

    auto tr_z_xyyz_yz = pbuffer.data(idx_dip_gd + 226);

    auto tr_z_xyyz_zz = pbuffer.data(idx_dip_gd + 227);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xyyz_xx, \
                             tr_z_xyyz_xy, \
                             tr_z_xyyz_xz, \
                             tr_z_xyyz_yy, \
                             tr_z_xyyz_yz, \
                             tr_z_xyyz_zz, \
                             tr_z_yyz_x,   \
                             tr_z_yyz_xx,  \
                             tr_z_yyz_xy,  \
                             tr_z_yyz_xz,  \
                             tr_z_yyz_y,   \
                             tr_z_yyz_yy,  \
                             tr_z_yyz_yz,  \
                             tr_z_yyz_z,   \
                             tr_z_yyz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_xx[i] = 2.0 * tr_z_yyz_x[i] * fe_0 + tr_z_yyz_xx[i] * pa_x[i];

        tr_z_xyyz_xy[i] = tr_z_yyz_y[i] * fe_0 + tr_z_yyz_xy[i] * pa_x[i];

        tr_z_xyyz_xz[i] = tr_z_yyz_z[i] * fe_0 + tr_z_yyz_xz[i] * pa_x[i];

        tr_z_xyyz_yy[i] = tr_z_yyz_yy[i] * pa_x[i];

        tr_z_xyyz_yz[i] = tr_z_yyz_yz[i] * pa_x[i];

        tr_z_xyyz_zz[i] = tr_z_yyz_zz[i] * pa_x[i];
    }

    // Set up 228-234 components of targeted buffer : GD

    auto tr_z_xyzz_xx = pbuffer.data(idx_dip_gd + 228);

    auto tr_z_xyzz_xy = pbuffer.data(idx_dip_gd + 229);

    auto tr_z_xyzz_xz = pbuffer.data(idx_dip_gd + 230);

    auto tr_z_xyzz_yy = pbuffer.data(idx_dip_gd + 231);

    auto tr_z_xyzz_yz = pbuffer.data(idx_dip_gd + 232);

    auto tr_z_xyzz_zz = pbuffer.data(idx_dip_gd + 233);

#pragma omp simd aligned(pa_x,             \
                             pa_y,         \
                             tr_z_xyzz_xx, \
                             tr_z_xyzz_xy, \
                             tr_z_xyzz_xz, \
                             tr_z_xyzz_yy, \
                             tr_z_xyzz_yz, \
                             tr_z_xyzz_zz, \
                             tr_z_xzz_xx,  \
                             tr_z_xzz_xz,  \
                             tr_z_yzz_xy,  \
                             tr_z_yzz_y,   \
                             tr_z_yzz_yy,  \
                             tr_z_yzz_yz,  \
                             tr_z_yzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzz_xx[i] = tr_z_xzz_xx[i] * pa_y[i];

        tr_z_xyzz_xy[i] = tr_z_yzz_y[i] * fe_0 + tr_z_yzz_xy[i] * pa_x[i];

        tr_z_xyzz_xz[i] = tr_z_xzz_xz[i] * pa_y[i];

        tr_z_xyzz_yy[i] = tr_z_yzz_yy[i] * pa_x[i];

        tr_z_xyzz_yz[i] = tr_z_yzz_yz[i] * pa_x[i];

        tr_z_xyzz_zz[i] = tr_z_yzz_zz[i] * pa_x[i];
    }

    // Set up 234-240 components of targeted buffer : GD

    auto tr_z_xzzz_xx = pbuffer.data(idx_dip_gd + 234);

    auto tr_z_xzzz_xy = pbuffer.data(idx_dip_gd + 235);

    auto tr_z_xzzz_xz = pbuffer.data(idx_dip_gd + 236);

    auto tr_z_xzzz_yy = pbuffer.data(idx_dip_gd + 237);

    auto tr_z_xzzz_yz = pbuffer.data(idx_dip_gd + 238);

    auto tr_z_xzzz_zz = pbuffer.data(idx_dip_gd + 239);

#pragma omp simd aligned(pa_x,             \
                             tr_z_xzzz_xx, \
                             tr_z_xzzz_xy, \
                             tr_z_xzzz_xz, \
                             tr_z_xzzz_yy, \
                             tr_z_xzzz_yz, \
                             tr_z_xzzz_zz, \
                             tr_z_zzz_x,   \
                             tr_z_zzz_xx,  \
                             tr_z_zzz_xy,  \
                             tr_z_zzz_xz,  \
                             tr_z_zzz_y,   \
                             tr_z_zzz_yy,  \
                             tr_z_zzz_yz,  \
                             tr_z_zzz_z,   \
                             tr_z_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_xx[i] = 2.0 * tr_z_zzz_x[i] * fe_0 + tr_z_zzz_xx[i] * pa_x[i];

        tr_z_xzzz_xy[i] = tr_z_zzz_y[i] * fe_0 + tr_z_zzz_xy[i] * pa_x[i];

        tr_z_xzzz_xz[i] = tr_z_zzz_z[i] * fe_0 + tr_z_zzz_xz[i] * pa_x[i];

        tr_z_xzzz_yy[i] = tr_z_zzz_yy[i] * pa_x[i];

        tr_z_xzzz_yz[i] = tr_z_zzz_yz[i] * pa_x[i];

        tr_z_xzzz_zz[i] = tr_z_zzz_zz[i] * pa_x[i];
    }

    // Set up 240-246 components of targeted buffer : GD

    auto tr_z_yyyy_xx = pbuffer.data(idx_dip_gd + 240);

    auto tr_z_yyyy_xy = pbuffer.data(idx_dip_gd + 241);

    auto tr_z_yyyy_xz = pbuffer.data(idx_dip_gd + 242);

    auto tr_z_yyyy_yy = pbuffer.data(idx_dip_gd + 243);

    auto tr_z_yyyy_yz = pbuffer.data(idx_dip_gd + 244);

    auto tr_z_yyyy_zz = pbuffer.data(idx_dip_gd + 245);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yy_xx,   \
                             tr_z_yy_xy,   \
                             tr_z_yy_xz,   \
                             tr_z_yy_yy,   \
                             tr_z_yy_yz,   \
                             tr_z_yy_zz,   \
                             tr_z_yyy_x,   \
                             tr_z_yyy_xx,  \
                             tr_z_yyy_xy,  \
                             tr_z_yyy_xz,  \
                             tr_z_yyy_y,   \
                             tr_z_yyy_yy,  \
                             tr_z_yyy_yz,  \
                             tr_z_yyy_z,   \
                             tr_z_yyy_zz,  \
                             tr_z_yyyy_xx, \
                             tr_z_yyyy_xy, \
                             tr_z_yyyy_xz, \
                             tr_z_yyyy_yy, \
                             tr_z_yyyy_yz, \
                             tr_z_yyyy_zz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_xx[i] = 3.0 * tr_z_yy_xx[i] * fe_0 + tr_z_yyy_xx[i] * pa_y[i];

        tr_z_yyyy_xy[i] = 3.0 * tr_z_yy_xy[i] * fe_0 + tr_z_yyy_x[i] * fe_0 + tr_z_yyy_xy[i] * pa_y[i];

        tr_z_yyyy_xz[i] = 3.0 * tr_z_yy_xz[i] * fe_0 + tr_z_yyy_xz[i] * pa_y[i];

        tr_z_yyyy_yy[i] = 3.0 * tr_z_yy_yy[i] * fe_0 + 2.0 * tr_z_yyy_y[i] * fe_0 + tr_z_yyy_yy[i] * pa_y[i];

        tr_z_yyyy_yz[i] = 3.0 * tr_z_yy_yz[i] * fe_0 + tr_z_yyy_z[i] * fe_0 + tr_z_yyy_yz[i] * pa_y[i];

        tr_z_yyyy_zz[i] = 3.0 * tr_z_yy_zz[i] * fe_0 + tr_z_yyy_zz[i] * pa_y[i];
    }

    // Set up 246-252 components of targeted buffer : GD

    auto tr_z_yyyz_xx = pbuffer.data(idx_dip_gd + 246);

    auto tr_z_yyyz_xy = pbuffer.data(idx_dip_gd + 247);

    auto tr_z_yyyz_xz = pbuffer.data(idx_dip_gd + 248);

    auto tr_z_yyyz_yy = pbuffer.data(idx_dip_gd + 249);

    auto tr_z_yyyz_yz = pbuffer.data(idx_dip_gd + 250);

    auto tr_z_yyyz_zz = pbuffer.data(idx_dip_gd + 251);

#pragma omp simd aligned(pa_y,             \
                             pa_z,         \
                             tr_z_yyy_xy,  \
                             tr_z_yyy_yy,  \
                             tr_z_yyyz_xx, \
                             tr_z_yyyz_xy, \
                             tr_z_yyyz_xz, \
                             tr_z_yyyz_yy, \
                             tr_z_yyyz_yz, \
                             tr_z_yyyz_zz, \
                             tr_z_yyz_xx,  \
                             tr_z_yyz_xz,  \
                             tr_z_yyz_yz,  \
                             tr_z_yyz_z,   \
                             tr_z_yyz_zz,  \
                             tr_z_yz_xx,   \
                             tr_z_yz_xz,   \
                             tr_z_yz_yz,   \
                             tr_z_yz_zz,   \
                             ts_yyy_xy,    \
                             ts_yyy_yy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_xx[i] = 2.0 * tr_z_yz_xx[i] * fe_0 + tr_z_yyz_xx[i] * pa_y[i];

        tr_z_yyyz_xy[i] = ts_yyy_xy[i] * fe_0 + tr_z_yyy_xy[i] * pa_z[i];

        tr_z_yyyz_xz[i] = 2.0 * tr_z_yz_xz[i] * fe_0 + tr_z_yyz_xz[i] * pa_y[i];

        tr_z_yyyz_yy[i] = ts_yyy_yy[i] * fe_0 + tr_z_yyy_yy[i] * pa_z[i];

        tr_z_yyyz_yz[i] = 2.0 * tr_z_yz_yz[i] * fe_0 + tr_z_yyz_z[i] * fe_0 + tr_z_yyz_yz[i] * pa_y[i];

        tr_z_yyyz_zz[i] = 2.0 * tr_z_yz_zz[i] * fe_0 + tr_z_yyz_zz[i] * pa_y[i];
    }

    // Set up 252-258 components of targeted buffer : GD

    auto tr_z_yyzz_xx = pbuffer.data(idx_dip_gd + 252);

    auto tr_z_yyzz_xy = pbuffer.data(idx_dip_gd + 253);

    auto tr_z_yyzz_xz = pbuffer.data(idx_dip_gd + 254);

    auto tr_z_yyzz_yy = pbuffer.data(idx_dip_gd + 255);

    auto tr_z_yyzz_yz = pbuffer.data(idx_dip_gd + 256);

    auto tr_z_yyzz_zz = pbuffer.data(idx_dip_gd + 257);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yyzz_xx, \
                             tr_z_yyzz_xy, \
                             tr_z_yyzz_xz, \
                             tr_z_yyzz_yy, \
                             tr_z_yyzz_yz, \
                             tr_z_yyzz_zz, \
                             tr_z_yzz_x,   \
                             tr_z_yzz_xx,  \
                             tr_z_yzz_xy,  \
                             tr_z_yzz_xz,  \
                             tr_z_yzz_y,   \
                             tr_z_yzz_yy,  \
                             tr_z_yzz_yz,  \
                             tr_z_yzz_z,   \
                             tr_z_yzz_zz,  \
                             tr_z_zz_xx,   \
                             tr_z_zz_xy,   \
                             tr_z_zz_xz,   \
                             tr_z_zz_yy,   \
                             tr_z_zz_yz,   \
                             tr_z_zz_zz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_xx[i] = tr_z_zz_xx[i] * fe_0 + tr_z_yzz_xx[i] * pa_y[i];

        tr_z_yyzz_xy[i] = tr_z_zz_xy[i] * fe_0 + tr_z_yzz_x[i] * fe_0 + tr_z_yzz_xy[i] * pa_y[i];

        tr_z_yyzz_xz[i] = tr_z_zz_xz[i] * fe_0 + tr_z_yzz_xz[i] * pa_y[i];

        tr_z_yyzz_yy[i] = tr_z_zz_yy[i] * fe_0 + 2.0 * tr_z_yzz_y[i] * fe_0 + tr_z_yzz_yy[i] * pa_y[i];

        tr_z_yyzz_yz[i] = tr_z_zz_yz[i] * fe_0 + tr_z_yzz_z[i] * fe_0 + tr_z_yzz_yz[i] * pa_y[i];

        tr_z_yyzz_zz[i] = tr_z_zz_zz[i] * fe_0 + tr_z_yzz_zz[i] * pa_y[i];
    }

    // Set up 258-264 components of targeted buffer : GD

    auto tr_z_yzzz_xx = pbuffer.data(idx_dip_gd + 258);

    auto tr_z_yzzz_xy = pbuffer.data(idx_dip_gd + 259);

    auto tr_z_yzzz_xz = pbuffer.data(idx_dip_gd + 260);

    auto tr_z_yzzz_yy = pbuffer.data(idx_dip_gd + 261);

    auto tr_z_yzzz_yz = pbuffer.data(idx_dip_gd + 262);

    auto tr_z_yzzz_zz = pbuffer.data(idx_dip_gd + 263);

#pragma omp simd aligned(pa_y,             \
                             tr_z_yzzz_xx, \
                             tr_z_yzzz_xy, \
                             tr_z_yzzz_xz, \
                             tr_z_yzzz_yy, \
                             tr_z_yzzz_yz, \
                             tr_z_yzzz_zz, \
                             tr_z_zzz_x,   \
                             tr_z_zzz_xx,  \
                             tr_z_zzz_xy,  \
                             tr_z_zzz_xz,  \
                             tr_z_zzz_y,   \
                             tr_z_zzz_yy,  \
                             tr_z_zzz_yz,  \
                             tr_z_zzz_z,   \
                             tr_z_zzz_zz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_xx[i] = tr_z_zzz_xx[i] * pa_y[i];

        tr_z_yzzz_xy[i] = tr_z_zzz_x[i] * fe_0 + tr_z_zzz_xy[i] * pa_y[i];

        tr_z_yzzz_xz[i] = tr_z_zzz_xz[i] * pa_y[i];

        tr_z_yzzz_yy[i] = 2.0 * tr_z_zzz_y[i] * fe_0 + tr_z_zzz_yy[i] * pa_y[i];

        tr_z_yzzz_yz[i] = tr_z_zzz_z[i] * fe_0 + tr_z_zzz_yz[i] * pa_y[i];

        tr_z_yzzz_zz[i] = tr_z_zzz_zz[i] * pa_y[i];
    }

    // Set up 264-270 components of targeted buffer : GD

    auto tr_z_zzzz_xx = pbuffer.data(idx_dip_gd + 264);

    auto tr_z_zzzz_xy = pbuffer.data(idx_dip_gd + 265);

    auto tr_z_zzzz_xz = pbuffer.data(idx_dip_gd + 266);

    auto tr_z_zzzz_yy = pbuffer.data(idx_dip_gd + 267);

    auto tr_z_zzzz_yz = pbuffer.data(idx_dip_gd + 268);

    auto tr_z_zzzz_zz = pbuffer.data(idx_dip_gd + 269);

#pragma omp simd aligned(pa_z,             \
                             tr_z_zz_xx,   \
                             tr_z_zz_xy,   \
                             tr_z_zz_xz,   \
                             tr_z_zz_yy,   \
                             tr_z_zz_yz,   \
                             tr_z_zz_zz,   \
                             tr_z_zzz_x,   \
                             tr_z_zzz_xx,  \
                             tr_z_zzz_xy,  \
                             tr_z_zzz_xz,  \
                             tr_z_zzz_y,   \
                             tr_z_zzz_yy,  \
                             tr_z_zzz_yz,  \
                             tr_z_zzz_z,   \
                             tr_z_zzz_zz,  \
                             tr_z_zzzz_xx, \
                             tr_z_zzzz_xy, \
                             tr_z_zzzz_xz, \
                             tr_z_zzzz_yy, \
                             tr_z_zzzz_yz, \
                             tr_z_zzzz_zz, \
                             ts_zzz_xx,    \
                             ts_zzz_xy,    \
                             ts_zzz_xz,    \
                             ts_zzz_yy,    \
                             ts_zzz_yz,    \
                             ts_zzz_zz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_xx[i] = 3.0 * tr_z_zz_xx[i] * fe_0 + ts_zzz_xx[i] * fe_0 + tr_z_zzz_xx[i] * pa_z[i];

        tr_z_zzzz_xy[i] = 3.0 * tr_z_zz_xy[i] * fe_0 + ts_zzz_xy[i] * fe_0 + tr_z_zzz_xy[i] * pa_z[i];

        tr_z_zzzz_xz[i] = 3.0 * tr_z_zz_xz[i] * fe_0 + tr_z_zzz_x[i] * fe_0 + ts_zzz_xz[i] * fe_0 + tr_z_zzz_xz[i] * pa_z[i];

        tr_z_zzzz_yy[i] = 3.0 * tr_z_zz_yy[i] * fe_0 + ts_zzz_yy[i] * fe_0 + tr_z_zzz_yy[i] * pa_z[i];

        tr_z_zzzz_yz[i] = 3.0 * tr_z_zz_yz[i] * fe_0 + tr_z_zzz_y[i] * fe_0 + ts_zzz_yz[i] * fe_0 + tr_z_zzz_yz[i] * pa_z[i];

        tr_z_zzzz_zz[i] = 3.0 * tr_z_zz_zz[i] * fe_0 + 2.0 * tr_z_zzz_z[i] * fe_0 + ts_zzz_zz[i] * fe_0 + tr_z_zzz_zz[i] * pa_z[i];
    }
}

}  // namespace diprec
