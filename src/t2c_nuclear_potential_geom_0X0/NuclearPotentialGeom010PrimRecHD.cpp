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

#include "NuclearPotentialGeom010PrimRecHD.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_hd(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_hd,
                                        const size_t              idx_npot_geom_010_0_fd,
                                        const size_t              idx_npot_geom_010_1_fd,
                                        const size_t              idx_npot_geom_010_0_gp,
                                        const size_t              idx_npot_geom_010_1_gp,
                                        const size_t              idx_npot_1_gd,
                                        const size_t              idx_npot_geom_010_0_gd,
                                        const size_t              idx_npot_geom_010_1_gd,
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

    // Set up components of auxiliary buffer : FD

    auto ta1_x_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd);

    auto ta1_x_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 1);

    auto ta1_x_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 2);

    auto ta1_x_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 3);

    auto ta1_x_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 4);

    auto ta1_x_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 5);

    auto ta1_x_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 6);

    auto ta1_x_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 7);

    auto ta1_x_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 8);

    auto ta1_x_xxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 11);

    auto ta1_x_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 12);

    auto ta1_x_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 13);

    auto ta1_x_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 14);

    auto ta1_x_xxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 15);

    auto ta1_x_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 17);

    auto ta1_x_xyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 18);

    auto ta1_x_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 19);

    auto ta1_x_xyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 20);

    auto ta1_x_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 21);

    auto ta1_x_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 22);

    auto ta1_x_xyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 26);

    auto ta1_x_xzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 30);

    auto ta1_x_xzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 31);

    auto ta1_x_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 32);

    auto ta1_x_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 34);

    auto ta1_x_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 35);

    auto ta1_x_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 36);

    auto ta1_x_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 37);

    auto ta1_x_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 38);

    auto ta1_x_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 39);

    auto ta1_x_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 40);

    auto ta1_x_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 41);

    auto ta1_x_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 43);

    auto ta1_x_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 44);

    auto ta1_x_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 45);

    auto ta1_x_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 47);

    auto ta1_x_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 48);

    auto ta1_x_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 50);

    auto ta1_x_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 52);

    auto ta1_x_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 53);

    auto ta1_x_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 54);

    auto ta1_x_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 55);

    auto ta1_x_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 56);

    auto ta1_x_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 57);

    auto ta1_x_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 58);

    auto ta1_x_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 59);

    auto ta1_y_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 60);

    auto ta1_y_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 61);

    auto ta1_y_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 62);

    auto ta1_y_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 63);

    auto ta1_y_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 64);

    auto ta1_y_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 65);

    auto ta1_y_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 66);

    auto ta1_y_xxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 67);

    auto ta1_y_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 68);

    auto ta1_y_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 69);

    auto ta1_y_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 70);

    auto ta1_y_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 72);

    auto ta1_y_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 73);

    auto ta1_y_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 76);

    auto ta1_y_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 77);

    auto ta1_y_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 79);

    auto ta1_y_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 81);

    auto ta1_y_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 82);

    auto ta1_y_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 83);

    auto ta1_y_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 88);

    auto ta1_y_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 92);

    auto ta1_y_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 93);

    auto ta1_y_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 94);

    auto ta1_y_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 95);

    auto ta1_y_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 96);

    auto ta1_y_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 97);

    auto ta1_y_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 98);

    auto ta1_y_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 99);

    auto ta1_y_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 100);

    auto ta1_y_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 101);

    auto ta1_y_yyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 102);

    auto ta1_y_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 103);

    auto ta1_y_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 105);

    auto ta1_y_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 106);

    auto ta1_y_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 107);

    auto ta1_y_yzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 109);

    auto ta1_y_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 110);

    auto ta1_y_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 111);

    auto ta1_y_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 112);

    auto ta1_y_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 113);

    auto ta1_y_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 114);

    auto ta1_y_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 115);

    auto ta1_y_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 116);

    auto ta1_y_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 117);

    auto ta1_y_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 118);

    auto ta1_y_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 119);

    auto ta1_z_xxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 120);

    auto ta1_z_xxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 121);

    auto ta1_z_xxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 122);

    auto ta1_z_xxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 123);

    auto ta1_z_xxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 124);

    auto ta1_z_xxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 125);

    auto ta1_z_xxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 126);

    auto ta1_z_xxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 128);

    auto ta1_z_xxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 129);

    auto ta1_z_xxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 130);

    auto ta1_z_xxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 132);

    auto ta1_z_xxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 133);

    auto ta1_z_xxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 134);

    auto ta1_z_xxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 136);

    auto ta1_z_xxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 137);

    auto ta1_z_xyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 139);

    auto ta1_z_xyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 141);

    auto ta1_z_xyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 142);

    auto ta1_z_xyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 143);

    auto ta1_z_xyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 148);

    auto ta1_z_xzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 152);

    auto ta1_z_xzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 153);

    auto ta1_z_xzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 154);

    auto ta1_z_xzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 155);

    auto ta1_z_yyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 156);

    auto ta1_z_yyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 157);

    auto ta1_z_yyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 158);

    auto ta1_z_yyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 159);

    auto ta1_z_yyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 160);

    auto ta1_z_yyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 161);

    auto ta1_z_yyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 163);

    auto ta1_z_yyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 164);

    auto ta1_z_yyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 165);

    auto ta1_z_yyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 166);

    auto ta1_z_yyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 167);

    auto ta1_z_yzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 168);

    auto ta1_z_yzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 170);

    auto ta1_z_yzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 171);

    auto ta1_z_yzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 172);

    auto ta1_z_yzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 173);

    auto ta1_z_zzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_fd + 174);

    auto ta1_z_zzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 175);

    auto ta1_z_zzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 176);

    auto ta1_z_zzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_fd + 177);

    auto ta1_z_zzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 178);

    auto ta1_z_zzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_fd + 179);

    // Set up components of auxiliary buffer : FD

    auto ta1_x_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd);

    auto ta1_x_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 1);

    auto ta1_x_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 2);

    auto ta1_x_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 3);

    auto ta1_x_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 4);

    auto ta1_x_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 5);

    auto ta1_x_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 6);

    auto ta1_x_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 7);

    auto ta1_x_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 8);

    auto ta1_x_xxy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 11);

    auto ta1_x_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 12);

    auto ta1_x_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 13);

    auto ta1_x_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 14);

    auto ta1_x_xxz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 15);

    auto ta1_x_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 17);

    auto ta1_x_xyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 18);

    auto ta1_x_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 19);

    auto ta1_x_xyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 20);

    auto ta1_x_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 21);

    auto ta1_x_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 22);

    auto ta1_x_xyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 26);

    auto ta1_x_xzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 30);

    auto ta1_x_xzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 31);

    auto ta1_x_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 32);

    auto ta1_x_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 34);

    auto ta1_x_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 35);

    auto ta1_x_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 36);

    auto ta1_x_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 37);

    auto ta1_x_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 38);

    auto ta1_x_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 39);

    auto ta1_x_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 40);

    auto ta1_x_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 41);

    auto ta1_x_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 43);

    auto ta1_x_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 44);

    auto ta1_x_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 45);

    auto ta1_x_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 47);

    auto ta1_x_yzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 48);

    auto ta1_x_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 50);

    auto ta1_x_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 52);

    auto ta1_x_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 53);

    auto ta1_x_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 54);

    auto ta1_x_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 55);

    auto ta1_x_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 56);

    auto ta1_x_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 57);

    auto ta1_x_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 58);

    auto ta1_x_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 59);

    auto ta1_y_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 60);

    auto ta1_y_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 61);

    auto ta1_y_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 62);

    auto ta1_y_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 63);

    auto ta1_y_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 64);

    auto ta1_y_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 65);

    auto ta1_y_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 66);

    auto ta1_y_xxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 67);

    auto ta1_y_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 68);

    auto ta1_y_xxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 69);

    auto ta1_y_xxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 70);

    auto ta1_y_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 72);

    auto ta1_y_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 73);

    auto ta1_y_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 76);

    auto ta1_y_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 77);

    auto ta1_y_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 79);

    auto ta1_y_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 81);

    auto ta1_y_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 82);

    auto ta1_y_xyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 83);

    auto ta1_y_xyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 88);

    auto ta1_y_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 92);

    auto ta1_y_xzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 93);

    auto ta1_y_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 94);

    auto ta1_y_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 95);

    auto ta1_y_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 96);

    auto ta1_y_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 97);

    auto ta1_y_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 98);

    auto ta1_y_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 99);

    auto ta1_y_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 100);

    auto ta1_y_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 101);

    auto ta1_y_yyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 102);

    auto ta1_y_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 103);

    auto ta1_y_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 105);

    auto ta1_y_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 106);

    auto ta1_y_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 107);

    auto ta1_y_yzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 109);

    auto ta1_y_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 110);

    auto ta1_y_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 111);

    auto ta1_y_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 112);

    auto ta1_y_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 113);

    auto ta1_y_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 114);

    auto ta1_y_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 115);

    auto ta1_y_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 116);

    auto ta1_y_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 117);

    auto ta1_y_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 118);

    auto ta1_y_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 119);

    auto ta1_z_xxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 120);

    auto ta1_z_xxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 121);

    auto ta1_z_xxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 122);

    auto ta1_z_xxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 123);

    auto ta1_z_xxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 124);

    auto ta1_z_xxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 125);

    auto ta1_z_xxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 126);

    auto ta1_z_xxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 128);

    auto ta1_z_xxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 129);

    auto ta1_z_xxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 130);

    auto ta1_z_xxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 132);

    auto ta1_z_xxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 133);

    auto ta1_z_xxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 134);

    auto ta1_z_xxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 136);

    auto ta1_z_xxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 137);

    auto ta1_z_xyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 139);

    auto ta1_z_xyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 141);

    auto ta1_z_xyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 142);

    auto ta1_z_xyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 143);

    auto ta1_z_xyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 148);

    auto ta1_z_xzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 152);

    auto ta1_z_xzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 153);

    auto ta1_z_xzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 154);

    auto ta1_z_xzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 155);

    auto ta1_z_yyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 156);

    auto ta1_z_yyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 157);

    auto ta1_z_yyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 158);

    auto ta1_z_yyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 159);

    auto ta1_z_yyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 160);

    auto ta1_z_yyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 161);

    auto ta1_z_yyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 163);

    auto ta1_z_yyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 164);

    auto ta1_z_yyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 165);

    auto ta1_z_yyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 166);

    auto ta1_z_yyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 167);

    auto ta1_z_yzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 168);

    auto ta1_z_yzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 170);

    auto ta1_z_yzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 171);

    auto ta1_z_yzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 172);

    auto ta1_z_yzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 173);

    auto ta1_z_zzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_fd + 174);

    auto ta1_z_zzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 175);

    auto ta1_z_zzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 176);

    auto ta1_z_zzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_fd + 177);

    auto ta1_z_zzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 178);

    auto ta1_z_zzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_fd + 179);

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp);

    auto ta1_x_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 1);

    auto ta1_x_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 2);

    auto ta1_x_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 3);

    auto ta1_x_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 6);

    auto ta1_x_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 8);

    auto ta1_x_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 9);

    auto ta1_x_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 10);

    auto ta1_x_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 15);

    auto ta1_x_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 16);

    auto ta1_x_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 17);

    auto ta1_x_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 27);

    auto ta1_x_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 30);

    auto ta1_x_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 31);

    auto ta1_x_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 32);

    auto ta1_x_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 38);

    auto ta1_x_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 41);

    auto ta1_x_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 42);

    auto ta1_x_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 43);

    auto ta1_x_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 44);

    auto ta1_y_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 45);

    auto ta1_y_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 46);

    auto ta1_y_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 47);

    auto ta1_y_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 54);

    auto ta1_y_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 55);

    auto ta1_y_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 62);

    auto ta1_y_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 64);

    auto ta1_y_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 74);

    auto ta1_y_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 75);

    auto ta1_y_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 76);

    auto ta1_y_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 77);

    auto ta1_y_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 79);

    auto ta1_y_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 80);

    auto ta1_y_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 81);

    auto ta1_y_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 82);

    auto ta1_y_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 83);

    auto ta1_y_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 85);

    auto ta1_y_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 87);

    auto ta1_y_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 88);

    auto ta1_y_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 89);

    auto ta1_z_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 90);

    auto ta1_z_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 91);

    auto ta1_z_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 92);

    auto ta1_z_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 100);

    auto ta1_z_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 105);

    auto ta1_z_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 107);

    auto ta1_z_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 109);

    auto ta1_z_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 119);

    auto ta1_z_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 120);

    auto ta1_z_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 121);

    auto ta1_z_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 122);

    auto ta1_z_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 125);

    auto ta1_z_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 126);

    auto ta1_z_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 127);

    auto ta1_z_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 128);

    auto ta1_z_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 130);

    auto ta1_z_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 131);

    auto ta1_z_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 132);

    auto ta1_z_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 133);

    auto ta1_z_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 134);

    // Set up components of auxiliary buffer : GP

    auto ta1_x_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp);

    auto ta1_x_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 1);

    auto ta1_x_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 2);

    auto ta1_x_xxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 3);

    auto ta1_x_xxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 6);

    auto ta1_x_xxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 8);

    auto ta1_x_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 9);

    auto ta1_x_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 10);

    auto ta1_x_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 15);

    auto ta1_x_xxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 16);

    auto ta1_x_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 17);

    auto ta1_x_xzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 27);

    auto ta1_x_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 30);

    auto ta1_x_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 31);

    auto ta1_x_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 32);

    auto ta1_x_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 38);

    auto ta1_x_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 41);

    auto ta1_x_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 42);

    auto ta1_x_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 43);

    auto ta1_x_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 44);

    auto ta1_y_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 45);

    auto ta1_y_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 46);

    auto ta1_y_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 47);

    auto ta1_y_xxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 54);

    auto ta1_y_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 55);

    auto ta1_y_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 62);

    auto ta1_y_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 64);

    auto ta1_y_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 74);

    auto ta1_y_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 75);

    auto ta1_y_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 76);

    auto ta1_y_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 77);

    auto ta1_y_yyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 79);

    auto ta1_y_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 80);

    auto ta1_y_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 81);

    auto ta1_y_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 82);

    auto ta1_y_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 83);

    auto ta1_y_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 85);

    auto ta1_y_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 87);

    auto ta1_y_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 88);

    auto ta1_y_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 89);

    auto ta1_z_xxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 90);

    auto ta1_z_xxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 91);

    auto ta1_z_xxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 92);

    auto ta1_z_xxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 100);

    auto ta1_z_xxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 105);

    auto ta1_z_xxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 107);

    auto ta1_z_xyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 109);

    auto ta1_z_xzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 119);

    auto ta1_z_yyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 120);

    auto ta1_z_yyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 121);

    auto ta1_z_yyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 122);

    auto ta1_z_yyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 125);

    auto ta1_z_yyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 126);

    auto ta1_z_yyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 127);

    auto ta1_z_yyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 128);

    auto ta1_z_yzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 130);

    auto ta1_z_yzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 131);

    auto ta1_z_zzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_gp + 132);

    auto ta1_z_zzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_gp + 133);

    auto ta1_z_zzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_gp + 134);

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_1 = pbuffer.data(idx_npot_1_gd);

    auto ta_xxxx_xy_1 = pbuffer.data(idx_npot_1_gd + 1);

    auto ta_xxxx_xz_1 = pbuffer.data(idx_npot_1_gd + 2);

    auto ta_xxxx_yy_1 = pbuffer.data(idx_npot_1_gd + 3);

    auto ta_xxxx_yz_1 = pbuffer.data(idx_npot_1_gd + 4);

    auto ta_xxxx_zz_1 = pbuffer.data(idx_npot_1_gd + 5);

    auto ta_xxxy_xx_1 = pbuffer.data(idx_npot_1_gd + 6);

    auto ta_xxxy_xy_1 = pbuffer.data(idx_npot_1_gd + 7);

    auto ta_xxxy_xz_1 = pbuffer.data(idx_npot_1_gd + 8);

    auto ta_xxxy_yy_1 = pbuffer.data(idx_npot_1_gd + 9);

    auto ta_xxxz_xx_1 = pbuffer.data(idx_npot_1_gd + 12);

    auto ta_xxxz_xy_1 = pbuffer.data(idx_npot_1_gd + 13);

    auto ta_xxxz_xz_1 = pbuffer.data(idx_npot_1_gd + 14);

    auto ta_xxxz_zz_1 = pbuffer.data(idx_npot_1_gd + 17);

    auto ta_xxyy_xx_1 = pbuffer.data(idx_npot_1_gd + 18);

    auto ta_xxyy_xy_1 = pbuffer.data(idx_npot_1_gd + 19);

    auto ta_xxyy_xz_1 = pbuffer.data(idx_npot_1_gd + 20);

    auto ta_xxyy_yy_1 = pbuffer.data(idx_npot_1_gd + 21);

    auto ta_xxyy_yz_1 = pbuffer.data(idx_npot_1_gd + 22);

    auto ta_xxzz_xx_1 = pbuffer.data(idx_npot_1_gd + 30);

    auto ta_xxzz_xy_1 = pbuffer.data(idx_npot_1_gd + 31);

    auto ta_xxzz_xz_1 = pbuffer.data(idx_npot_1_gd + 32);

    auto ta_xxzz_yz_1 = pbuffer.data(idx_npot_1_gd + 34);

    auto ta_xxzz_zz_1 = pbuffer.data(idx_npot_1_gd + 35);

    auto ta_xyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 36);

    auto ta_xyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 37);

    auto ta_xyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 39);

    auto ta_xyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 40);

    auto ta_xzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 54);

    auto ta_xzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 56);

    auto ta_xzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 58);

    auto ta_xzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 59);

    auto ta_yyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 60);

    auto ta_yyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 61);

    auto ta_yyyy_xz_1 = pbuffer.data(idx_npot_1_gd + 62);

    auto ta_yyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 63);

    auto ta_yyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 64);

    auto ta_yyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 65);

    auto ta_yyyz_xy_1 = pbuffer.data(idx_npot_1_gd + 67);

    auto ta_yyyz_yy_1 = pbuffer.data(idx_npot_1_gd + 69);

    auto ta_yyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 70);

    auto ta_yyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 71);

    auto ta_yyzz_xy_1 = pbuffer.data(idx_npot_1_gd + 73);

    auto ta_yyzz_xz_1 = pbuffer.data(idx_npot_1_gd + 74);

    auto ta_yyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 75);

    auto ta_yyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 76);

    auto ta_yyzz_zz_1 = pbuffer.data(idx_npot_1_gd + 77);

    auto ta_yzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 80);

    auto ta_yzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 81);

    auto ta_yzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 82);

    auto ta_yzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 83);

    auto ta_zzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 84);

    auto ta_zzzz_xy_1 = pbuffer.data(idx_npot_1_gd + 85);

    auto ta_zzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 86);

    auto ta_zzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 87);

    auto ta_zzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 88);

    auto ta_zzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto ta1_x_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd);

    auto ta1_x_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 1);

    auto ta1_x_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 2);

    auto ta1_x_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 3);

    auto ta1_x_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 4);

    auto ta1_x_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 5);

    auto ta1_x_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 6);

    auto ta1_x_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 7);

    auto ta1_x_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 8);

    auto ta1_x_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 9);

    auto ta1_x_xxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 11);

    auto ta1_x_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 12);

    auto ta1_x_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 13);

    auto ta1_x_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 14);

    auto ta1_x_xxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 15);

    auto ta1_x_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 16);

    auto ta1_x_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 17);

    auto ta1_x_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 18);

    auto ta1_x_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 19);

    auto ta1_x_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 20);

    auto ta1_x_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 21);

    auto ta1_x_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 22);

    auto ta1_x_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 23);

    auto ta1_x_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 26);

    auto ta1_x_xxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 29);

    auto ta1_x_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 30);

    auto ta1_x_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 31);

    auto ta1_x_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 32);

    auto ta1_x_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 33);

    auto ta1_x_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 34);

    auto ta1_x_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 35);

    auto ta1_x_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 36);

    auto ta1_x_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 37);

    auto ta1_x_xyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 38);

    auto ta1_x_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 39);

    auto ta1_x_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 40);

    auto ta1_x_xyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 43);

    auto ta1_x_xyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 44);

    auto ta1_x_xyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 48);

    auto ta1_x_xyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 50);

    auto ta1_x_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 54);

    auto ta1_x_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 55);

    auto ta1_x_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 56);

    auto ta1_x_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 58);

    auto ta1_x_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 59);

    auto ta1_x_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 60);

    auto ta1_x_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 61);

    auto ta1_x_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 62);

    auto ta1_x_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 63);

    auto ta1_x_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 64);

    auto ta1_x_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 65);

    auto ta1_x_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 67);

    auto ta1_x_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 68);

    auto ta1_x_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 69);

    auto ta1_x_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 70);

    auto ta1_x_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 71);

    auto ta1_x_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 72);

    auto ta1_x_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 73);

    auto ta1_x_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 74);

    auto ta1_x_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 75);

    auto ta1_x_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 76);

    auto ta1_x_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 77);

    auto ta1_x_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 78);

    auto ta1_x_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 80);

    auto ta1_x_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 81);

    auto ta1_x_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 82);

    auto ta1_x_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 83);

    auto ta1_x_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 84);

    auto ta1_x_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 85);

    auto ta1_x_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 86);

    auto ta1_x_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 87);

    auto ta1_x_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 88);

    auto ta1_x_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 89);

    auto ta1_y_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 90);

    auto ta1_y_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 91);

    auto ta1_y_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 92);

    auto ta1_y_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 93);

    auto ta1_y_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 94);

    auto ta1_y_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 95);

    auto ta1_y_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 96);

    auto ta1_y_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 97);

    auto ta1_y_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 98);

    auto ta1_y_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 99);

    auto ta1_y_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 100);

    auto ta1_y_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 102);

    auto ta1_y_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 103);

    auto ta1_y_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 104);

    auto ta1_y_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 106);

    auto ta1_y_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 107);

    auto ta1_y_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 108);

    auto ta1_y_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 109);

    auto ta1_y_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 110);

    auto ta1_y_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 111);

    auto ta1_y_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 112);

    auto ta1_y_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 113);

    auto ta1_y_xxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 115);

    auto ta1_y_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 118);

    auto ta1_y_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 120);

    auto ta1_y_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 121);

    auto ta1_y_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 122);

    auto ta1_y_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 123);

    auto ta1_y_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 124);

    auto ta1_y_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 125);

    auto ta1_y_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 126);

    auto ta1_y_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 127);

    auto ta1_y_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 129);

    auto ta1_y_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 130);

    auto ta1_y_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 131);

    auto ta1_y_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 136);

    auto ta1_y_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 137);

    auto ta1_y_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 141);

    auto ta1_y_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 142);

    auto ta1_y_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 144);

    auto ta1_y_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 146);

    auto ta1_y_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 147);

    auto ta1_y_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 148);

    auto ta1_y_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 149);

    auto ta1_y_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 150);

    auto ta1_y_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 151);

    auto ta1_y_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 152);

    auto ta1_y_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 153);

    auto ta1_y_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 154);

    auto ta1_y_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 155);

    auto ta1_y_yyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 156);

    auto ta1_y_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 157);

    auto ta1_y_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 158);

    auto ta1_y_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 159);

    auto ta1_y_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 160);

    auto ta1_y_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 161);

    auto ta1_y_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 162);

    auto ta1_y_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 163);

    auto ta1_y_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 164);

    auto ta1_y_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 165);

    auto ta1_y_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 166);

    auto ta1_y_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 167);

    auto ta1_y_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 169);

    auto ta1_y_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 170);

    auto ta1_y_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 171);

    auto ta1_y_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 172);

    auto ta1_y_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 173);

    auto ta1_y_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 174);

    auto ta1_y_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 175);

    auto ta1_y_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 176);

    auto ta1_y_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 177);

    auto ta1_y_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 178);

    auto ta1_y_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 179);

    auto ta1_z_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 180);

    auto ta1_z_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 181);

    auto ta1_z_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 182);

    auto ta1_z_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 183);

    auto ta1_z_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 184);

    auto ta1_z_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 185);

    auto ta1_z_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 186);

    auto ta1_z_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 187);

    auto ta1_z_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 188);

    auto ta1_z_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 189);

    auto ta1_z_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 190);

    auto ta1_z_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 192);

    auto ta1_z_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 193);

    auto ta1_z_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 194);

    auto ta1_z_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 196);

    auto ta1_z_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 197);

    auto ta1_z_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 198);

    auto ta1_z_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 199);

    auto ta1_z_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 200);

    auto ta1_z_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 201);

    auto ta1_z_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 202);

    auto ta1_z_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 203);

    auto ta1_z_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 206);

    auto ta1_z_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 208);

    auto ta1_z_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 210);

    auto ta1_z_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 211);

    auto ta1_z_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 212);

    auto ta1_z_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 213);

    auto ta1_z_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 214);

    auto ta1_z_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 215);

    auto ta1_z_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 216);

    auto ta1_z_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 217);

    auto ta1_z_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 219);

    auto ta1_z_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 220);

    auto ta1_z_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 221);

    auto ta1_z_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 226);

    auto ta1_z_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 227);

    auto ta1_z_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 231);

    auto ta1_z_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 232);

    auto ta1_z_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 234);

    auto ta1_z_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 236);

    auto ta1_z_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 237);

    auto ta1_z_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 238);

    auto ta1_z_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 239);

    auto ta1_z_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 240);

    auto ta1_z_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 241);

    auto ta1_z_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 242);

    auto ta1_z_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 243);

    auto ta1_z_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 244);

    auto ta1_z_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 245);

    auto ta1_z_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 247);

    auto ta1_z_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 248);

    auto ta1_z_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 249);

    auto ta1_z_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 250);

    auto ta1_z_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 251);

    auto ta1_z_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 252);

    auto ta1_z_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 253);

    auto ta1_z_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 254);

    auto ta1_z_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 255);

    auto ta1_z_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 256);

    auto ta1_z_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 257);

    auto ta1_z_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 258);

    auto ta1_z_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 259);

    auto ta1_z_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 260);

    auto ta1_z_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 261);

    auto ta1_z_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 262);

    auto ta1_z_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 263);

    auto ta1_z_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 264);

    auto ta1_z_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 265);

    auto ta1_z_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 266);

    auto ta1_z_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 267);

    auto ta1_z_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 268);

    auto ta1_z_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 269);

    // Set up components of auxiliary buffer : GD

    auto ta1_x_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd);

    auto ta1_x_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 1);

    auto ta1_x_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 2);

    auto ta1_x_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 3);

    auto ta1_x_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 4);

    auto ta1_x_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 5);

    auto ta1_x_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 6);

    auto ta1_x_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 7);

    auto ta1_x_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 8);

    auto ta1_x_xxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 9);

    auto ta1_x_xxxy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 11);

    auto ta1_x_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 12);

    auto ta1_x_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 13);

    auto ta1_x_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 14);

    auto ta1_x_xxxz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 15);

    auto ta1_x_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 16);

    auto ta1_x_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 17);

    auto ta1_x_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 18);

    auto ta1_x_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 19);

    auto ta1_x_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 20);

    auto ta1_x_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 21);

    auto ta1_x_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 22);

    auto ta1_x_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 23);

    auto ta1_x_xxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 26);

    auto ta1_x_xxyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 29);

    auto ta1_x_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 30);

    auto ta1_x_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 31);

    auto ta1_x_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 32);

    auto ta1_x_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 33);

    auto ta1_x_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 34);

    auto ta1_x_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 35);

    auto ta1_x_xyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 36);

    auto ta1_x_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 37);

    auto ta1_x_xyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 38);

    auto ta1_x_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 39);

    auto ta1_x_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 40);

    auto ta1_x_xyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 43);

    auto ta1_x_xyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 44);

    auto ta1_x_xyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 48);

    auto ta1_x_xyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 50);

    auto ta1_x_xzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 54);

    auto ta1_x_xzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 55);

    auto ta1_x_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 56);

    auto ta1_x_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 58);

    auto ta1_x_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 59);

    auto ta1_x_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 60);

    auto ta1_x_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 61);

    auto ta1_x_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 62);

    auto ta1_x_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 63);

    auto ta1_x_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 64);

    auto ta1_x_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 65);

    auto ta1_x_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 67);

    auto ta1_x_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 68);

    auto ta1_x_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 69);

    auto ta1_x_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 70);

    auto ta1_x_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 71);

    auto ta1_x_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 72);

    auto ta1_x_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 73);

    auto ta1_x_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 74);

    auto ta1_x_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 75);

    auto ta1_x_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 76);

    auto ta1_x_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 77);

    auto ta1_x_yzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 78);

    auto ta1_x_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 80);

    auto ta1_x_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 81);

    auto ta1_x_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 82);

    auto ta1_x_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 83);

    auto ta1_x_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 84);

    auto ta1_x_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 85);

    auto ta1_x_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 86);

    auto ta1_x_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 87);

    auto ta1_x_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 88);

    auto ta1_x_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 89);

    auto ta1_y_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 90);

    auto ta1_y_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 91);

    auto ta1_y_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 92);

    auto ta1_y_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 93);

    auto ta1_y_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 94);

    auto ta1_y_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 95);

    auto ta1_y_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 96);

    auto ta1_y_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 97);

    auto ta1_y_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 98);

    auto ta1_y_xxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 99);

    auto ta1_y_xxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 100);

    auto ta1_y_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 102);

    auto ta1_y_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 103);

    auto ta1_y_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 104);

    auto ta1_y_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 106);

    auto ta1_y_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 107);

    auto ta1_y_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 108);

    auto ta1_y_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 109);

    auto ta1_y_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 110);

    auto ta1_y_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 111);

    auto ta1_y_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 112);

    auto ta1_y_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 113);

    auto ta1_y_xxyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 115);

    auto ta1_y_xxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 118);

    auto ta1_y_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 120);

    auto ta1_y_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 121);

    auto ta1_y_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 122);

    auto ta1_y_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 123);

    auto ta1_y_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 124);

    auto ta1_y_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 125);

    auto ta1_y_xyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 126);

    auto ta1_y_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 127);

    auto ta1_y_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 129);

    auto ta1_y_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 130);

    auto ta1_y_xyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 131);

    auto ta1_y_xyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 136);

    auto ta1_y_xyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 137);

    auto ta1_y_xyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 141);

    auto ta1_y_xyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 142);

    auto ta1_y_xzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 144);

    auto ta1_y_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 146);

    auto ta1_y_xzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 147);

    auto ta1_y_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 148);

    auto ta1_y_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 149);

    auto ta1_y_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 150);

    auto ta1_y_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 151);

    auto ta1_y_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 152);

    auto ta1_y_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 153);

    auto ta1_y_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 154);

    auto ta1_y_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 155);

    auto ta1_y_yyyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 156);

    auto ta1_y_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 157);

    auto ta1_y_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 158);

    auto ta1_y_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 159);

    auto ta1_y_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 160);

    auto ta1_y_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 161);

    auto ta1_y_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 162);

    auto ta1_y_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 163);

    auto ta1_y_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 164);

    auto ta1_y_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 165);

    auto ta1_y_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 166);

    auto ta1_y_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 167);

    auto ta1_y_yzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 169);

    auto ta1_y_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 170);

    auto ta1_y_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 171);

    auto ta1_y_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 172);

    auto ta1_y_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 173);

    auto ta1_y_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 174);

    auto ta1_y_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 175);

    auto ta1_y_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 176);

    auto ta1_y_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 177);

    auto ta1_y_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 178);

    auto ta1_y_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 179);

    auto ta1_z_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 180);

    auto ta1_z_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 181);

    auto ta1_z_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 182);

    auto ta1_z_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 183);

    auto ta1_z_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 184);

    auto ta1_z_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 185);

    auto ta1_z_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 186);

    auto ta1_z_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 187);

    auto ta1_z_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 188);

    auto ta1_z_xxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 189);

    auto ta1_z_xxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 190);

    auto ta1_z_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 192);

    auto ta1_z_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 193);

    auto ta1_z_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 194);

    auto ta1_z_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 196);

    auto ta1_z_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 197);

    auto ta1_z_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 198);

    auto ta1_z_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 199);

    auto ta1_z_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 200);

    auto ta1_z_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 201);

    auto ta1_z_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 202);

    auto ta1_z_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 203);

    auto ta1_z_xxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 206);

    auto ta1_z_xxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 208);

    auto ta1_z_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 210);

    auto ta1_z_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 211);

    auto ta1_z_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 212);

    auto ta1_z_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 213);

    auto ta1_z_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 214);

    auto ta1_z_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 215);

    auto ta1_z_xyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 216);

    auto ta1_z_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 217);

    auto ta1_z_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 219);

    auto ta1_z_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 220);

    auto ta1_z_xyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 221);

    auto ta1_z_xyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 226);

    auto ta1_z_xyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 227);

    auto ta1_z_xyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 231);

    auto ta1_z_xyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 232);

    auto ta1_z_xzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 234);

    auto ta1_z_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 236);

    auto ta1_z_xzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 237);

    auto ta1_z_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 238);

    auto ta1_z_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 239);

    auto ta1_z_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 240);

    auto ta1_z_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 241);

    auto ta1_z_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 242);

    auto ta1_z_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 243);

    auto ta1_z_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 244);

    auto ta1_z_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 245);

    auto ta1_z_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 247);

    auto ta1_z_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 248);

    auto ta1_z_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 249);

    auto ta1_z_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 250);

    auto ta1_z_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 251);

    auto ta1_z_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 252);

    auto ta1_z_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 253);

    auto ta1_z_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 254);

    auto ta1_z_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 255);

    auto ta1_z_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 256);

    auto ta1_z_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 257);

    auto ta1_z_yzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 258);

    auto ta1_z_yzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 259);

    auto ta1_z_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 260);

    auto ta1_z_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 261);

    auto ta1_z_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 262);

    auto ta1_z_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 263);

    auto ta1_z_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 264);

    auto ta1_z_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 265);

    auto ta1_z_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 266);

    auto ta1_z_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 267);

    auto ta1_z_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 268);

    auto ta1_z_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 269);

    // Set up 0-6 components of targeted buffer : HD

    auto ta1_x_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd);

    auto ta1_x_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 1);

    auto ta1_x_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 2);

    auto ta1_x_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 3);

    auto ta1_x_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 4);

    auto ta1_x_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 5);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_yy_0,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxx_yz_0,   \
                             ta1_x_xxx_yz_1,   \
                             ta1_x_xxx_zz_0,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_x_xxxx_x_0,   \
                             ta1_x_xxxx_x_1,   \
                             ta1_x_xxxx_xx_0,  \
                             ta1_x_xxxx_xx_1,  \
                             ta1_x_xxxx_xy_0,  \
                             ta1_x_xxxx_xy_1,  \
                             ta1_x_xxxx_xz_0,  \
                             ta1_x_xxxx_xz_1,  \
                             ta1_x_xxxx_y_0,   \
                             ta1_x_xxxx_y_1,   \
                             ta1_x_xxxx_yy_0,  \
                             ta1_x_xxxx_yy_1,  \
                             ta1_x_xxxx_yz_0,  \
                             ta1_x_xxxx_yz_1,  \
                             ta1_x_xxxx_z_0,   \
                             ta1_x_xxxx_z_1,   \
                             ta1_x_xxxx_zz_0,  \
                             ta1_x_xxxx_zz_1,  \
                             ta1_x_xxxxx_xx_0, \
                             ta1_x_xxxxx_xy_0, \
                             ta1_x_xxxxx_xz_0, \
                             ta1_x_xxxxx_yy_0, \
                             ta1_x_xxxxx_yz_0, \
                             ta1_x_xxxxx_zz_0, \
                             ta_xxxx_xx_1,     \
                             ta_xxxx_xy_1,     \
                             ta_xxxx_xz_1,     \
                             ta_xxxx_yy_1,     \
                             ta_xxxx_yz_1,     \
                             ta_xxxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxx_xx_0[i] = 4.0 * ta1_x_xxx_xx_0[i] * fe_0 - 4.0 * ta1_x_xxx_xx_1[i] * fe_0 + 2.0 * ta1_x_xxxx_x_0[i] * fe_0 -
                              2.0 * ta1_x_xxxx_x_1[i] * fe_0 + ta_xxxx_xx_1[i] + ta1_x_xxxx_xx_0[i] * pa_x[i] - ta1_x_xxxx_xx_1[i] * pc_x[i];

        ta1_x_xxxxx_xy_0[i] = 4.0 * ta1_x_xxx_xy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xy_1[i] * fe_0 + ta1_x_xxxx_y_0[i] * fe_0 - ta1_x_xxxx_y_1[i] * fe_0 +
                              ta_xxxx_xy_1[i] + ta1_x_xxxx_xy_0[i] * pa_x[i] - ta1_x_xxxx_xy_1[i] * pc_x[i];

        ta1_x_xxxxx_xz_0[i] = 4.0 * ta1_x_xxx_xz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xz_1[i] * fe_0 + ta1_x_xxxx_z_0[i] * fe_0 - ta1_x_xxxx_z_1[i] * fe_0 +
                              ta_xxxx_xz_1[i] + ta1_x_xxxx_xz_0[i] * pa_x[i] - ta1_x_xxxx_xz_1[i] * pc_x[i];

        ta1_x_xxxxx_yy_0[i] = 4.0 * ta1_x_xxx_yy_0[i] * fe_0 - 4.0 * ta1_x_xxx_yy_1[i] * fe_0 + ta_xxxx_yy_1[i] + ta1_x_xxxx_yy_0[i] * pa_x[i] -
                              ta1_x_xxxx_yy_1[i] * pc_x[i];

        ta1_x_xxxxx_yz_0[i] = 4.0 * ta1_x_xxx_yz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yz_1[i] * fe_0 + ta_xxxx_yz_1[i] + ta1_x_xxxx_yz_0[i] * pa_x[i] -
                              ta1_x_xxxx_yz_1[i] * pc_x[i];

        ta1_x_xxxxx_zz_0[i] = 4.0 * ta1_x_xxx_zz_0[i] * fe_0 - 4.0 * ta1_x_xxx_zz_1[i] * fe_0 + ta_xxxx_zz_1[i] + ta1_x_xxxx_zz_0[i] * pa_x[i] -
                              ta1_x_xxxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : HD

    auto ta1_x_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 6);

    auto ta1_x_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 7);

    auto ta1_x_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 8);

    auto ta1_x_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 9);

    auto ta1_x_xxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 10);

    auto ta1_x_xxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 11);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_xxxx_x_0,   \
                             ta1_x_xxxx_x_1,   \
                             ta1_x_xxxx_xx_0,  \
                             ta1_x_xxxx_xx_1,  \
                             ta1_x_xxxx_xy_0,  \
                             ta1_x_xxxx_xy_1,  \
                             ta1_x_xxxx_xz_0,  \
                             ta1_x_xxxx_xz_1,  \
                             ta1_x_xxxx_y_0,   \
                             ta1_x_xxxx_y_1,   \
                             ta1_x_xxxx_yy_0,  \
                             ta1_x_xxxx_yy_1,  \
                             ta1_x_xxxx_yz_0,  \
                             ta1_x_xxxx_yz_1,  \
                             ta1_x_xxxx_z_0,   \
                             ta1_x_xxxx_z_1,   \
                             ta1_x_xxxx_zz_0,  \
                             ta1_x_xxxx_zz_1,  \
                             ta1_x_xxxxy_xx_0, \
                             ta1_x_xxxxy_xy_0, \
                             ta1_x_xxxxy_xz_0, \
                             ta1_x_xxxxy_yy_0, \
                             ta1_x_xxxxy_yz_0, \
                             ta1_x_xxxxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxy_xx_0[i] = ta1_x_xxxx_xx_0[i] * pa_y[i] - ta1_x_xxxx_xx_1[i] * pc_y[i];

        ta1_x_xxxxy_xy_0[i] = ta1_x_xxxx_x_0[i] * fe_0 - ta1_x_xxxx_x_1[i] * fe_0 + ta1_x_xxxx_xy_0[i] * pa_y[i] - ta1_x_xxxx_xy_1[i] * pc_y[i];

        ta1_x_xxxxy_xz_0[i] = ta1_x_xxxx_xz_0[i] * pa_y[i] - ta1_x_xxxx_xz_1[i] * pc_y[i];

        ta1_x_xxxxy_yy_0[i] =
            2.0 * ta1_x_xxxx_y_0[i] * fe_0 - 2.0 * ta1_x_xxxx_y_1[i] * fe_0 + ta1_x_xxxx_yy_0[i] * pa_y[i] - ta1_x_xxxx_yy_1[i] * pc_y[i];

        ta1_x_xxxxy_yz_0[i] = ta1_x_xxxx_z_0[i] * fe_0 - ta1_x_xxxx_z_1[i] * fe_0 + ta1_x_xxxx_yz_0[i] * pa_y[i] - ta1_x_xxxx_yz_1[i] * pc_y[i];

        ta1_x_xxxxy_zz_0[i] = ta1_x_xxxx_zz_0[i] * pa_y[i] - ta1_x_xxxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : HD

    auto ta1_x_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 12);

    auto ta1_x_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 13);

    auto ta1_x_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 14);

    auto ta1_x_xxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 15);

    auto ta1_x_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 16);

    auto ta1_x_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 17);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_xxxx_x_0,   \
                             ta1_x_xxxx_x_1,   \
                             ta1_x_xxxx_xx_0,  \
                             ta1_x_xxxx_xx_1,  \
                             ta1_x_xxxx_xy_0,  \
                             ta1_x_xxxx_xy_1,  \
                             ta1_x_xxxx_xz_0,  \
                             ta1_x_xxxx_xz_1,  \
                             ta1_x_xxxx_y_0,   \
                             ta1_x_xxxx_y_1,   \
                             ta1_x_xxxx_yy_0,  \
                             ta1_x_xxxx_yy_1,  \
                             ta1_x_xxxx_yz_0,  \
                             ta1_x_xxxx_yz_1,  \
                             ta1_x_xxxx_z_0,   \
                             ta1_x_xxxx_z_1,   \
                             ta1_x_xxxx_zz_0,  \
                             ta1_x_xxxx_zz_1,  \
                             ta1_x_xxxxz_xx_0, \
                             ta1_x_xxxxz_xy_0, \
                             ta1_x_xxxxz_xz_0, \
                             ta1_x_xxxxz_yy_0, \
                             ta1_x_xxxxz_yz_0, \
                             ta1_x_xxxxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxz_xx_0[i] = ta1_x_xxxx_xx_0[i] * pa_z[i] - ta1_x_xxxx_xx_1[i] * pc_z[i];

        ta1_x_xxxxz_xy_0[i] = ta1_x_xxxx_xy_0[i] * pa_z[i] - ta1_x_xxxx_xy_1[i] * pc_z[i];

        ta1_x_xxxxz_xz_0[i] = ta1_x_xxxx_x_0[i] * fe_0 - ta1_x_xxxx_x_1[i] * fe_0 + ta1_x_xxxx_xz_0[i] * pa_z[i] - ta1_x_xxxx_xz_1[i] * pc_z[i];

        ta1_x_xxxxz_yy_0[i] = ta1_x_xxxx_yy_0[i] * pa_z[i] - ta1_x_xxxx_yy_1[i] * pc_z[i];

        ta1_x_xxxxz_yz_0[i] = ta1_x_xxxx_y_0[i] * fe_0 - ta1_x_xxxx_y_1[i] * fe_0 + ta1_x_xxxx_yz_0[i] * pa_z[i] - ta1_x_xxxx_yz_1[i] * pc_z[i];

        ta1_x_xxxxz_zz_0[i] =
            2.0 * ta1_x_xxxx_z_0[i] * fe_0 - 2.0 * ta1_x_xxxx_z_1[i] * fe_0 + ta1_x_xxxx_zz_0[i] * pa_z[i] - ta1_x_xxxx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : HD

    auto ta1_x_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 18);

    auto ta1_x_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 19);

    auto ta1_x_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 20);

    auto ta1_x_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 21);

    auto ta1_x_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 22);

    auto ta1_x_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 23);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_zz_0,   \
                             ta1_x_xxx_zz_1,   \
                             ta1_x_xxxy_x_0,   \
                             ta1_x_xxxy_x_1,   \
                             ta1_x_xxxy_xx_0,  \
                             ta1_x_xxxy_xx_1,  \
                             ta1_x_xxxy_xy_0,  \
                             ta1_x_xxxy_xy_1,  \
                             ta1_x_xxxy_xz_0,  \
                             ta1_x_xxxy_xz_1,  \
                             ta1_x_xxxy_zz_0,  \
                             ta1_x_xxxy_zz_1,  \
                             ta1_x_xxxyy_xx_0, \
                             ta1_x_xxxyy_xy_0, \
                             ta1_x_xxxyy_xz_0, \
                             ta1_x_xxxyy_yy_0, \
                             ta1_x_xxxyy_yz_0, \
                             ta1_x_xxxyy_zz_0, \
                             ta1_x_xxyy_yy_0,  \
                             ta1_x_xxyy_yy_1,  \
                             ta1_x_xxyy_yz_0,  \
                             ta1_x_xxyy_yz_1,  \
                             ta1_x_xyy_yy_0,   \
                             ta1_x_xyy_yy_1,   \
                             ta1_x_xyy_yz_0,   \
                             ta1_x_xyy_yz_1,   \
                             ta_xxyy_yy_1,     \
                             ta_xxyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyy_xx_0[i] = ta1_x_xxx_xx_0[i] * fe_0 - ta1_x_xxx_xx_1[i] * fe_0 + ta1_x_xxxy_xx_0[i] * pa_y[i] - ta1_x_xxxy_xx_1[i] * pc_y[i];

        ta1_x_xxxyy_xy_0[i] = ta1_x_xxx_xy_0[i] * fe_0 - ta1_x_xxx_xy_1[i] * fe_0 + ta1_x_xxxy_x_0[i] * fe_0 - ta1_x_xxxy_x_1[i] * fe_0 +
                              ta1_x_xxxy_xy_0[i] * pa_y[i] - ta1_x_xxxy_xy_1[i] * pc_y[i];

        ta1_x_xxxyy_xz_0[i] = ta1_x_xxx_xz_0[i] * fe_0 - ta1_x_xxx_xz_1[i] * fe_0 + ta1_x_xxxy_xz_0[i] * pa_y[i] - ta1_x_xxxy_xz_1[i] * pc_y[i];

        ta1_x_xxxyy_yy_0[i] = 2.0 * ta1_x_xyy_yy_0[i] * fe_0 - 2.0 * ta1_x_xyy_yy_1[i] * fe_0 + ta_xxyy_yy_1[i] + ta1_x_xxyy_yy_0[i] * pa_x[i] -
                              ta1_x_xxyy_yy_1[i] * pc_x[i];

        ta1_x_xxxyy_yz_0[i] = 2.0 * ta1_x_xyy_yz_0[i] * fe_0 - 2.0 * ta1_x_xyy_yz_1[i] * fe_0 + ta_xxyy_yz_1[i] + ta1_x_xxyy_yz_0[i] * pa_x[i] -
                              ta1_x_xxyy_yz_1[i] * pc_x[i];

        ta1_x_xxxyy_zz_0[i] = ta1_x_xxx_zz_0[i] * fe_0 - ta1_x_xxx_zz_1[i] * fe_0 + ta1_x_xxxy_zz_0[i] * pa_y[i] - ta1_x_xxxy_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : HD

    auto ta1_x_xxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 24);

    auto ta1_x_xxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 25);

    auto ta1_x_xxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 26);

    auto ta1_x_xxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 27);

    auto ta1_x_xxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 28);

    auto ta1_x_xxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 29);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxxy_xy_0,  \
                             ta1_x_xxxy_xy_1,  \
                             ta1_x_xxxy_yy_0,  \
                             ta1_x_xxxy_yy_1,  \
                             ta1_x_xxxyz_xx_0, \
                             ta1_x_xxxyz_xy_0, \
                             ta1_x_xxxyz_xz_0, \
                             ta1_x_xxxyz_yy_0, \
                             ta1_x_xxxyz_yz_0, \
                             ta1_x_xxxyz_zz_0, \
                             ta1_x_xxxz_xx_0,  \
                             ta1_x_xxxz_xx_1,  \
                             ta1_x_xxxz_xz_0,  \
                             ta1_x_xxxz_xz_1,  \
                             ta1_x_xxxz_yz_0,  \
                             ta1_x_xxxz_yz_1,  \
                             ta1_x_xxxz_z_0,   \
                             ta1_x_xxxz_z_1,   \
                             ta1_x_xxxz_zz_0,  \
                             ta1_x_xxxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyz_xx_0[i] = ta1_x_xxxz_xx_0[i] * pa_y[i] - ta1_x_xxxz_xx_1[i] * pc_y[i];

        ta1_x_xxxyz_xy_0[i] = ta1_x_xxxy_xy_0[i] * pa_z[i] - ta1_x_xxxy_xy_1[i] * pc_z[i];

        ta1_x_xxxyz_xz_0[i] = ta1_x_xxxz_xz_0[i] * pa_y[i] - ta1_x_xxxz_xz_1[i] * pc_y[i];

        ta1_x_xxxyz_yy_0[i] = ta1_x_xxxy_yy_0[i] * pa_z[i] - ta1_x_xxxy_yy_1[i] * pc_z[i];

        ta1_x_xxxyz_yz_0[i] = ta1_x_xxxz_z_0[i] * fe_0 - ta1_x_xxxz_z_1[i] * fe_0 + ta1_x_xxxz_yz_0[i] * pa_y[i] - ta1_x_xxxz_yz_1[i] * pc_y[i];

        ta1_x_xxxyz_zz_0[i] = ta1_x_xxxz_zz_0[i] * pa_y[i] - ta1_x_xxxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : HD

    auto ta1_x_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 30);

    auto ta1_x_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 31);

    auto ta1_x_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 32);

    auto ta1_x_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 33);

    auto ta1_x_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 34);

    auto ta1_x_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 35);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xxx_xx_0,   \
                             ta1_x_xxx_xx_1,   \
                             ta1_x_xxx_xy_0,   \
                             ta1_x_xxx_xy_1,   \
                             ta1_x_xxx_xz_0,   \
                             ta1_x_xxx_xz_1,   \
                             ta1_x_xxx_yy_0,   \
                             ta1_x_xxx_yy_1,   \
                             ta1_x_xxxz_x_0,   \
                             ta1_x_xxxz_x_1,   \
                             ta1_x_xxxz_xx_0,  \
                             ta1_x_xxxz_xx_1,  \
                             ta1_x_xxxz_xy_0,  \
                             ta1_x_xxxz_xy_1,  \
                             ta1_x_xxxz_xz_0,  \
                             ta1_x_xxxz_xz_1,  \
                             ta1_x_xxxz_yy_0,  \
                             ta1_x_xxxz_yy_1,  \
                             ta1_x_xxxzz_xx_0, \
                             ta1_x_xxxzz_xy_0, \
                             ta1_x_xxxzz_xz_0, \
                             ta1_x_xxxzz_yy_0, \
                             ta1_x_xxxzz_yz_0, \
                             ta1_x_xxxzz_zz_0, \
                             ta1_x_xxzz_yz_0,  \
                             ta1_x_xxzz_yz_1,  \
                             ta1_x_xxzz_zz_0,  \
                             ta1_x_xxzz_zz_1,  \
                             ta1_x_xzz_yz_0,   \
                             ta1_x_xzz_yz_1,   \
                             ta1_x_xzz_zz_0,   \
                             ta1_x_xzz_zz_1,   \
                             ta_xxzz_yz_1,     \
                             ta_xxzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzz_xx_0[i] = ta1_x_xxx_xx_0[i] * fe_0 - ta1_x_xxx_xx_1[i] * fe_0 + ta1_x_xxxz_xx_0[i] * pa_z[i] - ta1_x_xxxz_xx_1[i] * pc_z[i];

        ta1_x_xxxzz_xy_0[i] = ta1_x_xxx_xy_0[i] * fe_0 - ta1_x_xxx_xy_1[i] * fe_0 + ta1_x_xxxz_xy_0[i] * pa_z[i] - ta1_x_xxxz_xy_1[i] * pc_z[i];

        ta1_x_xxxzz_xz_0[i] = ta1_x_xxx_xz_0[i] * fe_0 - ta1_x_xxx_xz_1[i] * fe_0 + ta1_x_xxxz_x_0[i] * fe_0 - ta1_x_xxxz_x_1[i] * fe_0 +
                              ta1_x_xxxz_xz_0[i] * pa_z[i] - ta1_x_xxxz_xz_1[i] * pc_z[i];

        ta1_x_xxxzz_yy_0[i] = ta1_x_xxx_yy_0[i] * fe_0 - ta1_x_xxx_yy_1[i] * fe_0 + ta1_x_xxxz_yy_0[i] * pa_z[i] - ta1_x_xxxz_yy_1[i] * pc_z[i];

        ta1_x_xxxzz_yz_0[i] = 2.0 * ta1_x_xzz_yz_0[i] * fe_0 - 2.0 * ta1_x_xzz_yz_1[i] * fe_0 + ta_xxzz_yz_1[i] + ta1_x_xxzz_yz_0[i] * pa_x[i] -
                              ta1_x_xxzz_yz_1[i] * pc_x[i];

        ta1_x_xxxzz_zz_0[i] = 2.0 * ta1_x_xzz_zz_0[i] * fe_0 - 2.0 * ta1_x_xzz_zz_1[i] * fe_0 + ta_xxzz_zz_1[i] + ta1_x_xxzz_zz_0[i] * pa_x[i] -
                              ta1_x_xxzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : HD

    auto ta1_x_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 36);

    auto ta1_x_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 37);

    auto ta1_x_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 38);

    auto ta1_x_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 39);

    auto ta1_x_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 40);

    auto ta1_x_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 41);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xxy_xx_0,   \
                             ta1_x_xxy_xx_1,   \
                             ta1_x_xxy_xy_0,   \
                             ta1_x_xxy_xy_1,   \
                             ta1_x_xxy_xz_0,   \
                             ta1_x_xxy_xz_1,   \
                             ta1_x_xxy_zz_0,   \
                             ta1_x_xxy_zz_1,   \
                             ta1_x_xxyy_x_0,   \
                             ta1_x_xxyy_x_1,   \
                             ta1_x_xxyy_xx_0,  \
                             ta1_x_xxyy_xx_1,  \
                             ta1_x_xxyy_xy_0,  \
                             ta1_x_xxyy_xy_1,  \
                             ta1_x_xxyy_xz_0,  \
                             ta1_x_xxyy_xz_1,  \
                             ta1_x_xxyy_zz_0,  \
                             ta1_x_xxyy_zz_1,  \
                             ta1_x_xxyyy_xx_0, \
                             ta1_x_xxyyy_xy_0, \
                             ta1_x_xxyyy_xz_0, \
                             ta1_x_xxyyy_yy_0, \
                             ta1_x_xxyyy_yz_0, \
                             ta1_x_xxyyy_zz_0, \
                             ta1_x_xyyy_yy_0,  \
                             ta1_x_xyyy_yy_1,  \
                             ta1_x_xyyy_yz_0,  \
                             ta1_x_xyyy_yz_1,  \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yz_0,   \
                             ta1_x_yyy_yz_1,   \
                             ta_xyyy_yy_1,     \
                             ta_xyyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyy_xx_0[i] =
            2.0 * ta1_x_xxy_xx_0[i] * fe_0 - 2.0 * ta1_x_xxy_xx_1[i] * fe_0 + ta1_x_xxyy_xx_0[i] * pa_y[i] - ta1_x_xxyy_xx_1[i] * pc_y[i];

        ta1_x_xxyyy_xy_0[i] = 2.0 * ta1_x_xxy_xy_0[i] * fe_0 - 2.0 * ta1_x_xxy_xy_1[i] * fe_0 + ta1_x_xxyy_x_0[i] * fe_0 - ta1_x_xxyy_x_1[i] * fe_0 +
                              ta1_x_xxyy_xy_0[i] * pa_y[i] - ta1_x_xxyy_xy_1[i] * pc_y[i];

        ta1_x_xxyyy_xz_0[i] =
            2.0 * ta1_x_xxy_xz_0[i] * fe_0 - 2.0 * ta1_x_xxy_xz_1[i] * fe_0 + ta1_x_xxyy_xz_0[i] * pa_y[i] - ta1_x_xxyy_xz_1[i] * pc_y[i];

        ta1_x_xxyyy_yy_0[i] =
            ta1_x_yyy_yy_0[i] * fe_0 - ta1_x_yyy_yy_1[i] * fe_0 + ta_xyyy_yy_1[i] + ta1_x_xyyy_yy_0[i] * pa_x[i] - ta1_x_xyyy_yy_1[i] * pc_x[i];

        ta1_x_xxyyy_yz_0[i] =
            ta1_x_yyy_yz_0[i] * fe_0 - ta1_x_yyy_yz_1[i] * fe_0 + ta_xyyy_yz_1[i] + ta1_x_xyyy_yz_0[i] * pa_x[i] - ta1_x_xyyy_yz_1[i] * pc_x[i];

        ta1_x_xxyyy_zz_0[i] =
            2.0 * ta1_x_xxy_zz_0[i] * fe_0 - 2.0 * ta1_x_xxy_zz_1[i] * fe_0 + ta1_x_xxyy_zz_0[i] * pa_y[i] - ta1_x_xxyy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : HD

    auto ta1_x_xxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 42);

    auto ta1_x_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 43);

    auto ta1_x_xxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 44);

    auto ta1_x_xxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 45);

    auto ta1_x_xxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 46);

    auto ta1_x_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 47);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xxyy_xx_0,  \
                             ta1_x_xxyy_xx_1,  \
                             ta1_x_xxyy_xy_0,  \
                             ta1_x_xxyy_xy_1,  \
                             ta1_x_xxyy_y_0,   \
                             ta1_x_xxyy_y_1,   \
                             ta1_x_xxyy_yy_0,  \
                             ta1_x_xxyy_yy_1,  \
                             ta1_x_xxyy_yz_0,  \
                             ta1_x_xxyy_yz_1,  \
                             ta1_x_xxyyz_xx_0, \
                             ta1_x_xxyyz_xy_0, \
                             ta1_x_xxyyz_xz_0, \
                             ta1_x_xxyyz_yy_0, \
                             ta1_x_xxyyz_yz_0, \
                             ta1_x_xxyyz_zz_0, \
                             ta1_x_xxyz_xz_0,  \
                             ta1_x_xxyz_xz_1,  \
                             ta1_x_xxyz_zz_0,  \
                             ta1_x_xxyz_zz_1,  \
                             ta1_x_xxz_xz_0,   \
                             ta1_x_xxz_xz_1,   \
                             ta1_x_xxz_zz_0,   \
                             ta1_x_xxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyz_xx_0[i] = ta1_x_xxyy_xx_0[i] * pa_z[i] - ta1_x_xxyy_xx_1[i] * pc_z[i];

        ta1_x_xxyyz_xy_0[i] = ta1_x_xxyy_xy_0[i] * pa_z[i] - ta1_x_xxyy_xy_1[i] * pc_z[i];

        ta1_x_xxyyz_xz_0[i] = ta1_x_xxz_xz_0[i] * fe_0 - ta1_x_xxz_xz_1[i] * fe_0 + ta1_x_xxyz_xz_0[i] * pa_y[i] - ta1_x_xxyz_xz_1[i] * pc_y[i];

        ta1_x_xxyyz_yy_0[i] = ta1_x_xxyy_yy_0[i] * pa_z[i] - ta1_x_xxyy_yy_1[i] * pc_z[i];

        ta1_x_xxyyz_yz_0[i] = ta1_x_xxyy_y_0[i] * fe_0 - ta1_x_xxyy_y_1[i] * fe_0 + ta1_x_xxyy_yz_0[i] * pa_z[i] - ta1_x_xxyy_yz_1[i] * pc_z[i];

        ta1_x_xxyyz_zz_0[i] = ta1_x_xxz_zz_0[i] * fe_0 - ta1_x_xxz_zz_1[i] * fe_0 + ta1_x_xxyz_zz_0[i] * pa_y[i] - ta1_x_xxyz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : HD

    auto ta1_x_xxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 48);

    auto ta1_x_xxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 49);

    auto ta1_x_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 50);

    auto ta1_x_xxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 51);

    auto ta1_x_xxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 52);

    auto ta1_x_xxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 53);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_xxyzz_xx_0, \
                             ta1_x_xxyzz_xy_0, \
                             ta1_x_xxyzz_xz_0, \
                             ta1_x_xxyzz_yy_0, \
                             ta1_x_xxyzz_yz_0, \
                             ta1_x_xxyzz_zz_0, \
                             ta1_x_xxzz_x_0,   \
                             ta1_x_xxzz_x_1,   \
                             ta1_x_xxzz_xx_0,  \
                             ta1_x_xxzz_xx_1,  \
                             ta1_x_xxzz_xy_0,  \
                             ta1_x_xxzz_xy_1,  \
                             ta1_x_xxzz_xz_0,  \
                             ta1_x_xxzz_xz_1,  \
                             ta1_x_xxzz_y_0,   \
                             ta1_x_xxzz_y_1,   \
                             ta1_x_xxzz_yy_0,  \
                             ta1_x_xxzz_yy_1,  \
                             ta1_x_xxzz_yz_0,  \
                             ta1_x_xxzz_yz_1,  \
                             ta1_x_xxzz_z_0,   \
                             ta1_x_xxzz_z_1,   \
                             ta1_x_xxzz_zz_0,  \
                             ta1_x_xxzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzz_xx_0[i] = ta1_x_xxzz_xx_0[i] * pa_y[i] - ta1_x_xxzz_xx_1[i] * pc_y[i];

        ta1_x_xxyzz_xy_0[i] = ta1_x_xxzz_x_0[i] * fe_0 - ta1_x_xxzz_x_1[i] * fe_0 + ta1_x_xxzz_xy_0[i] * pa_y[i] - ta1_x_xxzz_xy_1[i] * pc_y[i];

        ta1_x_xxyzz_xz_0[i] = ta1_x_xxzz_xz_0[i] * pa_y[i] - ta1_x_xxzz_xz_1[i] * pc_y[i];

        ta1_x_xxyzz_yy_0[i] =
            2.0 * ta1_x_xxzz_y_0[i] * fe_0 - 2.0 * ta1_x_xxzz_y_1[i] * fe_0 + ta1_x_xxzz_yy_0[i] * pa_y[i] - ta1_x_xxzz_yy_1[i] * pc_y[i];

        ta1_x_xxyzz_yz_0[i] = ta1_x_xxzz_z_0[i] * fe_0 - ta1_x_xxzz_z_1[i] * fe_0 + ta1_x_xxzz_yz_0[i] * pa_y[i] - ta1_x_xxzz_yz_1[i] * pc_y[i];

        ta1_x_xxyzz_zz_0[i] = ta1_x_xxzz_zz_0[i] * pa_y[i] - ta1_x_xxzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : HD

    auto ta1_x_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 54);

    auto ta1_x_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 55);

    auto ta1_x_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 56);

    auto ta1_x_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 57);

    auto ta1_x_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 58);

    auto ta1_x_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 59);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xxz_xx_0,   \
                             ta1_x_xxz_xx_1,   \
                             ta1_x_xxz_xy_0,   \
                             ta1_x_xxz_xy_1,   \
                             ta1_x_xxz_xz_0,   \
                             ta1_x_xxz_xz_1,   \
                             ta1_x_xxz_yy_0,   \
                             ta1_x_xxz_yy_1,   \
                             ta1_x_xxzz_x_0,   \
                             ta1_x_xxzz_x_1,   \
                             ta1_x_xxzz_xx_0,  \
                             ta1_x_xxzz_xx_1,  \
                             ta1_x_xxzz_xy_0,  \
                             ta1_x_xxzz_xy_1,  \
                             ta1_x_xxzz_xz_0,  \
                             ta1_x_xxzz_xz_1,  \
                             ta1_x_xxzz_yy_0,  \
                             ta1_x_xxzz_yy_1,  \
                             ta1_x_xxzzz_xx_0, \
                             ta1_x_xxzzz_xy_0, \
                             ta1_x_xxzzz_xz_0, \
                             ta1_x_xxzzz_yy_0, \
                             ta1_x_xxzzz_yz_0, \
                             ta1_x_xxzzz_zz_0, \
                             ta1_x_xzzz_yz_0,  \
                             ta1_x_xzzz_yz_1,  \
                             ta1_x_xzzz_zz_0,  \
                             ta1_x_xzzz_zz_1,  \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             ta_xzzz_yz_1,     \
                             ta_xzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzz_xx_0[i] =
            2.0 * ta1_x_xxz_xx_0[i] * fe_0 - 2.0 * ta1_x_xxz_xx_1[i] * fe_0 + ta1_x_xxzz_xx_0[i] * pa_z[i] - ta1_x_xxzz_xx_1[i] * pc_z[i];

        ta1_x_xxzzz_xy_0[i] =
            2.0 * ta1_x_xxz_xy_0[i] * fe_0 - 2.0 * ta1_x_xxz_xy_1[i] * fe_0 + ta1_x_xxzz_xy_0[i] * pa_z[i] - ta1_x_xxzz_xy_1[i] * pc_z[i];

        ta1_x_xxzzz_xz_0[i] = 2.0 * ta1_x_xxz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xz_1[i] * fe_0 + ta1_x_xxzz_x_0[i] * fe_0 - ta1_x_xxzz_x_1[i] * fe_0 +
                              ta1_x_xxzz_xz_0[i] * pa_z[i] - ta1_x_xxzz_xz_1[i] * pc_z[i];

        ta1_x_xxzzz_yy_0[i] =
            2.0 * ta1_x_xxz_yy_0[i] * fe_0 - 2.0 * ta1_x_xxz_yy_1[i] * fe_0 + ta1_x_xxzz_yy_0[i] * pa_z[i] - ta1_x_xxzz_yy_1[i] * pc_z[i];

        ta1_x_xxzzz_yz_0[i] =
            ta1_x_zzz_yz_0[i] * fe_0 - ta1_x_zzz_yz_1[i] * fe_0 + ta_xzzz_yz_1[i] + ta1_x_xzzz_yz_0[i] * pa_x[i] - ta1_x_xzzz_yz_1[i] * pc_x[i];

        ta1_x_xxzzz_zz_0[i] =
            ta1_x_zzz_zz_0[i] * fe_0 - ta1_x_zzz_zz_1[i] * fe_0 + ta_xzzz_zz_1[i] + ta1_x_xzzz_zz_0[i] * pa_x[i] - ta1_x_xzzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : HD

    auto ta1_x_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 60);

    auto ta1_x_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 61);

    auto ta1_x_xyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 62);

    auto ta1_x_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 63);

    auto ta1_x_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 64);

    auto ta1_x_xyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 65);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xyy_xx_0,   \
                             ta1_x_xyy_xx_1,   \
                             ta1_x_xyy_xz_0,   \
                             ta1_x_xyy_xz_1,   \
                             ta1_x_xyyy_xx_0,  \
                             ta1_x_xyyy_xx_1,  \
                             ta1_x_xyyy_xz_0,  \
                             ta1_x_xyyy_xz_1,  \
                             ta1_x_xyyyy_xx_0, \
                             ta1_x_xyyyy_xy_0, \
                             ta1_x_xyyyy_xz_0, \
                             ta1_x_xyyyy_yy_0, \
                             ta1_x_xyyyy_yz_0, \
                             ta1_x_xyyyy_zz_0, \
                             ta1_x_yyyy_xy_0,  \
                             ta1_x_yyyy_xy_1,  \
                             ta1_x_yyyy_y_0,   \
                             ta1_x_yyyy_y_1,   \
                             ta1_x_yyyy_yy_0,  \
                             ta1_x_yyyy_yy_1,  \
                             ta1_x_yyyy_yz_0,  \
                             ta1_x_yyyy_yz_1,  \
                             ta1_x_yyyy_zz_0,  \
                             ta1_x_yyyy_zz_1,  \
                             ta_yyyy_xy_1,     \
                             ta_yyyy_yy_1,     \
                             ta_yyyy_yz_1,     \
                             ta_yyyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyy_xx_0[i] =
            3.0 * ta1_x_xyy_xx_0[i] * fe_0 - 3.0 * ta1_x_xyy_xx_1[i] * fe_0 + ta1_x_xyyy_xx_0[i] * pa_y[i] - ta1_x_xyyy_xx_1[i] * pc_y[i];

        ta1_x_xyyyy_xy_0[i] =
            ta1_x_yyyy_y_0[i] * fe_0 - ta1_x_yyyy_y_1[i] * fe_0 + ta_yyyy_xy_1[i] + ta1_x_yyyy_xy_0[i] * pa_x[i] - ta1_x_yyyy_xy_1[i] * pc_x[i];

        ta1_x_xyyyy_xz_0[i] =
            3.0 * ta1_x_xyy_xz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xz_1[i] * fe_0 + ta1_x_xyyy_xz_0[i] * pa_y[i] - ta1_x_xyyy_xz_1[i] * pc_y[i];

        ta1_x_xyyyy_yy_0[i] = ta_yyyy_yy_1[i] + ta1_x_yyyy_yy_0[i] * pa_x[i] - ta1_x_yyyy_yy_1[i] * pc_x[i];

        ta1_x_xyyyy_yz_0[i] = ta_yyyy_yz_1[i] + ta1_x_yyyy_yz_0[i] * pa_x[i] - ta1_x_yyyy_yz_1[i] * pc_x[i];

        ta1_x_xyyyy_zz_0[i] = ta_yyyy_zz_1[i] + ta1_x_yyyy_zz_0[i] * pa_x[i] - ta1_x_yyyy_zz_1[i] * pc_x[i];
    }

    // Set up 66-72 components of targeted buffer : HD

    auto ta1_x_xyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 66);

    auto ta1_x_xyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 67);

    auto ta1_x_xyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 68);

    auto ta1_x_xyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 69);

    auto ta1_x_xyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 70);

    auto ta1_x_xyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 71);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xyyy_xx_0,  \
                             ta1_x_xyyy_xx_1,  \
                             ta1_x_xyyy_xy_0,  \
                             ta1_x_xyyy_xy_1,  \
                             ta1_x_xyyy_yy_0,  \
                             ta1_x_xyyy_yy_1,  \
                             ta1_x_xyyyz_xx_0, \
                             ta1_x_xyyyz_xy_0, \
                             ta1_x_xyyyz_xz_0, \
                             ta1_x_xyyyz_yy_0, \
                             ta1_x_xyyyz_yz_0, \
                             ta1_x_xyyyz_zz_0, \
                             ta1_x_xyyz_xz_0,  \
                             ta1_x_xyyz_xz_1,  \
                             ta1_x_xyz_xz_0,   \
                             ta1_x_xyz_xz_1,   \
                             ta1_x_yyyz_yz_0,  \
                             ta1_x_yyyz_yz_1,  \
                             ta1_x_yyyz_zz_0,  \
                             ta1_x_yyyz_zz_1,  \
                             ta_yyyz_yz_1,     \
                             ta_yyyz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyz_xx_0[i] = ta1_x_xyyy_xx_0[i] * pa_z[i] - ta1_x_xyyy_xx_1[i] * pc_z[i];

        ta1_x_xyyyz_xy_0[i] = ta1_x_xyyy_xy_0[i] * pa_z[i] - ta1_x_xyyy_xy_1[i] * pc_z[i];

        ta1_x_xyyyz_xz_0[i] =
            2.0 * ta1_x_xyz_xz_0[i] * fe_0 - 2.0 * ta1_x_xyz_xz_1[i] * fe_0 + ta1_x_xyyz_xz_0[i] * pa_y[i] - ta1_x_xyyz_xz_1[i] * pc_y[i];

        ta1_x_xyyyz_yy_0[i] = ta1_x_xyyy_yy_0[i] * pa_z[i] - ta1_x_xyyy_yy_1[i] * pc_z[i];

        ta1_x_xyyyz_yz_0[i] = ta_yyyz_yz_1[i] + ta1_x_yyyz_yz_0[i] * pa_x[i] - ta1_x_yyyz_yz_1[i] * pc_x[i];

        ta1_x_xyyyz_zz_0[i] = ta_yyyz_zz_1[i] + ta1_x_yyyz_zz_0[i] * pa_x[i] - ta1_x_yyyz_zz_1[i] * pc_x[i];
    }

    // Set up 72-78 components of targeted buffer : HD

    auto ta1_x_xyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 72);

    auto ta1_x_xyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 73);

    auto ta1_x_xyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 74);

    auto ta1_x_xyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 75);

    auto ta1_x_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 76);

    auto ta1_x_xyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 77);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_xyy_xy_0,   \
                             ta1_x_xyy_xy_1,   \
                             ta1_x_xyyz_xy_0,  \
                             ta1_x_xyyz_xy_1,  \
                             ta1_x_xyyzz_xx_0, \
                             ta1_x_xyyzz_xy_0, \
                             ta1_x_xyyzz_xz_0, \
                             ta1_x_xyyzz_yy_0, \
                             ta1_x_xyyzz_yz_0, \
                             ta1_x_xyyzz_zz_0, \
                             ta1_x_xyzz_xx_0,  \
                             ta1_x_xyzz_xx_1,  \
                             ta1_x_xyzz_xz_0,  \
                             ta1_x_xyzz_xz_1,  \
                             ta1_x_xzz_xx_0,   \
                             ta1_x_xzz_xx_1,   \
                             ta1_x_xzz_xz_0,   \
                             ta1_x_xzz_xz_1,   \
                             ta1_x_yyzz_yy_0,  \
                             ta1_x_yyzz_yy_1,  \
                             ta1_x_yyzz_yz_0,  \
                             ta1_x_yyzz_yz_1,  \
                             ta1_x_yyzz_zz_0,  \
                             ta1_x_yyzz_zz_1,  \
                             ta_yyzz_yy_1,     \
                             ta_yyzz_yz_1,     \
                             ta_yyzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzz_xx_0[i] = ta1_x_xzz_xx_0[i] * fe_0 - ta1_x_xzz_xx_1[i] * fe_0 + ta1_x_xyzz_xx_0[i] * pa_y[i] - ta1_x_xyzz_xx_1[i] * pc_y[i];

        ta1_x_xyyzz_xy_0[i] = ta1_x_xyy_xy_0[i] * fe_0 - ta1_x_xyy_xy_1[i] * fe_0 + ta1_x_xyyz_xy_0[i] * pa_z[i] - ta1_x_xyyz_xy_1[i] * pc_z[i];

        ta1_x_xyyzz_xz_0[i] = ta1_x_xzz_xz_0[i] * fe_0 - ta1_x_xzz_xz_1[i] * fe_0 + ta1_x_xyzz_xz_0[i] * pa_y[i] - ta1_x_xyzz_xz_1[i] * pc_y[i];

        ta1_x_xyyzz_yy_0[i] = ta_yyzz_yy_1[i] + ta1_x_yyzz_yy_0[i] * pa_x[i] - ta1_x_yyzz_yy_1[i] * pc_x[i];

        ta1_x_xyyzz_yz_0[i] = ta_yyzz_yz_1[i] + ta1_x_yyzz_yz_0[i] * pa_x[i] - ta1_x_yyzz_yz_1[i] * pc_x[i];

        ta1_x_xyyzz_zz_0[i] = ta_yyzz_zz_1[i] + ta1_x_yyzz_zz_0[i] * pa_x[i] - ta1_x_yyzz_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : HD

    auto ta1_x_xyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 78);

    auto ta1_x_xyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 79);

    auto ta1_x_xyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 80);

    auto ta1_x_xyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 81);

    auto ta1_x_xyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 82);

    auto ta1_x_xyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 83);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_x_xyzzz_xx_0, \
                             ta1_x_xyzzz_xy_0, \
                             ta1_x_xyzzz_xz_0, \
                             ta1_x_xyzzz_yy_0, \
                             ta1_x_xyzzz_yz_0, \
                             ta1_x_xyzzz_zz_0, \
                             ta1_x_xzzz_x_0,   \
                             ta1_x_xzzz_x_1,   \
                             ta1_x_xzzz_xx_0,  \
                             ta1_x_xzzz_xx_1,  \
                             ta1_x_xzzz_xy_0,  \
                             ta1_x_xzzz_xy_1,  \
                             ta1_x_xzzz_xz_0,  \
                             ta1_x_xzzz_xz_1,  \
                             ta1_x_xzzz_zz_0,  \
                             ta1_x_xzzz_zz_1,  \
                             ta1_x_yzzz_yy_0,  \
                             ta1_x_yzzz_yy_1,  \
                             ta1_x_yzzz_yz_0,  \
                             ta1_x_yzzz_yz_1,  \
                             ta_yzzz_yy_1,     \
                             ta_yzzz_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzz_xx_0[i] = ta1_x_xzzz_xx_0[i] * pa_y[i] - ta1_x_xzzz_xx_1[i] * pc_y[i];

        ta1_x_xyzzz_xy_0[i] = ta1_x_xzzz_x_0[i] * fe_0 - ta1_x_xzzz_x_1[i] * fe_0 + ta1_x_xzzz_xy_0[i] * pa_y[i] - ta1_x_xzzz_xy_1[i] * pc_y[i];

        ta1_x_xyzzz_xz_0[i] = ta1_x_xzzz_xz_0[i] * pa_y[i] - ta1_x_xzzz_xz_1[i] * pc_y[i];

        ta1_x_xyzzz_yy_0[i] = ta_yzzz_yy_1[i] + ta1_x_yzzz_yy_0[i] * pa_x[i] - ta1_x_yzzz_yy_1[i] * pc_x[i];

        ta1_x_xyzzz_yz_0[i] = ta_yzzz_yz_1[i] + ta1_x_yzzz_yz_0[i] * pa_x[i] - ta1_x_yzzz_yz_1[i] * pc_x[i];

        ta1_x_xyzzz_zz_0[i] = ta1_x_xzzz_zz_0[i] * pa_y[i] - ta1_x_xzzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : HD

    auto ta1_x_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 84);

    auto ta1_x_xzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 85);

    auto ta1_x_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 86);

    auto ta1_x_xzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 87);

    auto ta1_x_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 88);

    auto ta1_x_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 89);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_x_xzz_xx_0,   \
                             ta1_x_xzz_xx_1,   \
                             ta1_x_xzz_xy_0,   \
                             ta1_x_xzz_xy_1,   \
                             ta1_x_xzzz_xx_0,  \
                             ta1_x_xzzz_xx_1,  \
                             ta1_x_xzzz_xy_0,  \
                             ta1_x_xzzz_xy_1,  \
                             ta1_x_xzzzz_xx_0, \
                             ta1_x_xzzzz_xy_0, \
                             ta1_x_xzzzz_xz_0, \
                             ta1_x_xzzzz_yy_0, \
                             ta1_x_xzzzz_yz_0, \
                             ta1_x_xzzzz_zz_0, \
                             ta1_x_zzzz_xz_0,  \
                             ta1_x_zzzz_xz_1,  \
                             ta1_x_zzzz_yy_0,  \
                             ta1_x_zzzz_yy_1,  \
                             ta1_x_zzzz_yz_0,  \
                             ta1_x_zzzz_yz_1,  \
                             ta1_x_zzzz_z_0,   \
                             ta1_x_zzzz_z_1,   \
                             ta1_x_zzzz_zz_0,  \
                             ta1_x_zzzz_zz_1,  \
                             ta_zzzz_xz_1,     \
                             ta_zzzz_yy_1,     \
                             ta_zzzz_yz_1,     \
                             ta_zzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzz_xx_0[i] =
            3.0 * ta1_x_xzz_xx_0[i] * fe_0 - 3.0 * ta1_x_xzz_xx_1[i] * fe_0 + ta1_x_xzzz_xx_0[i] * pa_z[i] - ta1_x_xzzz_xx_1[i] * pc_z[i];

        ta1_x_xzzzz_xy_0[i] =
            3.0 * ta1_x_xzz_xy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xy_1[i] * fe_0 + ta1_x_xzzz_xy_0[i] * pa_z[i] - ta1_x_xzzz_xy_1[i] * pc_z[i];

        ta1_x_xzzzz_xz_0[i] =
            ta1_x_zzzz_z_0[i] * fe_0 - ta1_x_zzzz_z_1[i] * fe_0 + ta_zzzz_xz_1[i] + ta1_x_zzzz_xz_0[i] * pa_x[i] - ta1_x_zzzz_xz_1[i] * pc_x[i];

        ta1_x_xzzzz_yy_0[i] = ta_zzzz_yy_1[i] + ta1_x_zzzz_yy_0[i] * pa_x[i] - ta1_x_zzzz_yy_1[i] * pc_x[i];

        ta1_x_xzzzz_yz_0[i] = ta_zzzz_yz_1[i] + ta1_x_zzzz_yz_0[i] * pa_x[i] - ta1_x_zzzz_yz_1[i] * pc_x[i];

        ta1_x_xzzzz_zz_0[i] = ta_zzzz_zz_1[i] + ta1_x_zzzz_zz_0[i] * pa_x[i] - ta1_x_zzzz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : HD

    auto ta1_x_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 90);

    auto ta1_x_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 91);

    auto ta1_x_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 92);

    auto ta1_x_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 93);

    auto ta1_x_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 94);

    auto ta1_x_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 95);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yyy_xx_0,   \
                             ta1_x_yyy_xx_1,   \
                             ta1_x_yyy_xy_0,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_xz_0,   \
                             ta1_x_yyy_xz_1,   \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyy_yz_0,   \
                             ta1_x_yyy_yz_1,   \
                             ta1_x_yyy_zz_0,   \
                             ta1_x_yyy_zz_1,   \
                             ta1_x_yyyy_x_0,   \
                             ta1_x_yyyy_x_1,   \
                             ta1_x_yyyy_xx_0,  \
                             ta1_x_yyyy_xx_1,  \
                             ta1_x_yyyy_xy_0,  \
                             ta1_x_yyyy_xy_1,  \
                             ta1_x_yyyy_xz_0,  \
                             ta1_x_yyyy_xz_1,  \
                             ta1_x_yyyy_y_0,   \
                             ta1_x_yyyy_y_1,   \
                             ta1_x_yyyy_yy_0,  \
                             ta1_x_yyyy_yy_1,  \
                             ta1_x_yyyy_yz_0,  \
                             ta1_x_yyyy_yz_1,  \
                             ta1_x_yyyy_z_0,   \
                             ta1_x_yyyy_z_1,   \
                             ta1_x_yyyy_zz_0,  \
                             ta1_x_yyyy_zz_1,  \
                             ta1_x_yyyyy_xx_0, \
                             ta1_x_yyyyy_xy_0, \
                             ta1_x_yyyyy_xz_0, \
                             ta1_x_yyyyy_yy_0, \
                             ta1_x_yyyyy_yz_0, \
                             ta1_x_yyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyy_xx_0[i] =
            4.0 * ta1_x_yyy_xx_0[i] * fe_0 - 4.0 * ta1_x_yyy_xx_1[i] * fe_0 + ta1_x_yyyy_xx_0[i] * pa_y[i] - ta1_x_yyyy_xx_1[i] * pc_y[i];

        ta1_x_yyyyy_xy_0[i] = 4.0 * ta1_x_yyy_xy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xy_1[i] * fe_0 + ta1_x_yyyy_x_0[i] * fe_0 - ta1_x_yyyy_x_1[i] * fe_0 +
                              ta1_x_yyyy_xy_0[i] * pa_y[i] - ta1_x_yyyy_xy_1[i] * pc_y[i];

        ta1_x_yyyyy_xz_0[i] =
            4.0 * ta1_x_yyy_xz_0[i] * fe_0 - 4.0 * ta1_x_yyy_xz_1[i] * fe_0 + ta1_x_yyyy_xz_0[i] * pa_y[i] - ta1_x_yyyy_xz_1[i] * pc_y[i];

        ta1_x_yyyyy_yy_0[i] = 4.0 * ta1_x_yyy_yy_0[i] * fe_0 - 4.0 * ta1_x_yyy_yy_1[i] * fe_0 + 2.0 * ta1_x_yyyy_y_0[i] * fe_0 -
                              2.0 * ta1_x_yyyy_y_1[i] * fe_0 + ta1_x_yyyy_yy_0[i] * pa_y[i] - ta1_x_yyyy_yy_1[i] * pc_y[i];

        ta1_x_yyyyy_yz_0[i] = 4.0 * ta1_x_yyy_yz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yz_1[i] * fe_0 + ta1_x_yyyy_z_0[i] * fe_0 - ta1_x_yyyy_z_1[i] * fe_0 +
                              ta1_x_yyyy_yz_0[i] * pa_y[i] - ta1_x_yyyy_yz_1[i] * pc_y[i];

        ta1_x_yyyyy_zz_0[i] =
            4.0 * ta1_x_yyy_zz_0[i] * fe_0 - 4.0 * ta1_x_yyy_zz_1[i] * fe_0 + ta1_x_yyyy_zz_0[i] * pa_y[i] - ta1_x_yyyy_zz_1[i] * pc_y[i];
    }

    // Set up 96-102 components of targeted buffer : HD

    auto ta1_x_yyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 96);

    auto ta1_x_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 97);

    auto ta1_x_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 98);

    auto ta1_x_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 99);

    auto ta1_x_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 100);

    auto ta1_x_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 101);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyyy_xx_0,  \
                             ta1_x_yyyy_xx_1,  \
                             ta1_x_yyyy_xy_0,  \
                             ta1_x_yyyy_xy_1,  \
                             ta1_x_yyyy_y_0,   \
                             ta1_x_yyyy_y_1,   \
                             ta1_x_yyyy_yy_0,  \
                             ta1_x_yyyy_yy_1,  \
                             ta1_x_yyyy_yz_0,  \
                             ta1_x_yyyy_yz_1,  \
                             ta1_x_yyyyz_xx_0, \
                             ta1_x_yyyyz_xy_0, \
                             ta1_x_yyyyz_xz_0, \
                             ta1_x_yyyyz_yy_0, \
                             ta1_x_yyyyz_yz_0, \
                             ta1_x_yyyyz_zz_0, \
                             ta1_x_yyyz_xz_0,  \
                             ta1_x_yyyz_xz_1,  \
                             ta1_x_yyyz_zz_0,  \
                             ta1_x_yyyz_zz_1,  \
                             ta1_x_yyz_xz_0,   \
                             ta1_x_yyz_xz_1,   \
                             ta1_x_yyz_zz_0,   \
                             ta1_x_yyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyz_xx_0[i] = ta1_x_yyyy_xx_0[i] * pa_z[i] - ta1_x_yyyy_xx_1[i] * pc_z[i];

        ta1_x_yyyyz_xy_0[i] = ta1_x_yyyy_xy_0[i] * pa_z[i] - ta1_x_yyyy_xy_1[i] * pc_z[i];

        ta1_x_yyyyz_xz_0[i] =
            3.0 * ta1_x_yyz_xz_0[i] * fe_0 - 3.0 * ta1_x_yyz_xz_1[i] * fe_0 + ta1_x_yyyz_xz_0[i] * pa_y[i] - ta1_x_yyyz_xz_1[i] * pc_y[i];

        ta1_x_yyyyz_yy_0[i] = ta1_x_yyyy_yy_0[i] * pa_z[i] - ta1_x_yyyy_yy_1[i] * pc_z[i];

        ta1_x_yyyyz_yz_0[i] = ta1_x_yyyy_y_0[i] * fe_0 - ta1_x_yyyy_y_1[i] * fe_0 + ta1_x_yyyy_yz_0[i] * pa_z[i] - ta1_x_yyyy_yz_1[i] * pc_z[i];

        ta1_x_yyyyz_zz_0[i] =
            3.0 * ta1_x_yyz_zz_0[i] * fe_0 - 3.0 * ta1_x_yyz_zz_1[i] * fe_0 + ta1_x_yyyz_zz_0[i] * pa_y[i] - ta1_x_yyyz_zz_1[i] * pc_y[i];
    }

    // Set up 102-108 components of targeted buffer : HD

    auto ta1_x_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 102);

    auto ta1_x_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 103);

    auto ta1_x_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 104);

    auto ta1_x_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 105);

    auto ta1_x_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 106);

    auto ta1_x_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 107);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyy_xy_0,   \
                             ta1_x_yyy_xy_1,   \
                             ta1_x_yyy_yy_0,   \
                             ta1_x_yyy_yy_1,   \
                             ta1_x_yyyz_xy_0,  \
                             ta1_x_yyyz_xy_1,  \
                             ta1_x_yyyz_yy_0,  \
                             ta1_x_yyyz_yy_1,  \
                             ta1_x_yyyzz_xx_0, \
                             ta1_x_yyyzz_xy_0, \
                             ta1_x_yyyzz_xz_0, \
                             ta1_x_yyyzz_yy_0, \
                             ta1_x_yyyzz_yz_0, \
                             ta1_x_yyyzz_zz_0, \
                             ta1_x_yyzz_xx_0,  \
                             ta1_x_yyzz_xx_1,  \
                             ta1_x_yyzz_xz_0,  \
                             ta1_x_yyzz_xz_1,  \
                             ta1_x_yyzz_yz_0,  \
                             ta1_x_yyzz_yz_1,  \
                             ta1_x_yyzz_z_0,   \
                             ta1_x_yyzz_z_1,   \
                             ta1_x_yyzz_zz_0,  \
                             ta1_x_yyzz_zz_1,  \
                             ta1_x_yzz_xx_0,   \
                             ta1_x_yzz_xx_1,   \
                             ta1_x_yzz_xz_0,   \
                             ta1_x_yzz_xz_1,   \
                             ta1_x_yzz_yz_0,   \
                             ta1_x_yzz_yz_1,   \
                             ta1_x_yzz_zz_0,   \
                             ta1_x_yzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzz_xx_0[i] =
            2.0 * ta1_x_yzz_xx_0[i] * fe_0 - 2.0 * ta1_x_yzz_xx_1[i] * fe_0 + ta1_x_yyzz_xx_0[i] * pa_y[i] - ta1_x_yyzz_xx_1[i] * pc_y[i];

        ta1_x_yyyzz_xy_0[i] = ta1_x_yyy_xy_0[i] * fe_0 - ta1_x_yyy_xy_1[i] * fe_0 + ta1_x_yyyz_xy_0[i] * pa_z[i] - ta1_x_yyyz_xy_1[i] * pc_z[i];

        ta1_x_yyyzz_xz_0[i] =
            2.0 * ta1_x_yzz_xz_0[i] * fe_0 - 2.0 * ta1_x_yzz_xz_1[i] * fe_0 + ta1_x_yyzz_xz_0[i] * pa_y[i] - ta1_x_yyzz_xz_1[i] * pc_y[i];

        ta1_x_yyyzz_yy_0[i] = ta1_x_yyy_yy_0[i] * fe_0 - ta1_x_yyy_yy_1[i] * fe_0 + ta1_x_yyyz_yy_0[i] * pa_z[i] - ta1_x_yyyz_yy_1[i] * pc_z[i];

        ta1_x_yyyzz_yz_0[i] = 2.0 * ta1_x_yzz_yz_0[i] * fe_0 - 2.0 * ta1_x_yzz_yz_1[i] * fe_0 + ta1_x_yyzz_z_0[i] * fe_0 - ta1_x_yyzz_z_1[i] * fe_0 +
                              ta1_x_yyzz_yz_0[i] * pa_y[i] - ta1_x_yyzz_yz_1[i] * pc_y[i];

        ta1_x_yyyzz_zz_0[i] =
            2.0 * ta1_x_yzz_zz_0[i] * fe_0 - 2.0 * ta1_x_yzz_zz_1[i] * fe_0 + ta1_x_yyzz_zz_0[i] * pa_y[i] - ta1_x_yyzz_zz_1[i] * pc_y[i];
    }

    // Set up 108-114 components of targeted buffer : HD

    auto ta1_x_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 108);

    auto ta1_x_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 109);

    auto ta1_x_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 110);

    auto ta1_x_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 111);

    auto ta1_x_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 112);

    auto ta1_x_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 113);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_x_yyz_xy_0,   \
                             ta1_x_yyz_xy_1,   \
                             ta1_x_yyz_yy_0,   \
                             ta1_x_yyz_yy_1,   \
                             ta1_x_yyzz_xy_0,  \
                             ta1_x_yyzz_xy_1,  \
                             ta1_x_yyzz_yy_0,  \
                             ta1_x_yyzz_yy_1,  \
                             ta1_x_yyzzz_xx_0, \
                             ta1_x_yyzzz_xy_0, \
                             ta1_x_yyzzz_xz_0, \
                             ta1_x_yyzzz_yy_0, \
                             ta1_x_yyzzz_yz_0, \
                             ta1_x_yyzzz_zz_0, \
                             ta1_x_yzzz_xx_0,  \
                             ta1_x_yzzz_xx_1,  \
                             ta1_x_yzzz_xz_0,  \
                             ta1_x_yzzz_xz_1,  \
                             ta1_x_yzzz_yz_0,  \
                             ta1_x_yzzz_yz_1,  \
                             ta1_x_yzzz_z_0,   \
                             ta1_x_yzzz_z_1,   \
                             ta1_x_yzzz_zz_0,  \
                             ta1_x_yzzz_zz_1,  \
                             ta1_x_zzz_xx_0,   \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xz_0,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzz_xx_0[i] = ta1_x_zzz_xx_0[i] * fe_0 - ta1_x_zzz_xx_1[i] * fe_0 + ta1_x_yzzz_xx_0[i] * pa_y[i] - ta1_x_yzzz_xx_1[i] * pc_y[i];

        ta1_x_yyzzz_xy_0[i] =
            2.0 * ta1_x_yyz_xy_0[i] * fe_0 - 2.0 * ta1_x_yyz_xy_1[i] * fe_0 + ta1_x_yyzz_xy_0[i] * pa_z[i] - ta1_x_yyzz_xy_1[i] * pc_z[i];

        ta1_x_yyzzz_xz_0[i] = ta1_x_zzz_xz_0[i] * fe_0 - ta1_x_zzz_xz_1[i] * fe_0 + ta1_x_yzzz_xz_0[i] * pa_y[i] - ta1_x_yzzz_xz_1[i] * pc_y[i];

        ta1_x_yyzzz_yy_0[i] =
            2.0 * ta1_x_yyz_yy_0[i] * fe_0 - 2.0 * ta1_x_yyz_yy_1[i] * fe_0 + ta1_x_yyzz_yy_0[i] * pa_z[i] - ta1_x_yyzz_yy_1[i] * pc_z[i];

        ta1_x_yyzzz_yz_0[i] = ta1_x_zzz_yz_0[i] * fe_0 - ta1_x_zzz_yz_1[i] * fe_0 + ta1_x_yzzz_z_0[i] * fe_0 - ta1_x_yzzz_z_1[i] * fe_0 +
                              ta1_x_yzzz_yz_0[i] * pa_y[i] - ta1_x_yzzz_yz_1[i] * pc_y[i];

        ta1_x_yyzzz_zz_0[i] = ta1_x_zzz_zz_0[i] * fe_0 - ta1_x_zzz_zz_1[i] * fe_0 + ta1_x_yzzz_zz_0[i] * pa_y[i] - ta1_x_yzzz_zz_1[i] * pc_y[i];
    }

    // Set up 114-120 components of targeted buffer : HD

    auto ta1_x_yzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 114);

    auto ta1_x_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 115);

    auto ta1_x_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 116);

    auto ta1_x_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 117);

    auto ta1_x_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 118);

    auto ta1_x_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 119);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_x_yzzzz_xx_0, \
                             ta1_x_yzzzz_xy_0, \
                             ta1_x_yzzzz_xz_0, \
                             ta1_x_yzzzz_yy_0, \
                             ta1_x_yzzzz_yz_0, \
                             ta1_x_yzzzz_zz_0, \
                             ta1_x_zzzz_x_0,   \
                             ta1_x_zzzz_x_1,   \
                             ta1_x_zzzz_xx_0,  \
                             ta1_x_zzzz_xx_1,  \
                             ta1_x_zzzz_xy_0,  \
                             ta1_x_zzzz_xy_1,  \
                             ta1_x_zzzz_xz_0,  \
                             ta1_x_zzzz_xz_1,  \
                             ta1_x_zzzz_y_0,   \
                             ta1_x_zzzz_y_1,   \
                             ta1_x_zzzz_yy_0,  \
                             ta1_x_zzzz_yy_1,  \
                             ta1_x_zzzz_yz_0,  \
                             ta1_x_zzzz_yz_1,  \
                             ta1_x_zzzz_z_0,   \
                             ta1_x_zzzz_z_1,   \
                             ta1_x_zzzz_zz_0,  \
                             ta1_x_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzz_xx_0[i] = ta1_x_zzzz_xx_0[i] * pa_y[i] - ta1_x_zzzz_xx_1[i] * pc_y[i];

        ta1_x_yzzzz_xy_0[i] = ta1_x_zzzz_x_0[i] * fe_0 - ta1_x_zzzz_x_1[i] * fe_0 + ta1_x_zzzz_xy_0[i] * pa_y[i] - ta1_x_zzzz_xy_1[i] * pc_y[i];

        ta1_x_yzzzz_xz_0[i] = ta1_x_zzzz_xz_0[i] * pa_y[i] - ta1_x_zzzz_xz_1[i] * pc_y[i];

        ta1_x_yzzzz_yy_0[i] =
            2.0 * ta1_x_zzzz_y_0[i] * fe_0 - 2.0 * ta1_x_zzzz_y_1[i] * fe_0 + ta1_x_zzzz_yy_0[i] * pa_y[i] - ta1_x_zzzz_yy_1[i] * pc_y[i];

        ta1_x_yzzzz_yz_0[i] = ta1_x_zzzz_z_0[i] * fe_0 - ta1_x_zzzz_z_1[i] * fe_0 + ta1_x_zzzz_yz_0[i] * pa_y[i] - ta1_x_zzzz_yz_1[i] * pc_y[i];

        ta1_x_yzzzz_zz_0[i] = ta1_x_zzzz_zz_0[i] * pa_y[i] - ta1_x_zzzz_zz_1[i] * pc_y[i];
    }

    // Set up 120-126 components of targeted buffer : HD

    auto ta1_x_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 120);

    auto ta1_x_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 121);

    auto ta1_x_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 122);

    auto ta1_x_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 123);

    auto ta1_x_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 124);

    auto ta1_x_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 125);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_x_zzz_xx_0,   \
                             ta1_x_zzz_xx_1,   \
                             ta1_x_zzz_xy_0,   \
                             ta1_x_zzz_xy_1,   \
                             ta1_x_zzz_xz_0,   \
                             ta1_x_zzz_xz_1,   \
                             ta1_x_zzz_yy_0,   \
                             ta1_x_zzz_yy_1,   \
                             ta1_x_zzz_yz_0,   \
                             ta1_x_zzz_yz_1,   \
                             ta1_x_zzz_zz_0,   \
                             ta1_x_zzz_zz_1,   \
                             ta1_x_zzzz_x_0,   \
                             ta1_x_zzzz_x_1,   \
                             ta1_x_zzzz_xx_0,  \
                             ta1_x_zzzz_xx_1,  \
                             ta1_x_zzzz_xy_0,  \
                             ta1_x_zzzz_xy_1,  \
                             ta1_x_zzzz_xz_0,  \
                             ta1_x_zzzz_xz_1,  \
                             ta1_x_zzzz_y_0,   \
                             ta1_x_zzzz_y_1,   \
                             ta1_x_zzzz_yy_0,  \
                             ta1_x_zzzz_yy_1,  \
                             ta1_x_zzzz_yz_0,  \
                             ta1_x_zzzz_yz_1,  \
                             ta1_x_zzzz_z_0,   \
                             ta1_x_zzzz_z_1,   \
                             ta1_x_zzzz_zz_0,  \
                             ta1_x_zzzz_zz_1,  \
                             ta1_x_zzzzz_xx_0, \
                             ta1_x_zzzzz_xy_0, \
                             ta1_x_zzzzz_xz_0, \
                             ta1_x_zzzzz_yy_0, \
                             ta1_x_zzzzz_yz_0, \
                             ta1_x_zzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzz_xx_0[i] =
            4.0 * ta1_x_zzz_xx_0[i] * fe_0 - 4.0 * ta1_x_zzz_xx_1[i] * fe_0 + ta1_x_zzzz_xx_0[i] * pa_z[i] - ta1_x_zzzz_xx_1[i] * pc_z[i];

        ta1_x_zzzzz_xy_0[i] =
            4.0 * ta1_x_zzz_xy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xy_1[i] * fe_0 + ta1_x_zzzz_xy_0[i] * pa_z[i] - ta1_x_zzzz_xy_1[i] * pc_z[i];

        ta1_x_zzzzz_xz_0[i] = 4.0 * ta1_x_zzz_xz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xz_1[i] * fe_0 + ta1_x_zzzz_x_0[i] * fe_0 - ta1_x_zzzz_x_1[i] * fe_0 +
                              ta1_x_zzzz_xz_0[i] * pa_z[i] - ta1_x_zzzz_xz_1[i] * pc_z[i];

        ta1_x_zzzzz_yy_0[i] =
            4.0 * ta1_x_zzz_yy_0[i] * fe_0 - 4.0 * ta1_x_zzz_yy_1[i] * fe_0 + ta1_x_zzzz_yy_0[i] * pa_z[i] - ta1_x_zzzz_yy_1[i] * pc_z[i];

        ta1_x_zzzzz_yz_0[i] = 4.0 * ta1_x_zzz_yz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yz_1[i] * fe_0 + ta1_x_zzzz_y_0[i] * fe_0 - ta1_x_zzzz_y_1[i] * fe_0 +
                              ta1_x_zzzz_yz_0[i] * pa_z[i] - ta1_x_zzzz_yz_1[i] * pc_z[i];

        ta1_x_zzzzz_zz_0[i] = 4.0 * ta1_x_zzz_zz_0[i] * fe_0 - 4.0 * ta1_x_zzz_zz_1[i] * fe_0 + 2.0 * ta1_x_zzzz_z_0[i] * fe_0 -
                              2.0 * ta1_x_zzzz_z_1[i] * fe_0 + ta1_x_zzzz_zz_0[i] * pa_z[i] - ta1_x_zzzz_zz_1[i] * pc_z[i];
    }

    // Set up 126-132 components of targeted buffer : HD

    auto ta1_y_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 126);

    auto ta1_y_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 127);

    auto ta1_y_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 128);

    auto ta1_y_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 129);

    auto ta1_y_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 130);

    auto ta1_y_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 131);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xy_0,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxx_xz_0,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxx_yy_0,   \
                             ta1_y_xxx_yy_1,   \
                             ta1_y_xxx_yz_0,   \
                             ta1_y_xxx_yz_1,   \
                             ta1_y_xxx_zz_0,   \
                             ta1_y_xxx_zz_1,   \
                             ta1_y_xxxx_x_0,   \
                             ta1_y_xxxx_x_1,   \
                             ta1_y_xxxx_xx_0,  \
                             ta1_y_xxxx_xx_1,  \
                             ta1_y_xxxx_xy_0,  \
                             ta1_y_xxxx_xy_1,  \
                             ta1_y_xxxx_xz_0,  \
                             ta1_y_xxxx_xz_1,  \
                             ta1_y_xxxx_y_0,   \
                             ta1_y_xxxx_y_1,   \
                             ta1_y_xxxx_yy_0,  \
                             ta1_y_xxxx_yy_1,  \
                             ta1_y_xxxx_yz_0,  \
                             ta1_y_xxxx_yz_1,  \
                             ta1_y_xxxx_z_0,   \
                             ta1_y_xxxx_z_1,   \
                             ta1_y_xxxx_zz_0,  \
                             ta1_y_xxxx_zz_1,  \
                             ta1_y_xxxxx_xx_0, \
                             ta1_y_xxxxx_xy_0, \
                             ta1_y_xxxxx_xz_0, \
                             ta1_y_xxxxx_yy_0, \
                             ta1_y_xxxxx_yz_0, \
                             ta1_y_xxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxx_xx_0[i] = 4.0 * ta1_y_xxx_xx_0[i] * fe_0 - 4.0 * ta1_y_xxx_xx_1[i] * fe_0 + 2.0 * ta1_y_xxxx_x_0[i] * fe_0 -
                              2.0 * ta1_y_xxxx_x_1[i] * fe_0 + ta1_y_xxxx_xx_0[i] * pa_x[i] - ta1_y_xxxx_xx_1[i] * pc_x[i];

        ta1_y_xxxxx_xy_0[i] = 4.0 * ta1_y_xxx_xy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xy_1[i] * fe_0 + ta1_y_xxxx_y_0[i] * fe_0 - ta1_y_xxxx_y_1[i] * fe_0 +
                              ta1_y_xxxx_xy_0[i] * pa_x[i] - ta1_y_xxxx_xy_1[i] * pc_x[i];

        ta1_y_xxxxx_xz_0[i] = 4.0 * ta1_y_xxx_xz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xz_1[i] * fe_0 + ta1_y_xxxx_z_0[i] * fe_0 - ta1_y_xxxx_z_1[i] * fe_0 +
                              ta1_y_xxxx_xz_0[i] * pa_x[i] - ta1_y_xxxx_xz_1[i] * pc_x[i];

        ta1_y_xxxxx_yy_0[i] =
            4.0 * ta1_y_xxx_yy_0[i] * fe_0 - 4.0 * ta1_y_xxx_yy_1[i] * fe_0 + ta1_y_xxxx_yy_0[i] * pa_x[i] - ta1_y_xxxx_yy_1[i] * pc_x[i];

        ta1_y_xxxxx_yz_0[i] =
            4.0 * ta1_y_xxx_yz_0[i] * fe_0 - 4.0 * ta1_y_xxx_yz_1[i] * fe_0 + ta1_y_xxxx_yz_0[i] * pa_x[i] - ta1_y_xxxx_yz_1[i] * pc_x[i];

        ta1_y_xxxxx_zz_0[i] =
            4.0 * ta1_y_xxx_zz_0[i] * fe_0 - 4.0 * ta1_y_xxx_zz_1[i] * fe_0 + ta1_y_xxxx_zz_0[i] * pa_x[i] - ta1_y_xxxx_zz_1[i] * pc_x[i];
    }

    // Set up 132-138 components of targeted buffer : HD

    auto ta1_y_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 132);

    auto ta1_y_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 133);

    auto ta1_y_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 134);

    auto ta1_y_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 135);

    auto ta1_y_xxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 136);

    auto ta1_y_xxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 137);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxxx_x_0,   \
                             ta1_y_xxxx_x_1,   \
                             ta1_y_xxxx_xx_0,  \
                             ta1_y_xxxx_xx_1,  \
                             ta1_y_xxxx_xy_0,  \
                             ta1_y_xxxx_xy_1,  \
                             ta1_y_xxxx_xz_0,  \
                             ta1_y_xxxx_xz_1,  \
                             ta1_y_xxxx_zz_0,  \
                             ta1_y_xxxx_zz_1,  \
                             ta1_y_xxxxy_xx_0, \
                             ta1_y_xxxxy_xy_0, \
                             ta1_y_xxxxy_xz_0, \
                             ta1_y_xxxxy_yy_0, \
                             ta1_y_xxxxy_yz_0, \
                             ta1_y_xxxxy_zz_0, \
                             ta1_y_xxxy_yy_0,  \
                             ta1_y_xxxy_yy_1,  \
                             ta1_y_xxxy_yz_0,  \
                             ta1_y_xxxy_yz_1,  \
                             ta1_y_xxy_yy_0,   \
                             ta1_y_xxy_yy_1,   \
                             ta1_y_xxy_yz_0,   \
                             ta1_y_xxy_yz_1,   \
                             ta_xxxx_xx_1,     \
                             ta_xxxx_xy_1,     \
                             ta_xxxx_xz_1,     \
                             ta_xxxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxy_xx_0[i] = ta_xxxx_xx_1[i] + ta1_y_xxxx_xx_0[i] * pa_y[i] - ta1_y_xxxx_xx_1[i] * pc_y[i];

        ta1_y_xxxxy_xy_0[i] =
            ta1_y_xxxx_x_0[i] * fe_0 - ta1_y_xxxx_x_1[i] * fe_0 + ta_xxxx_xy_1[i] + ta1_y_xxxx_xy_0[i] * pa_y[i] - ta1_y_xxxx_xy_1[i] * pc_y[i];

        ta1_y_xxxxy_xz_0[i] = ta_xxxx_xz_1[i] + ta1_y_xxxx_xz_0[i] * pa_y[i] - ta1_y_xxxx_xz_1[i] * pc_y[i];

        ta1_y_xxxxy_yy_0[i] =
            3.0 * ta1_y_xxy_yy_0[i] * fe_0 - 3.0 * ta1_y_xxy_yy_1[i] * fe_0 + ta1_y_xxxy_yy_0[i] * pa_x[i] - ta1_y_xxxy_yy_1[i] * pc_x[i];

        ta1_y_xxxxy_yz_0[i] =
            3.0 * ta1_y_xxy_yz_0[i] * fe_0 - 3.0 * ta1_y_xxy_yz_1[i] * fe_0 + ta1_y_xxxy_yz_0[i] * pa_x[i] - ta1_y_xxxy_yz_1[i] * pc_x[i];

        ta1_y_xxxxy_zz_0[i] = ta_xxxx_zz_1[i] + ta1_y_xxxx_zz_0[i] * pa_y[i] - ta1_y_xxxx_zz_1[i] * pc_y[i];
    }

    // Set up 138-144 components of targeted buffer : HD

    auto ta1_y_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 138);

    auto ta1_y_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 139);

    auto ta1_y_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 140);

    auto ta1_y_xxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 141);

    auto ta1_y_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 142);

    auto ta1_y_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 143);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxxx_x_0,   \
                             ta1_y_xxxx_x_1,   \
                             ta1_y_xxxx_xx_0,  \
                             ta1_y_xxxx_xx_1,  \
                             ta1_y_xxxx_xy_0,  \
                             ta1_y_xxxx_xy_1,  \
                             ta1_y_xxxx_xz_0,  \
                             ta1_y_xxxx_xz_1,  \
                             ta1_y_xxxx_yy_0,  \
                             ta1_y_xxxx_yy_1,  \
                             ta1_y_xxxxz_xx_0, \
                             ta1_y_xxxxz_xy_0, \
                             ta1_y_xxxxz_xz_0, \
                             ta1_y_xxxxz_yy_0, \
                             ta1_y_xxxxz_yz_0, \
                             ta1_y_xxxxz_zz_0, \
                             ta1_y_xxxz_yz_0,  \
                             ta1_y_xxxz_yz_1,  \
                             ta1_y_xxxz_zz_0,  \
                             ta1_y_xxxz_zz_1,  \
                             ta1_y_xxz_yz_0,   \
                             ta1_y_xxz_yz_1,   \
                             ta1_y_xxz_zz_0,   \
                             ta1_y_xxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxz_xx_0[i] = ta1_y_xxxx_xx_0[i] * pa_z[i] - ta1_y_xxxx_xx_1[i] * pc_z[i];

        ta1_y_xxxxz_xy_0[i] = ta1_y_xxxx_xy_0[i] * pa_z[i] - ta1_y_xxxx_xy_1[i] * pc_z[i];

        ta1_y_xxxxz_xz_0[i] = ta1_y_xxxx_x_0[i] * fe_0 - ta1_y_xxxx_x_1[i] * fe_0 + ta1_y_xxxx_xz_0[i] * pa_z[i] - ta1_y_xxxx_xz_1[i] * pc_z[i];

        ta1_y_xxxxz_yy_0[i] = ta1_y_xxxx_yy_0[i] * pa_z[i] - ta1_y_xxxx_yy_1[i] * pc_z[i];

        ta1_y_xxxxz_yz_0[i] =
            3.0 * ta1_y_xxz_yz_0[i] * fe_0 - 3.0 * ta1_y_xxz_yz_1[i] * fe_0 + ta1_y_xxxz_yz_0[i] * pa_x[i] - ta1_y_xxxz_yz_1[i] * pc_x[i];

        ta1_y_xxxxz_zz_0[i] =
            3.0 * ta1_y_xxz_zz_0[i] * fe_0 - 3.0 * ta1_y_xxz_zz_1[i] * fe_0 + ta1_y_xxxz_zz_0[i] * pa_x[i] - ta1_y_xxxz_zz_1[i] * pc_x[i];
    }

    // Set up 144-150 components of targeted buffer : HD

    auto ta1_y_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 144);

    auto ta1_y_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 145);

    auto ta1_y_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 146);

    auto ta1_y_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 147);

    auto ta1_y_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 148);

    auto ta1_y_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 149);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xz_0,   \
                             ta1_y_xxx_xz_1,   \
                             ta1_y_xxxy_xx_0,  \
                             ta1_y_xxxy_xx_1,  \
                             ta1_y_xxxy_xz_0,  \
                             ta1_y_xxxy_xz_1,  \
                             ta1_y_xxxyy_xx_0, \
                             ta1_y_xxxyy_xy_0, \
                             ta1_y_xxxyy_xz_0, \
                             ta1_y_xxxyy_yy_0, \
                             ta1_y_xxxyy_yz_0, \
                             ta1_y_xxxyy_zz_0, \
                             ta1_y_xxyy_xy_0,  \
                             ta1_y_xxyy_xy_1,  \
                             ta1_y_xxyy_y_0,   \
                             ta1_y_xxyy_y_1,   \
                             ta1_y_xxyy_yy_0,  \
                             ta1_y_xxyy_yy_1,  \
                             ta1_y_xxyy_yz_0,  \
                             ta1_y_xxyy_yz_1,  \
                             ta1_y_xxyy_zz_0,  \
                             ta1_y_xxyy_zz_1,  \
                             ta1_y_xyy_xy_0,   \
                             ta1_y_xyy_xy_1,   \
                             ta1_y_xyy_yy_0,   \
                             ta1_y_xyy_yy_1,   \
                             ta1_y_xyy_yz_0,   \
                             ta1_y_xyy_yz_1,   \
                             ta1_y_xyy_zz_0,   \
                             ta1_y_xyy_zz_1,   \
                             ta_xxxy_xx_1,     \
                             ta_xxxy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyy_xx_0[i] =
            ta1_y_xxx_xx_0[i] * fe_0 - ta1_y_xxx_xx_1[i] * fe_0 + ta_xxxy_xx_1[i] + ta1_y_xxxy_xx_0[i] * pa_y[i] - ta1_y_xxxy_xx_1[i] * pc_y[i];

        ta1_y_xxxyy_xy_0[i] = 2.0 * ta1_y_xyy_xy_0[i] * fe_0 - 2.0 * ta1_y_xyy_xy_1[i] * fe_0 + ta1_y_xxyy_y_0[i] * fe_0 - ta1_y_xxyy_y_1[i] * fe_0 +
                              ta1_y_xxyy_xy_0[i] * pa_x[i] - ta1_y_xxyy_xy_1[i] * pc_x[i];

        ta1_y_xxxyy_xz_0[i] =
            ta1_y_xxx_xz_0[i] * fe_0 - ta1_y_xxx_xz_1[i] * fe_0 + ta_xxxy_xz_1[i] + ta1_y_xxxy_xz_0[i] * pa_y[i] - ta1_y_xxxy_xz_1[i] * pc_y[i];

        ta1_y_xxxyy_yy_0[i] =
            2.0 * ta1_y_xyy_yy_0[i] * fe_0 - 2.0 * ta1_y_xyy_yy_1[i] * fe_0 + ta1_y_xxyy_yy_0[i] * pa_x[i] - ta1_y_xxyy_yy_1[i] * pc_x[i];

        ta1_y_xxxyy_yz_0[i] =
            2.0 * ta1_y_xyy_yz_0[i] * fe_0 - 2.0 * ta1_y_xyy_yz_1[i] * fe_0 + ta1_y_xxyy_yz_0[i] * pa_x[i] - ta1_y_xxyy_yz_1[i] * pc_x[i];

        ta1_y_xxxyy_zz_0[i] =
            2.0 * ta1_y_xyy_zz_0[i] * fe_0 - 2.0 * ta1_y_xyy_zz_1[i] * fe_0 + ta1_y_xxyy_zz_0[i] * pa_x[i] - ta1_y_xxyy_zz_1[i] * pc_x[i];
    }

    // Set up 150-156 components of targeted buffer : HD

    auto ta1_y_xxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 150);

    auto ta1_y_xxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 151);

    auto ta1_y_xxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 152);

    auto ta1_y_xxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 153);

    auto ta1_y_xxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 154);

    auto ta1_y_xxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 155);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xxxy_xx_0,  \
                             ta1_y_xxxy_xx_1,  \
                             ta1_y_xxxy_xy_0,  \
                             ta1_y_xxxy_xy_1,  \
                             ta1_y_xxxy_yy_0,  \
                             ta1_y_xxxy_yy_1,  \
                             ta1_y_xxxyz_xx_0, \
                             ta1_y_xxxyz_xy_0, \
                             ta1_y_xxxyz_xz_0, \
                             ta1_y_xxxyz_yy_0, \
                             ta1_y_xxxyz_yz_0, \
                             ta1_y_xxxyz_zz_0, \
                             ta1_y_xxxz_xz_0,  \
                             ta1_y_xxxz_xz_1,  \
                             ta1_y_xxxz_zz_0,  \
                             ta1_y_xxxz_zz_1,  \
                             ta1_y_xxyz_yz_0,  \
                             ta1_y_xxyz_yz_1,  \
                             ta1_y_xyz_yz_0,   \
                             ta1_y_xyz_yz_1,   \
                             ta_xxxz_xz_1,     \
                             ta_xxxz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyz_xx_0[i] = ta1_y_xxxy_xx_0[i] * pa_z[i] - ta1_y_xxxy_xx_1[i] * pc_z[i];

        ta1_y_xxxyz_xy_0[i] = ta1_y_xxxy_xy_0[i] * pa_z[i] - ta1_y_xxxy_xy_1[i] * pc_z[i];

        ta1_y_xxxyz_xz_0[i] = ta_xxxz_xz_1[i] + ta1_y_xxxz_xz_0[i] * pa_y[i] - ta1_y_xxxz_xz_1[i] * pc_y[i];

        ta1_y_xxxyz_yy_0[i] = ta1_y_xxxy_yy_0[i] * pa_z[i] - ta1_y_xxxy_yy_1[i] * pc_z[i];

        ta1_y_xxxyz_yz_0[i] =
            2.0 * ta1_y_xyz_yz_0[i] * fe_0 - 2.0 * ta1_y_xyz_yz_1[i] * fe_0 + ta1_y_xxyz_yz_0[i] * pa_x[i] - ta1_y_xxyz_yz_1[i] * pc_x[i];

        ta1_y_xxxyz_zz_0[i] = ta_xxxz_zz_1[i] + ta1_y_xxxz_zz_0[i] * pa_y[i] - ta1_y_xxxz_zz_1[i] * pc_y[i];
    }

    // Set up 156-162 components of targeted buffer : HD

    auto ta1_y_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 156);

    auto ta1_y_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 157);

    auto ta1_y_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 158);

    auto ta1_y_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 159);

    auto ta1_y_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 160);

    auto ta1_y_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 161);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxx_xx_0,   \
                             ta1_y_xxx_xx_1,   \
                             ta1_y_xxx_xy_0,   \
                             ta1_y_xxx_xy_1,   \
                             ta1_y_xxxz_xx_0,  \
                             ta1_y_xxxz_xx_1,  \
                             ta1_y_xxxz_xy_0,  \
                             ta1_y_xxxz_xy_1,  \
                             ta1_y_xxxzz_xx_0, \
                             ta1_y_xxxzz_xy_0, \
                             ta1_y_xxxzz_xz_0, \
                             ta1_y_xxxzz_yy_0, \
                             ta1_y_xxxzz_yz_0, \
                             ta1_y_xxxzz_zz_0, \
                             ta1_y_xxzz_xz_0,  \
                             ta1_y_xxzz_xz_1,  \
                             ta1_y_xxzz_yy_0,  \
                             ta1_y_xxzz_yy_1,  \
                             ta1_y_xxzz_yz_0,  \
                             ta1_y_xxzz_yz_1,  \
                             ta1_y_xxzz_z_0,   \
                             ta1_y_xxzz_z_1,   \
                             ta1_y_xxzz_zz_0,  \
                             ta1_y_xxzz_zz_1,  \
                             ta1_y_xzz_xz_0,   \
                             ta1_y_xzz_xz_1,   \
                             ta1_y_xzz_yy_0,   \
                             ta1_y_xzz_yy_1,   \
                             ta1_y_xzz_yz_0,   \
                             ta1_y_xzz_yz_1,   \
                             ta1_y_xzz_zz_0,   \
                             ta1_y_xzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzz_xx_0[i] = ta1_y_xxx_xx_0[i] * fe_0 - ta1_y_xxx_xx_1[i] * fe_0 + ta1_y_xxxz_xx_0[i] * pa_z[i] - ta1_y_xxxz_xx_1[i] * pc_z[i];

        ta1_y_xxxzz_xy_0[i] = ta1_y_xxx_xy_0[i] * fe_0 - ta1_y_xxx_xy_1[i] * fe_0 + ta1_y_xxxz_xy_0[i] * pa_z[i] - ta1_y_xxxz_xy_1[i] * pc_z[i];

        ta1_y_xxxzz_xz_0[i] = 2.0 * ta1_y_xzz_xz_0[i] * fe_0 - 2.0 * ta1_y_xzz_xz_1[i] * fe_0 + ta1_y_xxzz_z_0[i] * fe_0 - ta1_y_xxzz_z_1[i] * fe_0 +
                              ta1_y_xxzz_xz_0[i] * pa_x[i] - ta1_y_xxzz_xz_1[i] * pc_x[i];

        ta1_y_xxxzz_yy_0[i] =
            2.0 * ta1_y_xzz_yy_0[i] * fe_0 - 2.0 * ta1_y_xzz_yy_1[i] * fe_0 + ta1_y_xxzz_yy_0[i] * pa_x[i] - ta1_y_xxzz_yy_1[i] * pc_x[i];

        ta1_y_xxxzz_yz_0[i] =
            2.0 * ta1_y_xzz_yz_0[i] * fe_0 - 2.0 * ta1_y_xzz_yz_1[i] * fe_0 + ta1_y_xxzz_yz_0[i] * pa_x[i] - ta1_y_xxzz_yz_1[i] * pc_x[i];

        ta1_y_xxxzz_zz_0[i] =
            2.0 * ta1_y_xzz_zz_0[i] * fe_0 - 2.0 * ta1_y_xzz_zz_1[i] * fe_0 + ta1_y_xxzz_zz_0[i] * pa_x[i] - ta1_y_xxzz_zz_1[i] * pc_x[i];
    }

    // Set up 162-168 components of targeted buffer : HD

    auto ta1_y_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 162);

    auto ta1_y_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 163);

    auto ta1_y_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 164);

    auto ta1_y_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 165);

    auto ta1_y_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 166);

    auto ta1_y_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 167);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xxy_xx_0,   \
                             ta1_y_xxy_xx_1,   \
                             ta1_y_xxy_xz_0,   \
                             ta1_y_xxy_xz_1,   \
                             ta1_y_xxyy_xx_0,  \
                             ta1_y_xxyy_xx_1,  \
                             ta1_y_xxyy_xz_0,  \
                             ta1_y_xxyy_xz_1,  \
                             ta1_y_xxyyy_xx_0, \
                             ta1_y_xxyyy_xy_0, \
                             ta1_y_xxyyy_xz_0, \
                             ta1_y_xxyyy_yy_0, \
                             ta1_y_xxyyy_yz_0, \
                             ta1_y_xxyyy_zz_0, \
                             ta1_y_xyyy_xy_0,  \
                             ta1_y_xyyy_xy_1,  \
                             ta1_y_xyyy_y_0,   \
                             ta1_y_xyyy_y_1,   \
                             ta1_y_xyyy_yy_0,  \
                             ta1_y_xyyy_yy_1,  \
                             ta1_y_xyyy_yz_0,  \
                             ta1_y_xyyy_yz_1,  \
                             ta1_y_xyyy_zz_0,  \
                             ta1_y_xyyy_zz_1,  \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_zz_0,   \
                             ta1_y_yyy_zz_1,   \
                             ta_xxyy_xx_1,     \
                             ta_xxyy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyy_xx_0[i] = 2.0 * ta1_y_xxy_xx_0[i] * fe_0 - 2.0 * ta1_y_xxy_xx_1[i] * fe_0 + ta_xxyy_xx_1[i] + ta1_y_xxyy_xx_0[i] * pa_y[i] -
                              ta1_y_xxyy_xx_1[i] * pc_y[i];

        ta1_y_xxyyy_xy_0[i] = ta1_y_yyy_xy_0[i] * fe_0 - ta1_y_yyy_xy_1[i] * fe_0 + ta1_y_xyyy_y_0[i] * fe_0 - ta1_y_xyyy_y_1[i] * fe_0 +
                              ta1_y_xyyy_xy_0[i] * pa_x[i] - ta1_y_xyyy_xy_1[i] * pc_x[i];

        ta1_y_xxyyy_xz_0[i] = 2.0 * ta1_y_xxy_xz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xz_1[i] * fe_0 + ta_xxyy_xz_1[i] + ta1_y_xxyy_xz_0[i] * pa_y[i] -
                              ta1_y_xxyy_xz_1[i] * pc_y[i];

        ta1_y_xxyyy_yy_0[i] = ta1_y_yyy_yy_0[i] * fe_0 - ta1_y_yyy_yy_1[i] * fe_0 + ta1_y_xyyy_yy_0[i] * pa_x[i] - ta1_y_xyyy_yy_1[i] * pc_x[i];

        ta1_y_xxyyy_yz_0[i] = ta1_y_yyy_yz_0[i] * fe_0 - ta1_y_yyy_yz_1[i] * fe_0 + ta1_y_xyyy_yz_0[i] * pa_x[i] - ta1_y_xyyy_yz_1[i] * pc_x[i];

        ta1_y_xxyyy_zz_0[i] = ta1_y_yyy_zz_0[i] * fe_0 - ta1_y_yyy_zz_1[i] * fe_0 + ta1_y_xyyy_zz_0[i] * pa_x[i] - ta1_y_xyyy_zz_1[i] * pc_x[i];
    }

    // Set up 168-174 components of targeted buffer : HD

    auto ta1_y_xxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 168);

    auto ta1_y_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 169);

    auto ta1_y_xxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 170);

    auto ta1_y_xxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 171);

    auto ta1_y_xxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 172);

    auto ta1_y_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 173);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxyy_x_0,   \
                             ta1_y_xxyy_x_1,   \
                             ta1_y_xxyy_xx_0,  \
                             ta1_y_xxyy_xx_1,  \
                             ta1_y_xxyy_xy_0,  \
                             ta1_y_xxyy_xy_1,  \
                             ta1_y_xxyy_xz_0,  \
                             ta1_y_xxyy_xz_1,  \
                             ta1_y_xxyy_yy_0,  \
                             ta1_y_xxyy_yy_1,  \
                             ta1_y_xxyyz_xx_0, \
                             ta1_y_xxyyz_xy_0, \
                             ta1_y_xxyyz_xz_0, \
                             ta1_y_xxyyz_yy_0, \
                             ta1_y_xxyyz_yz_0, \
                             ta1_y_xxyyz_zz_0, \
                             ta1_y_xyyz_yz_0,  \
                             ta1_y_xyyz_yz_1,  \
                             ta1_y_xyyz_zz_0,  \
                             ta1_y_xyyz_zz_1,  \
                             ta1_y_yyz_yz_0,   \
                             ta1_y_yyz_yz_1,   \
                             ta1_y_yyz_zz_0,   \
                             ta1_y_yyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyz_xx_0[i] = ta1_y_xxyy_xx_0[i] * pa_z[i] - ta1_y_xxyy_xx_1[i] * pc_z[i];

        ta1_y_xxyyz_xy_0[i] = ta1_y_xxyy_xy_0[i] * pa_z[i] - ta1_y_xxyy_xy_1[i] * pc_z[i];

        ta1_y_xxyyz_xz_0[i] = ta1_y_xxyy_x_0[i] * fe_0 - ta1_y_xxyy_x_1[i] * fe_0 + ta1_y_xxyy_xz_0[i] * pa_z[i] - ta1_y_xxyy_xz_1[i] * pc_z[i];

        ta1_y_xxyyz_yy_0[i] = ta1_y_xxyy_yy_0[i] * pa_z[i] - ta1_y_xxyy_yy_1[i] * pc_z[i];

        ta1_y_xxyyz_yz_0[i] = ta1_y_yyz_yz_0[i] * fe_0 - ta1_y_yyz_yz_1[i] * fe_0 + ta1_y_xyyz_yz_0[i] * pa_x[i] - ta1_y_xyyz_yz_1[i] * pc_x[i];

        ta1_y_xxyyz_zz_0[i] = ta1_y_yyz_zz_0[i] * fe_0 - ta1_y_yyz_zz_1[i] * fe_0 + ta1_y_xyyz_zz_0[i] * pa_x[i] - ta1_y_xyyz_zz_1[i] * pc_x[i];
    }

    // Set up 174-180 components of targeted buffer : HD

    auto ta1_y_xxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 174);

    auto ta1_y_xxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 175);

    auto ta1_y_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 176);

    auto ta1_y_xxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 177);

    auto ta1_y_xxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 178);

    auto ta1_y_xxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 179);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_xxy_xy_0,   \
                             ta1_y_xxy_xy_1,   \
                             ta1_y_xxyz_xy_0,  \
                             ta1_y_xxyz_xy_1,  \
                             ta1_y_xxyzz_xx_0, \
                             ta1_y_xxyzz_xy_0, \
                             ta1_y_xxyzz_xz_0, \
                             ta1_y_xxyzz_yy_0, \
                             ta1_y_xxyzz_yz_0, \
                             ta1_y_xxyzz_zz_0, \
                             ta1_y_xxzz_xx_0,  \
                             ta1_y_xxzz_xx_1,  \
                             ta1_y_xxzz_xz_0,  \
                             ta1_y_xxzz_xz_1,  \
                             ta1_y_xxzz_zz_0,  \
                             ta1_y_xxzz_zz_1,  \
                             ta1_y_xyzz_yy_0,  \
                             ta1_y_xyzz_yy_1,  \
                             ta1_y_xyzz_yz_0,  \
                             ta1_y_xyzz_yz_1,  \
                             ta1_y_yzz_yy_0,   \
                             ta1_y_yzz_yy_1,   \
                             ta1_y_yzz_yz_0,   \
                             ta1_y_yzz_yz_1,   \
                             ta_xxzz_xx_1,     \
                             ta_xxzz_xz_1,     \
                             ta_xxzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzz_xx_0[i] = ta_xxzz_xx_1[i] + ta1_y_xxzz_xx_0[i] * pa_y[i] - ta1_y_xxzz_xx_1[i] * pc_y[i];

        ta1_y_xxyzz_xy_0[i] = ta1_y_xxy_xy_0[i] * fe_0 - ta1_y_xxy_xy_1[i] * fe_0 + ta1_y_xxyz_xy_0[i] * pa_z[i] - ta1_y_xxyz_xy_1[i] * pc_z[i];

        ta1_y_xxyzz_xz_0[i] = ta_xxzz_xz_1[i] + ta1_y_xxzz_xz_0[i] * pa_y[i] - ta1_y_xxzz_xz_1[i] * pc_y[i];

        ta1_y_xxyzz_yy_0[i] = ta1_y_yzz_yy_0[i] * fe_0 - ta1_y_yzz_yy_1[i] * fe_0 + ta1_y_xyzz_yy_0[i] * pa_x[i] - ta1_y_xyzz_yy_1[i] * pc_x[i];

        ta1_y_xxyzz_yz_0[i] = ta1_y_yzz_yz_0[i] * fe_0 - ta1_y_yzz_yz_1[i] * fe_0 + ta1_y_xyzz_yz_0[i] * pa_x[i] - ta1_y_xyzz_yz_1[i] * pc_x[i];

        ta1_y_xxyzz_zz_0[i] = ta_xxzz_zz_1[i] + ta1_y_xxzz_zz_0[i] * pa_y[i] - ta1_y_xxzz_zz_1[i] * pc_y[i];
    }

    // Set up 180-186 components of targeted buffer : HD

    auto ta1_y_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 180);

    auto ta1_y_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 181);

    auto ta1_y_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 182);

    auto ta1_y_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 183);

    auto ta1_y_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 184);

    auto ta1_y_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 185);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xxz_xx_0,   \
                             ta1_y_xxz_xx_1,   \
                             ta1_y_xxz_xy_0,   \
                             ta1_y_xxz_xy_1,   \
                             ta1_y_xxzz_xx_0,  \
                             ta1_y_xxzz_xx_1,  \
                             ta1_y_xxzz_xy_0,  \
                             ta1_y_xxzz_xy_1,  \
                             ta1_y_xxzzz_xx_0, \
                             ta1_y_xxzzz_xy_0, \
                             ta1_y_xxzzz_xz_0, \
                             ta1_y_xxzzz_yy_0, \
                             ta1_y_xxzzz_yz_0, \
                             ta1_y_xxzzz_zz_0, \
                             ta1_y_xzzz_xz_0,  \
                             ta1_y_xzzz_xz_1,  \
                             ta1_y_xzzz_yy_0,  \
                             ta1_y_xzzz_yy_1,  \
                             ta1_y_xzzz_yz_0,  \
                             ta1_y_xzzz_yz_1,  \
                             ta1_y_xzzz_z_0,   \
                             ta1_y_xzzz_z_1,   \
                             ta1_y_xzzz_zz_0,  \
                             ta1_y_xzzz_zz_1,  \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_yy_0,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yz_0,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzz_xx_0[i] =
            2.0 * ta1_y_xxz_xx_0[i] * fe_0 - 2.0 * ta1_y_xxz_xx_1[i] * fe_0 + ta1_y_xxzz_xx_0[i] * pa_z[i] - ta1_y_xxzz_xx_1[i] * pc_z[i];

        ta1_y_xxzzz_xy_0[i] =
            2.0 * ta1_y_xxz_xy_0[i] * fe_0 - 2.0 * ta1_y_xxz_xy_1[i] * fe_0 + ta1_y_xxzz_xy_0[i] * pa_z[i] - ta1_y_xxzz_xy_1[i] * pc_z[i];

        ta1_y_xxzzz_xz_0[i] = ta1_y_zzz_xz_0[i] * fe_0 - ta1_y_zzz_xz_1[i] * fe_0 + ta1_y_xzzz_z_0[i] * fe_0 - ta1_y_xzzz_z_1[i] * fe_0 +
                              ta1_y_xzzz_xz_0[i] * pa_x[i] - ta1_y_xzzz_xz_1[i] * pc_x[i];

        ta1_y_xxzzz_yy_0[i] = ta1_y_zzz_yy_0[i] * fe_0 - ta1_y_zzz_yy_1[i] * fe_0 + ta1_y_xzzz_yy_0[i] * pa_x[i] - ta1_y_xzzz_yy_1[i] * pc_x[i];

        ta1_y_xxzzz_yz_0[i] = ta1_y_zzz_yz_0[i] * fe_0 - ta1_y_zzz_yz_1[i] * fe_0 + ta1_y_xzzz_yz_0[i] * pa_x[i] - ta1_y_xzzz_yz_1[i] * pc_x[i];

        ta1_y_xxzzz_zz_0[i] = ta1_y_zzz_zz_0[i] * fe_0 - ta1_y_zzz_zz_1[i] * fe_0 + ta1_y_xzzz_zz_0[i] * pa_x[i] - ta1_y_xzzz_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : HD

    auto ta1_y_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 186);

    auto ta1_y_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 187);

    auto ta1_y_xyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 188);

    auto ta1_y_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 189);

    auto ta1_y_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 190);

    auto ta1_y_xyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 191);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xyyyy_xx_0, \
                             ta1_y_xyyyy_xy_0, \
                             ta1_y_xyyyy_xz_0, \
                             ta1_y_xyyyy_yy_0, \
                             ta1_y_xyyyy_yz_0, \
                             ta1_y_xyyyy_zz_0, \
                             ta1_y_yyyy_x_0,   \
                             ta1_y_yyyy_x_1,   \
                             ta1_y_yyyy_xx_0,  \
                             ta1_y_yyyy_xx_1,  \
                             ta1_y_yyyy_xy_0,  \
                             ta1_y_yyyy_xy_1,  \
                             ta1_y_yyyy_xz_0,  \
                             ta1_y_yyyy_xz_1,  \
                             ta1_y_yyyy_y_0,   \
                             ta1_y_yyyy_y_1,   \
                             ta1_y_yyyy_yy_0,  \
                             ta1_y_yyyy_yy_1,  \
                             ta1_y_yyyy_yz_0,  \
                             ta1_y_yyyy_yz_1,  \
                             ta1_y_yyyy_z_0,   \
                             ta1_y_yyyy_z_1,   \
                             ta1_y_yyyy_zz_0,  \
                             ta1_y_yyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyy_xx_0[i] =
            2.0 * ta1_y_yyyy_x_0[i] * fe_0 - 2.0 * ta1_y_yyyy_x_1[i] * fe_0 + ta1_y_yyyy_xx_0[i] * pa_x[i] - ta1_y_yyyy_xx_1[i] * pc_x[i];

        ta1_y_xyyyy_xy_0[i] = ta1_y_yyyy_y_0[i] * fe_0 - ta1_y_yyyy_y_1[i] * fe_0 + ta1_y_yyyy_xy_0[i] * pa_x[i] - ta1_y_yyyy_xy_1[i] * pc_x[i];

        ta1_y_xyyyy_xz_0[i] = ta1_y_yyyy_z_0[i] * fe_0 - ta1_y_yyyy_z_1[i] * fe_0 + ta1_y_yyyy_xz_0[i] * pa_x[i] - ta1_y_yyyy_xz_1[i] * pc_x[i];

        ta1_y_xyyyy_yy_0[i] = ta1_y_yyyy_yy_0[i] * pa_x[i] - ta1_y_yyyy_yy_1[i] * pc_x[i];

        ta1_y_xyyyy_yz_0[i] = ta1_y_yyyy_yz_0[i] * pa_x[i] - ta1_y_yyyy_yz_1[i] * pc_x[i];

        ta1_y_xyyyy_zz_0[i] = ta1_y_yyyy_zz_0[i] * pa_x[i] - ta1_y_yyyy_zz_1[i] * pc_x[i];
    }

    // Set up 192-198 components of targeted buffer : HD

    auto ta1_y_xyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 192);

    auto ta1_y_xyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 193);

    auto ta1_y_xyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 194);

    auto ta1_y_xyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 195);

    auto ta1_y_xyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 196);

    auto ta1_y_xyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 197);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_y_xyyy_xx_0,  \
                             ta1_y_xyyy_xx_1,  \
                             ta1_y_xyyy_xy_0,  \
                             ta1_y_xyyy_xy_1,  \
                             ta1_y_xyyyz_xx_0, \
                             ta1_y_xyyyz_xy_0, \
                             ta1_y_xyyyz_xz_0, \
                             ta1_y_xyyyz_yy_0, \
                             ta1_y_xyyyz_yz_0, \
                             ta1_y_xyyyz_zz_0, \
                             ta1_y_yyyz_xz_0,  \
                             ta1_y_yyyz_xz_1,  \
                             ta1_y_yyyz_yy_0,  \
                             ta1_y_yyyz_yy_1,  \
                             ta1_y_yyyz_yz_0,  \
                             ta1_y_yyyz_yz_1,  \
                             ta1_y_yyyz_z_0,   \
                             ta1_y_yyyz_z_1,   \
                             ta1_y_yyyz_zz_0,  \
                             ta1_y_yyyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyz_xx_0[i] = ta1_y_xyyy_xx_0[i] * pa_z[i] - ta1_y_xyyy_xx_1[i] * pc_z[i];

        ta1_y_xyyyz_xy_0[i] = ta1_y_xyyy_xy_0[i] * pa_z[i] - ta1_y_xyyy_xy_1[i] * pc_z[i];

        ta1_y_xyyyz_xz_0[i] = ta1_y_yyyz_z_0[i] * fe_0 - ta1_y_yyyz_z_1[i] * fe_0 + ta1_y_yyyz_xz_0[i] * pa_x[i] - ta1_y_yyyz_xz_1[i] * pc_x[i];

        ta1_y_xyyyz_yy_0[i] = ta1_y_yyyz_yy_0[i] * pa_x[i] - ta1_y_yyyz_yy_1[i] * pc_x[i];

        ta1_y_xyyyz_yz_0[i] = ta1_y_yyyz_yz_0[i] * pa_x[i] - ta1_y_yyyz_yz_1[i] * pc_x[i];

        ta1_y_xyyyz_zz_0[i] = ta1_y_yyyz_zz_0[i] * pa_x[i] - ta1_y_yyyz_zz_1[i] * pc_x[i];
    }

    // Set up 198-204 components of targeted buffer : HD

    auto ta1_y_xyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 198);

    auto ta1_y_xyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 199);

    auto ta1_y_xyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 200);

    auto ta1_y_xyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 201);

    auto ta1_y_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 202);

    auto ta1_y_xyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 203);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xyyzz_xx_0, \
                             ta1_y_xyyzz_xy_0, \
                             ta1_y_xyyzz_xz_0, \
                             ta1_y_xyyzz_yy_0, \
                             ta1_y_xyyzz_yz_0, \
                             ta1_y_xyyzz_zz_0, \
                             ta1_y_yyzz_x_0,   \
                             ta1_y_yyzz_x_1,   \
                             ta1_y_yyzz_xx_0,  \
                             ta1_y_yyzz_xx_1,  \
                             ta1_y_yyzz_xy_0,  \
                             ta1_y_yyzz_xy_1,  \
                             ta1_y_yyzz_xz_0,  \
                             ta1_y_yyzz_xz_1,  \
                             ta1_y_yyzz_y_0,   \
                             ta1_y_yyzz_y_1,   \
                             ta1_y_yyzz_yy_0,  \
                             ta1_y_yyzz_yy_1,  \
                             ta1_y_yyzz_yz_0,  \
                             ta1_y_yyzz_yz_1,  \
                             ta1_y_yyzz_z_0,   \
                             ta1_y_yyzz_z_1,   \
                             ta1_y_yyzz_zz_0,  \
                             ta1_y_yyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzz_xx_0[i] =
            2.0 * ta1_y_yyzz_x_0[i] * fe_0 - 2.0 * ta1_y_yyzz_x_1[i] * fe_0 + ta1_y_yyzz_xx_0[i] * pa_x[i] - ta1_y_yyzz_xx_1[i] * pc_x[i];

        ta1_y_xyyzz_xy_0[i] = ta1_y_yyzz_y_0[i] * fe_0 - ta1_y_yyzz_y_1[i] * fe_0 + ta1_y_yyzz_xy_0[i] * pa_x[i] - ta1_y_yyzz_xy_1[i] * pc_x[i];

        ta1_y_xyyzz_xz_0[i] = ta1_y_yyzz_z_0[i] * fe_0 - ta1_y_yyzz_z_1[i] * fe_0 + ta1_y_yyzz_xz_0[i] * pa_x[i] - ta1_y_yyzz_xz_1[i] * pc_x[i];

        ta1_y_xyyzz_yy_0[i] = ta1_y_yyzz_yy_0[i] * pa_x[i] - ta1_y_yyzz_yy_1[i] * pc_x[i];

        ta1_y_xyyzz_yz_0[i] = ta1_y_yyzz_yz_0[i] * pa_x[i] - ta1_y_yyzz_yz_1[i] * pc_x[i];

        ta1_y_xyyzz_zz_0[i] = ta1_y_yyzz_zz_0[i] * pa_x[i] - ta1_y_yyzz_zz_1[i] * pc_x[i];
    }

    // Set up 204-210 components of targeted buffer : HD

    auto ta1_y_xyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 204);

    auto ta1_y_xyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 205);

    auto ta1_y_xyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 206);

    auto ta1_y_xyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 207);

    auto ta1_y_xyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 208);

    auto ta1_y_xyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 209);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_y_xyzzz_xx_0, \
                             ta1_y_xyzzz_xy_0, \
                             ta1_y_xyzzz_xz_0, \
                             ta1_y_xyzzz_yy_0, \
                             ta1_y_xyzzz_yz_0, \
                             ta1_y_xyzzz_zz_0, \
                             ta1_y_xzzz_xx_0,  \
                             ta1_y_xzzz_xx_1,  \
                             ta1_y_xzzz_xz_0,  \
                             ta1_y_xzzz_xz_1,  \
                             ta1_y_yzzz_xy_0,  \
                             ta1_y_yzzz_xy_1,  \
                             ta1_y_yzzz_y_0,   \
                             ta1_y_yzzz_y_1,   \
                             ta1_y_yzzz_yy_0,  \
                             ta1_y_yzzz_yy_1,  \
                             ta1_y_yzzz_yz_0,  \
                             ta1_y_yzzz_yz_1,  \
                             ta1_y_yzzz_zz_0,  \
                             ta1_y_yzzz_zz_1,  \
                             ta_xzzz_xx_1,     \
                             ta_xzzz_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzz_xx_0[i] = ta_xzzz_xx_1[i] + ta1_y_xzzz_xx_0[i] * pa_y[i] - ta1_y_xzzz_xx_1[i] * pc_y[i];

        ta1_y_xyzzz_xy_0[i] = ta1_y_yzzz_y_0[i] * fe_0 - ta1_y_yzzz_y_1[i] * fe_0 + ta1_y_yzzz_xy_0[i] * pa_x[i] - ta1_y_yzzz_xy_1[i] * pc_x[i];

        ta1_y_xyzzz_xz_0[i] = ta_xzzz_xz_1[i] + ta1_y_xzzz_xz_0[i] * pa_y[i] - ta1_y_xzzz_xz_1[i] * pc_y[i];

        ta1_y_xyzzz_yy_0[i] = ta1_y_yzzz_yy_0[i] * pa_x[i] - ta1_y_yzzz_yy_1[i] * pc_x[i];

        ta1_y_xyzzz_yz_0[i] = ta1_y_yzzz_yz_0[i] * pa_x[i] - ta1_y_yzzz_yz_1[i] * pc_x[i];

        ta1_y_xyzzz_zz_0[i] = ta1_y_yzzz_zz_0[i] * pa_x[i] - ta1_y_yzzz_zz_1[i] * pc_x[i];
    }

    // Set up 210-216 components of targeted buffer : HD

    auto ta1_y_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 210);

    auto ta1_y_xzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 211);

    auto ta1_y_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 212);

    auto ta1_y_xzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 213);

    auto ta1_y_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 214);

    auto ta1_y_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 215);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_y_xzzzz_xx_0, \
                             ta1_y_xzzzz_xy_0, \
                             ta1_y_xzzzz_xz_0, \
                             ta1_y_xzzzz_yy_0, \
                             ta1_y_xzzzz_yz_0, \
                             ta1_y_xzzzz_zz_0, \
                             ta1_y_zzzz_x_0,   \
                             ta1_y_zzzz_x_1,   \
                             ta1_y_zzzz_xx_0,  \
                             ta1_y_zzzz_xx_1,  \
                             ta1_y_zzzz_xy_0,  \
                             ta1_y_zzzz_xy_1,  \
                             ta1_y_zzzz_xz_0,  \
                             ta1_y_zzzz_xz_1,  \
                             ta1_y_zzzz_y_0,   \
                             ta1_y_zzzz_y_1,   \
                             ta1_y_zzzz_yy_0,  \
                             ta1_y_zzzz_yy_1,  \
                             ta1_y_zzzz_yz_0,  \
                             ta1_y_zzzz_yz_1,  \
                             ta1_y_zzzz_z_0,   \
                             ta1_y_zzzz_z_1,   \
                             ta1_y_zzzz_zz_0,  \
                             ta1_y_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzz_xx_0[i] =
            2.0 * ta1_y_zzzz_x_0[i] * fe_0 - 2.0 * ta1_y_zzzz_x_1[i] * fe_0 + ta1_y_zzzz_xx_0[i] * pa_x[i] - ta1_y_zzzz_xx_1[i] * pc_x[i];

        ta1_y_xzzzz_xy_0[i] = ta1_y_zzzz_y_0[i] * fe_0 - ta1_y_zzzz_y_1[i] * fe_0 + ta1_y_zzzz_xy_0[i] * pa_x[i] - ta1_y_zzzz_xy_1[i] * pc_x[i];

        ta1_y_xzzzz_xz_0[i] = ta1_y_zzzz_z_0[i] * fe_0 - ta1_y_zzzz_z_1[i] * fe_0 + ta1_y_zzzz_xz_0[i] * pa_x[i] - ta1_y_zzzz_xz_1[i] * pc_x[i];

        ta1_y_xzzzz_yy_0[i] = ta1_y_zzzz_yy_0[i] * pa_x[i] - ta1_y_zzzz_yy_1[i] * pc_x[i];

        ta1_y_xzzzz_yz_0[i] = ta1_y_zzzz_yz_0[i] * pa_x[i] - ta1_y_zzzz_yz_1[i] * pc_x[i];

        ta1_y_xzzzz_zz_0[i] = ta1_y_zzzz_zz_0[i] * pa_x[i] - ta1_y_zzzz_zz_1[i] * pc_x[i];
    }

    // Set up 216-222 components of targeted buffer : HD

    auto ta1_y_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 216);

    auto ta1_y_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 217);

    auto ta1_y_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 218);

    auto ta1_y_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 219);

    auto ta1_y_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 220);

    auto ta1_y_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 221);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_y_yyy_xx_0,   \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_xz_0,   \
                             ta1_y_yyy_xz_1,   \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyy_zz_0,   \
                             ta1_y_yyy_zz_1,   \
                             ta1_y_yyyy_x_0,   \
                             ta1_y_yyyy_x_1,   \
                             ta1_y_yyyy_xx_0,  \
                             ta1_y_yyyy_xx_1,  \
                             ta1_y_yyyy_xy_0,  \
                             ta1_y_yyyy_xy_1,  \
                             ta1_y_yyyy_xz_0,  \
                             ta1_y_yyyy_xz_1,  \
                             ta1_y_yyyy_y_0,   \
                             ta1_y_yyyy_y_1,   \
                             ta1_y_yyyy_yy_0,  \
                             ta1_y_yyyy_yy_1,  \
                             ta1_y_yyyy_yz_0,  \
                             ta1_y_yyyy_yz_1,  \
                             ta1_y_yyyy_z_0,   \
                             ta1_y_yyyy_z_1,   \
                             ta1_y_yyyy_zz_0,  \
                             ta1_y_yyyy_zz_1,  \
                             ta1_y_yyyyy_xx_0, \
                             ta1_y_yyyyy_xy_0, \
                             ta1_y_yyyyy_xz_0, \
                             ta1_y_yyyyy_yy_0, \
                             ta1_y_yyyyy_yz_0, \
                             ta1_y_yyyyy_zz_0, \
                             ta_yyyy_xx_1,     \
                             ta_yyyy_xy_1,     \
                             ta_yyyy_xz_1,     \
                             ta_yyyy_yy_1,     \
                             ta_yyyy_yz_1,     \
                             ta_yyyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyy_xx_0[i] = 4.0 * ta1_y_yyy_xx_0[i] * fe_0 - 4.0 * ta1_y_yyy_xx_1[i] * fe_0 + ta_yyyy_xx_1[i] + ta1_y_yyyy_xx_0[i] * pa_y[i] -
                              ta1_y_yyyy_xx_1[i] * pc_y[i];

        ta1_y_yyyyy_xy_0[i] = 4.0 * ta1_y_yyy_xy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xy_1[i] * fe_0 + ta1_y_yyyy_x_0[i] * fe_0 - ta1_y_yyyy_x_1[i] * fe_0 +
                              ta_yyyy_xy_1[i] + ta1_y_yyyy_xy_0[i] * pa_y[i] - ta1_y_yyyy_xy_1[i] * pc_y[i];

        ta1_y_yyyyy_xz_0[i] = 4.0 * ta1_y_yyy_xz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xz_1[i] * fe_0 + ta_yyyy_xz_1[i] + ta1_y_yyyy_xz_0[i] * pa_y[i] -
                              ta1_y_yyyy_xz_1[i] * pc_y[i];

        ta1_y_yyyyy_yy_0[i] = 4.0 * ta1_y_yyy_yy_0[i] * fe_0 - 4.0 * ta1_y_yyy_yy_1[i] * fe_0 + 2.0 * ta1_y_yyyy_y_0[i] * fe_0 -
                              2.0 * ta1_y_yyyy_y_1[i] * fe_0 + ta_yyyy_yy_1[i] + ta1_y_yyyy_yy_0[i] * pa_y[i] - ta1_y_yyyy_yy_1[i] * pc_y[i];

        ta1_y_yyyyy_yz_0[i] = 4.0 * ta1_y_yyy_yz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yz_1[i] * fe_0 + ta1_y_yyyy_z_0[i] * fe_0 - ta1_y_yyyy_z_1[i] * fe_0 +
                              ta_yyyy_yz_1[i] + ta1_y_yyyy_yz_0[i] * pa_y[i] - ta1_y_yyyy_yz_1[i] * pc_y[i];

        ta1_y_yyyyy_zz_0[i] = 4.0 * ta1_y_yyy_zz_0[i] * fe_0 - 4.0 * ta1_y_yyy_zz_1[i] * fe_0 + ta_yyyy_zz_1[i] + ta1_y_yyyy_zz_0[i] * pa_y[i] -
                              ta1_y_yyyy_zz_1[i] * pc_y[i];
    }

    // Set up 222-228 components of targeted buffer : HD

    auto ta1_y_yyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 222);

    auto ta1_y_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 223);

    auto ta1_y_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 224);

    auto ta1_y_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 225);

    auto ta1_y_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 226);

    auto ta1_y_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 227);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_yyyy_x_0,   \
                             ta1_y_yyyy_x_1,   \
                             ta1_y_yyyy_xx_0,  \
                             ta1_y_yyyy_xx_1,  \
                             ta1_y_yyyy_xy_0,  \
                             ta1_y_yyyy_xy_1,  \
                             ta1_y_yyyy_xz_0,  \
                             ta1_y_yyyy_xz_1,  \
                             ta1_y_yyyy_y_0,   \
                             ta1_y_yyyy_y_1,   \
                             ta1_y_yyyy_yy_0,  \
                             ta1_y_yyyy_yy_1,  \
                             ta1_y_yyyy_yz_0,  \
                             ta1_y_yyyy_yz_1,  \
                             ta1_y_yyyy_z_0,   \
                             ta1_y_yyyy_z_1,   \
                             ta1_y_yyyy_zz_0,  \
                             ta1_y_yyyy_zz_1,  \
                             ta1_y_yyyyz_xx_0, \
                             ta1_y_yyyyz_xy_0, \
                             ta1_y_yyyyz_xz_0, \
                             ta1_y_yyyyz_yy_0, \
                             ta1_y_yyyyz_yz_0, \
                             ta1_y_yyyyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyz_xx_0[i] = ta1_y_yyyy_xx_0[i] * pa_z[i] - ta1_y_yyyy_xx_1[i] * pc_z[i];

        ta1_y_yyyyz_xy_0[i] = ta1_y_yyyy_xy_0[i] * pa_z[i] - ta1_y_yyyy_xy_1[i] * pc_z[i];

        ta1_y_yyyyz_xz_0[i] = ta1_y_yyyy_x_0[i] * fe_0 - ta1_y_yyyy_x_1[i] * fe_0 + ta1_y_yyyy_xz_0[i] * pa_z[i] - ta1_y_yyyy_xz_1[i] * pc_z[i];

        ta1_y_yyyyz_yy_0[i] = ta1_y_yyyy_yy_0[i] * pa_z[i] - ta1_y_yyyy_yy_1[i] * pc_z[i];

        ta1_y_yyyyz_yz_0[i] = ta1_y_yyyy_y_0[i] * fe_0 - ta1_y_yyyy_y_1[i] * fe_0 + ta1_y_yyyy_yz_0[i] * pa_z[i] - ta1_y_yyyy_yz_1[i] * pc_z[i];

        ta1_y_yyyyz_zz_0[i] =
            2.0 * ta1_y_yyyy_z_0[i] * fe_0 - 2.0 * ta1_y_yyyy_z_1[i] * fe_0 + ta1_y_yyyy_zz_0[i] * pa_z[i] - ta1_y_yyyy_zz_1[i] * pc_z[i];
    }

    // Set up 228-234 components of targeted buffer : HD

    auto ta1_y_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 228);

    auto ta1_y_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 229);

    auto ta1_y_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 230);

    auto ta1_y_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 231);

    auto ta1_y_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 232);

    auto ta1_y_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 233);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yyy_xx_0,   \
                             ta1_y_yyy_xx_1,   \
                             ta1_y_yyy_xy_0,   \
                             ta1_y_yyy_xy_1,   \
                             ta1_y_yyy_yy_0,   \
                             ta1_y_yyy_yy_1,   \
                             ta1_y_yyy_yz_0,   \
                             ta1_y_yyy_yz_1,   \
                             ta1_y_yyyz_xx_0,  \
                             ta1_y_yyyz_xx_1,  \
                             ta1_y_yyyz_xy_0,  \
                             ta1_y_yyyz_xy_1,  \
                             ta1_y_yyyz_y_0,   \
                             ta1_y_yyyz_y_1,   \
                             ta1_y_yyyz_yy_0,  \
                             ta1_y_yyyz_yy_1,  \
                             ta1_y_yyyz_yz_0,  \
                             ta1_y_yyyz_yz_1,  \
                             ta1_y_yyyzz_xx_0, \
                             ta1_y_yyyzz_xy_0, \
                             ta1_y_yyyzz_xz_0, \
                             ta1_y_yyyzz_yy_0, \
                             ta1_y_yyyzz_yz_0, \
                             ta1_y_yyyzz_zz_0, \
                             ta1_y_yyzz_xz_0,  \
                             ta1_y_yyzz_xz_1,  \
                             ta1_y_yyzz_zz_0,  \
                             ta1_y_yyzz_zz_1,  \
                             ta1_y_yzz_xz_0,   \
                             ta1_y_yzz_xz_1,   \
                             ta1_y_yzz_zz_0,   \
                             ta1_y_yzz_zz_1,   \
                             ta_yyzz_xz_1,     \
                             ta_yyzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzz_xx_0[i] = ta1_y_yyy_xx_0[i] * fe_0 - ta1_y_yyy_xx_1[i] * fe_0 + ta1_y_yyyz_xx_0[i] * pa_z[i] - ta1_y_yyyz_xx_1[i] * pc_z[i];

        ta1_y_yyyzz_xy_0[i] = ta1_y_yyy_xy_0[i] * fe_0 - ta1_y_yyy_xy_1[i] * fe_0 + ta1_y_yyyz_xy_0[i] * pa_z[i] - ta1_y_yyyz_xy_1[i] * pc_z[i];

        ta1_y_yyyzz_xz_0[i] = 2.0 * ta1_y_yzz_xz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xz_1[i] * fe_0 + ta_yyzz_xz_1[i] + ta1_y_yyzz_xz_0[i] * pa_y[i] -
                              ta1_y_yyzz_xz_1[i] * pc_y[i];

        ta1_y_yyyzz_yy_0[i] = ta1_y_yyy_yy_0[i] * fe_0 - ta1_y_yyy_yy_1[i] * fe_0 + ta1_y_yyyz_yy_0[i] * pa_z[i] - ta1_y_yyyz_yy_1[i] * pc_z[i];

        ta1_y_yyyzz_yz_0[i] = ta1_y_yyy_yz_0[i] * fe_0 - ta1_y_yyy_yz_1[i] * fe_0 + ta1_y_yyyz_y_0[i] * fe_0 - ta1_y_yyyz_y_1[i] * fe_0 +
                              ta1_y_yyyz_yz_0[i] * pa_z[i] - ta1_y_yyyz_yz_1[i] * pc_z[i];

        ta1_y_yyyzz_zz_0[i] = 2.0 * ta1_y_yzz_zz_0[i] * fe_0 - 2.0 * ta1_y_yzz_zz_1[i] * fe_0 + ta_yyzz_zz_1[i] + ta1_y_yyzz_zz_0[i] * pa_y[i] -
                              ta1_y_yyzz_zz_1[i] * pc_y[i];
    }

    // Set up 234-240 components of targeted buffer : HD

    auto ta1_y_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 234);

    auto ta1_y_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 235);

    auto ta1_y_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 236);

    auto ta1_y_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 237);

    auto ta1_y_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 238);

    auto ta1_y_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 239);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yyz_xx_0,   \
                             ta1_y_yyz_xx_1,   \
                             ta1_y_yyz_xy_0,   \
                             ta1_y_yyz_xy_1,   \
                             ta1_y_yyz_yy_0,   \
                             ta1_y_yyz_yy_1,   \
                             ta1_y_yyz_yz_0,   \
                             ta1_y_yyz_yz_1,   \
                             ta1_y_yyzz_xx_0,  \
                             ta1_y_yyzz_xx_1,  \
                             ta1_y_yyzz_xy_0,  \
                             ta1_y_yyzz_xy_1,  \
                             ta1_y_yyzz_y_0,   \
                             ta1_y_yyzz_y_1,   \
                             ta1_y_yyzz_yy_0,  \
                             ta1_y_yyzz_yy_1,  \
                             ta1_y_yyzz_yz_0,  \
                             ta1_y_yyzz_yz_1,  \
                             ta1_y_yyzzz_xx_0, \
                             ta1_y_yyzzz_xy_0, \
                             ta1_y_yyzzz_xz_0, \
                             ta1_y_yyzzz_yy_0, \
                             ta1_y_yyzzz_yz_0, \
                             ta1_y_yyzzz_zz_0, \
                             ta1_y_yzzz_xz_0,  \
                             ta1_y_yzzz_xz_1,  \
                             ta1_y_yzzz_zz_0,  \
                             ta1_y_yzzz_zz_1,  \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             ta_yzzz_xz_1,     \
                             ta_yzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzz_xx_0[i] =
            2.0 * ta1_y_yyz_xx_0[i] * fe_0 - 2.0 * ta1_y_yyz_xx_1[i] * fe_0 + ta1_y_yyzz_xx_0[i] * pa_z[i] - ta1_y_yyzz_xx_1[i] * pc_z[i];

        ta1_y_yyzzz_xy_0[i] =
            2.0 * ta1_y_yyz_xy_0[i] * fe_0 - 2.0 * ta1_y_yyz_xy_1[i] * fe_0 + ta1_y_yyzz_xy_0[i] * pa_z[i] - ta1_y_yyzz_xy_1[i] * pc_z[i];

        ta1_y_yyzzz_xz_0[i] =
            ta1_y_zzz_xz_0[i] * fe_0 - ta1_y_zzz_xz_1[i] * fe_0 + ta_yzzz_xz_1[i] + ta1_y_yzzz_xz_0[i] * pa_y[i] - ta1_y_yzzz_xz_1[i] * pc_y[i];

        ta1_y_yyzzz_yy_0[i] =
            2.0 * ta1_y_yyz_yy_0[i] * fe_0 - 2.0 * ta1_y_yyz_yy_1[i] * fe_0 + ta1_y_yyzz_yy_0[i] * pa_z[i] - ta1_y_yyzz_yy_1[i] * pc_z[i];

        ta1_y_yyzzz_yz_0[i] = 2.0 * ta1_y_yyz_yz_0[i] * fe_0 - 2.0 * ta1_y_yyz_yz_1[i] * fe_0 + ta1_y_yyzz_y_0[i] * fe_0 - ta1_y_yyzz_y_1[i] * fe_0 +
                              ta1_y_yyzz_yz_0[i] * pa_z[i] - ta1_y_yyzz_yz_1[i] * pc_z[i];

        ta1_y_yyzzz_zz_0[i] =
            ta1_y_zzz_zz_0[i] * fe_0 - ta1_y_zzz_zz_1[i] * fe_0 + ta_yzzz_zz_1[i] + ta1_y_yzzz_zz_0[i] * pa_y[i] - ta1_y_yzzz_zz_1[i] * pc_y[i];
    }

    // Set up 240-246 components of targeted buffer : HD

    auto ta1_y_yzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 240);

    auto ta1_y_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 241);

    auto ta1_y_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 242);

    auto ta1_y_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 243);

    auto ta1_y_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 244);

    auto ta1_y_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 245);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_y_yzz_xy_0,   \
                             ta1_y_yzz_xy_1,   \
                             ta1_y_yzz_yy_0,   \
                             ta1_y_yzz_yy_1,   \
                             ta1_y_yzzz_xy_0,  \
                             ta1_y_yzzz_xy_1,  \
                             ta1_y_yzzz_yy_0,  \
                             ta1_y_yzzz_yy_1,  \
                             ta1_y_yzzzz_xx_0, \
                             ta1_y_yzzzz_xy_0, \
                             ta1_y_yzzzz_xz_0, \
                             ta1_y_yzzzz_yy_0, \
                             ta1_y_yzzzz_yz_0, \
                             ta1_y_yzzzz_zz_0, \
                             ta1_y_zzzz_xx_0,  \
                             ta1_y_zzzz_xx_1,  \
                             ta1_y_zzzz_xz_0,  \
                             ta1_y_zzzz_xz_1,  \
                             ta1_y_zzzz_yz_0,  \
                             ta1_y_zzzz_yz_1,  \
                             ta1_y_zzzz_z_0,   \
                             ta1_y_zzzz_z_1,   \
                             ta1_y_zzzz_zz_0,  \
                             ta1_y_zzzz_zz_1,  \
                             ta_zzzz_xx_1,     \
                             ta_zzzz_xz_1,     \
                             ta_zzzz_yz_1,     \
                             ta_zzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzz_xx_0[i] = ta_zzzz_xx_1[i] + ta1_y_zzzz_xx_0[i] * pa_y[i] - ta1_y_zzzz_xx_1[i] * pc_y[i];

        ta1_y_yzzzz_xy_0[i] =
            3.0 * ta1_y_yzz_xy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xy_1[i] * fe_0 + ta1_y_yzzz_xy_0[i] * pa_z[i] - ta1_y_yzzz_xy_1[i] * pc_z[i];

        ta1_y_yzzzz_xz_0[i] = ta_zzzz_xz_1[i] + ta1_y_zzzz_xz_0[i] * pa_y[i] - ta1_y_zzzz_xz_1[i] * pc_y[i];

        ta1_y_yzzzz_yy_0[i] =
            3.0 * ta1_y_yzz_yy_0[i] * fe_0 - 3.0 * ta1_y_yzz_yy_1[i] * fe_0 + ta1_y_yzzz_yy_0[i] * pa_z[i] - ta1_y_yzzz_yy_1[i] * pc_z[i];

        ta1_y_yzzzz_yz_0[i] =
            ta1_y_zzzz_z_0[i] * fe_0 - ta1_y_zzzz_z_1[i] * fe_0 + ta_zzzz_yz_1[i] + ta1_y_zzzz_yz_0[i] * pa_y[i] - ta1_y_zzzz_yz_1[i] * pc_y[i];

        ta1_y_yzzzz_zz_0[i] = ta_zzzz_zz_1[i] + ta1_y_zzzz_zz_0[i] * pa_y[i] - ta1_y_zzzz_zz_1[i] * pc_y[i];
    }

    // Set up 246-252 components of targeted buffer : HD

    auto ta1_y_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 246);

    auto ta1_y_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 247);

    auto ta1_y_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 248);

    auto ta1_y_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 249);

    auto ta1_y_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 250);

    auto ta1_y_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 251);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_y_zzz_xx_0,   \
                             ta1_y_zzz_xx_1,   \
                             ta1_y_zzz_xy_0,   \
                             ta1_y_zzz_xy_1,   \
                             ta1_y_zzz_xz_0,   \
                             ta1_y_zzz_xz_1,   \
                             ta1_y_zzz_yy_0,   \
                             ta1_y_zzz_yy_1,   \
                             ta1_y_zzz_yz_0,   \
                             ta1_y_zzz_yz_1,   \
                             ta1_y_zzz_zz_0,   \
                             ta1_y_zzz_zz_1,   \
                             ta1_y_zzzz_x_0,   \
                             ta1_y_zzzz_x_1,   \
                             ta1_y_zzzz_xx_0,  \
                             ta1_y_zzzz_xx_1,  \
                             ta1_y_zzzz_xy_0,  \
                             ta1_y_zzzz_xy_1,  \
                             ta1_y_zzzz_xz_0,  \
                             ta1_y_zzzz_xz_1,  \
                             ta1_y_zzzz_y_0,   \
                             ta1_y_zzzz_y_1,   \
                             ta1_y_zzzz_yy_0,  \
                             ta1_y_zzzz_yy_1,  \
                             ta1_y_zzzz_yz_0,  \
                             ta1_y_zzzz_yz_1,  \
                             ta1_y_zzzz_z_0,   \
                             ta1_y_zzzz_z_1,   \
                             ta1_y_zzzz_zz_0,  \
                             ta1_y_zzzz_zz_1,  \
                             ta1_y_zzzzz_xx_0, \
                             ta1_y_zzzzz_xy_0, \
                             ta1_y_zzzzz_xz_0, \
                             ta1_y_zzzzz_yy_0, \
                             ta1_y_zzzzz_yz_0, \
                             ta1_y_zzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzz_xx_0[i] =
            4.0 * ta1_y_zzz_xx_0[i] * fe_0 - 4.0 * ta1_y_zzz_xx_1[i] * fe_0 + ta1_y_zzzz_xx_0[i] * pa_z[i] - ta1_y_zzzz_xx_1[i] * pc_z[i];

        ta1_y_zzzzz_xy_0[i] =
            4.0 * ta1_y_zzz_xy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xy_1[i] * fe_0 + ta1_y_zzzz_xy_0[i] * pa_z[i] - ta1_y_zzzz_xy_1[i] * pc_z[i];

        ta1_y_zzzzz_xz_0[i] = 4.0 * ta1_y_zzz_xz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xz_1[i] * fe_0 + ta1_y_zzzz_x_0[i] * fe_0 - ta1_y_zzzz_x_1[i] * fe_0 +
                              ta1_y_zzzz_xz_0[i] * pa_z[i] - ta1_y_zzzz_xz_1[i] * pc_z[i];

        ta1_y_zzzzz_yy_0[i] =
            4.0 * ta1_y_zzz_yy_0[i] * fe_0 - 4.0 * ta1_y_zzz_yy_1[i] * fe_0 + ta1_y_zzzz_yy_0[i] * pa_z[i] - ta1_y_zzzz_yy_1[i] * pc_z[i];

        ta1_y_zzzzz_yz_0[i] = 4.0 * ta1_y_zzz_yz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yz_1[i] * fe_0 + ta1_y_zzzz_y_0[i] * fe_0 - ta1_y_zzzz_y_1[i] * fe_0 +
                              ta1_y_zzzz_yz_0[i] * pa_z[i] - ta1_y_zzzz_yz_1[i] * pc_z[i];

        ta1_y_zzzzz_zz_0[i] = 4.0 * ta1_y_zzz_zz_0[i] * fe_0 - 4.0 * ta1_y_zzz_zz_1[i] * fe_0 + 2.0 * ta1_y_zzzz_z_0[i] * fe_0 -
                              2.0 * ta1_y_zzzz_z_1[i] * fe_0 + ta1_y_zzzz_zz_0[i] * pa_z[i] - ta1_y_zzzz_zz_1[i] * pc_z[i];
    }

    // Set up 252-258 components of targeted buffer : HD

    auto ta1_z_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 252);

    auto ta1_z_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 253);

    auto ta1_z_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 254);

    auto ta1_z_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 255);

    auto ta1_z_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 256);

    auto ta1_z_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 257);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xy_0,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxx_xz_0,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxx_yy_0,   \
                             ta1_z_xxx_yy_1,   \
                             ta1_z_xxx_yz_0,   \
                             ta1_z_xxx_yz_1,   \
                             ta1_z_xxx_zz_0,   \
                             ta1_z_xxx_zz_1,   \
                             ta1_z_xxxx_x_0,   \
                             ta1_z_xxxx_x_1,   \
                             ta1_z_xxxx_xx_0,  \
                             ta1_z_xxxx_xx_1,  \
                             ta1_z_xxxx_xy_0,  \
                             ta1_z_xxxx_xy_1,  \
                             ta1_z_xxxx_xz_0,  \
                             ta1_z_xxxx_xz_1,  \
                             ta1_z_xxxx_y_0,   \
                             ta1_z_xxxx_y_1,   \
                             ta1_z_xxxx_yy_0,  \
                             ta1_z_xxxx_yy_1,  \
                             ta1_z_xxxx_yz_0,  \
                             ta1_z_xxxx_yz_1,  \
                             ta1_z_xxxx_z_0,   \
                             ta1_z_xxxx_z_1,   \
                             ta1_z_xxxx_zz_0,  \
                             ta1_z_xxxx_zz_1,  \
                             ta1_z_xxxxx_xx_0, \
                             ta1_z_xxxxx_xy_0, \
                             ta1_z_xxxxx_xz_0, \
                             ta1_z_xxxxx_yy_0, \
                             ta1_z_xxxxx_yz_0, \
                             ta1_z_xxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxx_xx_0[i] = 4.0 * ta1_z_xxx_xx_0[i] * fe_0 - 4.0 * ta1_z_xxx_xx_1[i] * fe_0 + 2.0 * ta1_z_xxxx_x_0[i] * fe_0 -
                              2.0 * ta1_z_xxxx_x_1[i] * fe_0 + ta1_z_xxxx_xx_0[i] * pa_x[i] - ta1_z_xxxx_xx_1[i] * pc_x[i];

        ta1_z_xxxxx_xy_0[i] = 4.0 * ta1_z_xxx_xy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xy_1[i] * fe_0 + ta1_z_xxxx_y_0[i] * fe_0 - ta1_z_xxxx_y_1[i] * fe_0 +
                              ta1_z_xxxx_xy_0[i] * pa_x[i] - ta1_z_xxxx_xy_1[i] * pc_x[i];

        ta1_z_xxxxx_xz_0[i] = 4.0 * ta1_z_xxx_xz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xz_1[i] * fe_0 + ta1_z_xxxx_z_0[i] * fe_0 - ta1_z_xxxx_z_1[i] * fe_0 +
                              ta1_z_xxxx_xz_0[i] * pa_x[i] - ta1_z_xxxx_xz_1[i] * pc_x[i];

        ta1_z_xxxxx_yy_0[i] =
            4.0 * ta1_z_xxx_yy_0[i] * fe_0 - 4.0 * ta1_z_xxx_yy_1[i] * fe_0 + ta1_z_xxxx_yy_0[i] * pa_x[i] - ta1_z_xxxx_yy_1[i] * pc_x[i];

        ta1_z_xxxxx_yz_0[i] =
            4.0 * ta1_z_xxx_yz_0[i] * fe_0 - 4.0 * ta1_z_xxx_yz_1[i] * fe_0 + ta1_z_xxxx_yz_0[i] * pa_x[i] - ta1_z_xxxx_yz_1[i] * pc_x[i];

        ta1_z_xxxxx_zz_0[i] =
            4.0 * ta1_z_xxx_zz_0[i] * fe_0 - 4.0 * ta1_z_xxx_zz_1[i] * fe_0 + ta1_z_xxxx_zz_0[i] * pa_x[i] - ta1_z_xxxx_zz_1[i] * pc_x[i];
    }

    // Set up 258-264 components of targeted buffer : HD

    auto ta1_z_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 258);

    auto ta1_z_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 259);

    auto ta1_z_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 260);

    auto ta1_z_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 261);

    auto ta1_z_xxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 262);

    auto ta1_z_xxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 263);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxxx_x_0,   \
                             ta1_z_xxxx_x_1,   \
                             ta1_z_xxxx_xx_0,  \
                             ta1_z_xxxx_xx_1,  \
                             ta1_z_xxxx_xy_0,  \
                             ta1_z_xxxx_xy_1,  \
                             ta1_z_xxxx_xz_0,  \
                             ta1_z_xxxx_xz_1,  \
                             ta1_z_xxxx_zz_0,  \
                             ta1_z_xxxx_zz_1,  \
                             ta1_z_xxxxy_xx_0, \
                             ta1_z_xxxxy_xy_0, \
                             ta1_z_xxxxy_xz_0, \
                             ta1_z_xxxxy_yy_0, \
                             ta1_z_xxxxy_yz_0, \
                             ta1_z_xxxxy_zz_0, \
                             ta1_z_xxxy_yy_0,  \
                             ta1_z_xxxy_yy_1,  \
                             ta1_z_xxxy_yz_0,  \
                             ta1_z_xxxy_yz_1,  \
                             ta1_z_xxy_yy_0,   \
                             ta1_z_xxy_yy_1,   \
                             ta1_z_xxy_yz_0,   \
                             ta1_z_xxy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxy_xx_0[i] = ta1_z_xxxx_xx_0[i] * pa_y[i] - ta1_z_xxxx_xx_1[i] * pc_y[i];

        ta1_z_xxxxy_xy_0[i] = ta1_z_xxxx_x_0[i] * fe_0 - ta1_z_xxxx_x_1[i] * fe_0 + ta1_z_xxxx_xy_0[i] * pa_y[i] - ta1_z_xxxx_xy_1[i] * pc_y[i];

        ta1_z_xxxxy_xz_0[i] = ta1_z_xxxx_xz_0[i] * pa_y[i] - ta1_z_xxxx_xz_1[i] * pc_y[i];

        ta1_z_xxxxy_yy_0[i] =
            3.0 * ta1_z_xxy_yy_0[i] * fe_0 - 3.0 * ta1_z_xxy_yy_1[i] * fe_0 + ta1_z_xxxy_yy_0[i] * pa_x[i] - ta1_z_xxxy_yy_1[i] * pc_x[i];

        ta1_z_xxxxy_yz_0[i] =
            3.0 * ta1_z_xxy_yz_0[i] * fe_0 - 3.0 * ta1_z_xxy_yz_1[i] * fe_0 + ta1_z_xxxy_yz_0[i] * pa_x[i] - ta1_z_xxxy_yz_1[i] * pc_x[i];

        ta1_z_xxxxy_zz_0[i] = ta1_z_xxxx_zz_0[i] * pa_y[i] - ta1_z_xxxx_zz_1[i] * pc_y[i];
    }

    // Set up 264-270 components of targeted buffer : HD

    auto ta1_z_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 264);

    auto ta1_z_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 265);

    auto ta1_z_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 266);

    auto ta1_z_xxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 267);

    auto ta1_z_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 268);

    auto ta1_z_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 269);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxxx_x_0,   \
                             ta1_z_xxxx_x_1,   \
                             ta1_z_xxxx_xx_0,  \
                             ta1_z_xxxx_xx_1,  \
                             ta1_z_xxxx_xy_0,  \
                             ta1_z_xxxx_xy_1,  \
                             ta1_z_xxxx_xz_0,  \
                             ta1_z_xxxx_xz_1,  \
                             ta1_z_xxxx_yy_0,  \
                             ta1_z_xxxx_yy_1,  \
                             ta1_z_xxxxz_xx_0, \
                             ta1_z_xxxxz_xy_0, \
                             ta1_z_xxxxz_xz_0, \
                             ta1_z_xxxxz_yy_0, \
                             ta1_z_xxxxz_yz_0, \
                             ta1_z_xxxxz_zz_0, \
                             ta1_z_xxxz_yz_0,  \
                             ta1_z_xxxz_yz_1,  \
                             ta1_z_xxxz_zz_0,  \
                             ta1_z_xxxz_zz_1,  \
                             ta1_z_xxz_yz_0,   \
                             ta1_z_xxz_yz_1,   \
                             ta1_z_xxz_zz_0,   \
                             ta1_z_xxz_zz_1,   \
                             ta_xxxx_xx_1,     \
                             ta_xxxx_xy_1,     \
                             ta_xxxx_xz_1,     \
                             ta_xxxx_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxz_xx_0[i] = ta_xxxx_xx_1[i] + ta1_z_xxxx_xx_0[i] * pa_z[i] - ta1_z_xxxx_xx_1[i] * pc_z[i];

        ta1_z_xxxxz_xy_0[i] = ta_xxxx_xy_1[i] + ta1_z_xxxx_xy_0[i] * pa_z[i] - ta1_z_xxxx_xy_1[i] * pc_z[i];

        ta1_z_xxxxz_xz_0[i] =
            ta1_z_xxxx_x_0[i] * fe_0 - ta1_z_xxxx_x_1[i] * fe_0 + ta_xxxx_xz_1[i] + ta1_z_xxxx_xz_0[i] * pa_z[i] - ta1_z_xxxx_xz_1[i] * pc_z[i];

        ta1_z_xxxxz_yy_0[i] = ta_xxxx_yy_1[i] + ta1_z_xxxx_yy_0[i] * pa_z[i] - ta1_z_xxxx_yy_1[i] * pc_z[i];

        ta1_z_xxxxz_yz_0[i] =
            3.0 * ta1_z_xxz_yz_0[i] * fe_0 - 3.0 * ta1_z_xxz_yz_1[i] * fe_0 + ta1_z_xxxz_yz_0[i] * pa_x[i] - ta1_z_xxxz_yz_1[i] * pc_x[i];

        ta1_z_xxxxz_zz_0[i] =
            3.0 * ta1_z_xxz_zz_0[i] * fe_0 - 3.0 * ta1_z_xxz_zz_1[i] * fe_0 + ta1_z_xxxz_zz_0[i] * pa_x[i] - ta1_z_xxxz_zz_1[i] * pc_x[i];
    }

    // Set up 270-276 components of targeted buffer : HD

    auto ta1_z_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 270);

    auto ta1_z_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 271);

    auto ta1_z_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 272);

    auto ta1_z_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 273);

    auto ta1_z_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 274);

    auto ta1_z_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 275);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xz_0,   \
                             ta1_z_xxx_xz_1,   \
                             ta1_z_xxxy_xx_0,  \
                             ta1_z_xxxy_xx_1,  \
                             ta1_z_xxxy_xz_0,  \
                             ta1_z_xxxy_xz_1,  \
                             ta1_z_xxxyy_xx_0, \
                             ta1_z_xxxyy_xy_0, \
                             ta1_z_xxxyy_xz_0, \
                             ta1_z_xxxyy_yy_0, \
                             ta1_z_xxxyy_yz_0, \
                             ta1_z_xxxyy_zz_0, \
                             ta1_z_xxyy_xy_0,  \
                             ta1_z_xxyy_xy_1,  \
                             ta1_z_xxyy_y_0,   \
                             ta1_z_xxyy_y_1,   \
                             ta1_z_xxyy_yy_0,  \
                             ta1_z_xxyy_yy_1,  \
                             ta1_z_xxyy_yz_0,  \
                             ta1_z_xxyy_yz_1,  \
                             ta1_z_xxyy_zz_0,  \
                             ta1_z_xxyy_zz_1,  \
                             ta1_z_xyy_xy_0,   \
                             ta1_z_xyy_xy_1,   \
                             ta1_z_xyy_yy_0,   \
                             ta1_z_xyy_yy_1,   \
                             ta1_z_xyy_yz_0,   \
                             ta1_z_xyy_yz_1,   \
                             ta1_z_xyy_zz_0,   \
                             ta1_z_xyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyy_xx_0[i] = ta1_z_xxx_xx_0[i] * fe_0 - ta1_z_xxx_xx_1[i] * fe_0 + ta1_z_xxxy_xx_0[i] * pa_y[i] - ta1_z_xxxy_xx_1[i] * pc_y[i];

        ta1_z_xxxyy_xy_0[i] = 2.0 * ta1_z_xyy_xy_0[i] * fe_0 - 2.0 * ta1_z_xyy_xy_1[i] * fe_0 + ta1_z_xxyy_y_0[i] * fe_0 - ta1_z_xxyy_y_1[i] * fe_0 +
                              ta1_z_xxyy_xy_0[i] * pa_x[i] - ta1_z_xxyy_xy_1[i] * pc_x[i];

        ta1_z_xxxyy_xz_0[i] = ta1_z_xxx_xz_0[i] * fe_0 - ta1_z_xxx_xz_1[i] * fe_0 + ta1_z_xxxy_xz_0[i] * pa_y[i] - ta1_z_xxxy_xz_1[i] * pc_y[i];

        ta1_z_xxxyy_yy_0[i] =
            2.0 * ta1_z_xyy_yy_0[i] * fe_0 - 2.0 * ta1_z_xyy_yy_1[i] * fe_0 + ta1_z_xxyy_yy_0[i] * pa_x[i] - ta1_z_xxyy_yy_1[i] * pc_x[i];

        ta1_z_xxxyy_yz_0[i] =
            2.0 * ta1_z_xyy_yz_0[i] * fe_0 - 2.0 * ta1_z_xyy_yz_1[i] * fe_0 + ta1_z_xxyy_yz_0[i] * pa_x[i] - ta1_z_xxyy_yz_1[i] * pc_x[i];

        ta1_z_xxxyy_zz_0[i] =
            2.0 * ta1_z_xyy_zz_0[i] * fe_0 - 2.0 * ta1_z_xyy_zz_1[i] * fe_0 + ta1_z_xxyy_zz_0[i] * pa_x[i] - ta1_z_xxyy_zz_1[i] * pc_x[i];
    }

    // Set up 276-282 components of targeted buffer : HD

    auto ta1_z_xxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 276);

    auto ta1_z_xxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 277);

    auto ta1_z_xxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 278);

    auto ta1_z_xxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 279);

    auto ta1_z_xxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 280);

    auto ta1_z_xxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 281);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_xxxy_xy_0,  \
                             ta1_z_xxxy_xy_1,  \
                             ta1_z_xxxy_yy_0,  \
                             ta1_z_xxxy_yy_1,  \
                             ta1_z_xxxyz_xx_0, \
                             ta1_z_xxxyz_xy_0, \
                             ta1_z_xxxyz_xz_0, \
                             ta1_z_xxxyz_yy_0, \
                             ta1_z_xxxyz_yz_0, \
                             ta1_z_xxxyz_zz_0, \
                             ta1_z_xxxz_xx_0,  \
                             ta1_z_xxxz_xx_1,  \
                             ta1_z_xxxz_xz_0,  \
                             ta1_z_xxxz_xz_1,  \
                             ta1_z_xxxz_zz_0,  \
                             ta1_z_xxxz_zz_1,  \
                             ta1_z_xxyz_yz_0,  \
                             ta1_z_xxyz_yz_1,  \
                             ta1_z_xyz_yz_0,   \
                             ta1_z_xyz_yz_1,   \
                             ta_xxxy_xy_1,     \
                             ta_xxxy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyz_xx_0[i] = ta1_z_xxxz_xx_0[i] * pa_y[i] - ta1_z_xxxz_xx_1[i] * pc_y[i];

        ta1_z_xxxyz_xy_0[i] = ta_xxxy_xy_1[i] + ta1_z_xxxy_xy_0[i] * pa_z[i] - ta1_z_xxxy_xy_1[i] * pc_z[i];

        ta1_z_xxxyz_xz_0[i] = ta1_z_xxxz_xz_0[i] * pa_y[i] - ta1_z_xxxz_xz_1[i] * pc_y[i];

        ta1_z_xxxyz_yy_0[i] = ta_xxxy_yy_1[i] + ta1_z_xxxy_yy_0[i] * pa_z[i] - ta1_z_xxxy_yy_1[i] * pc_z[i];

        ta1_z_xxxyz_yz_0[i] =
            2.0 * ta1_z_xyz_yz_0[i] * fe_0 - 2.0 * ta1_z_xyz_yz_1[i] * fe_0 + ta1_z_xxyz_yz_0[i] * pa_x[i] - ta1_z_xxyz_yz_1[i] * pc_x[i];

        ta1_z_xxxyz_zz_0[i] = ta1_z_xxxz_zz_0[i] * pa_y[i] - ta1_z_xxxz_zz_1[i] * pc_y[i];
    }

    // Set up 282-288 components of targeted buffer : HD

    auto ta1_z_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 282);

    auto ta1_z_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 283);

    auto ta1_z_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 284);

    auto ta1_z_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 285);

    auto ta1_z_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 286);

    auto ta1_z_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 287);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxx_xx_0,   \
                             ta1_z_xxx_xx_1,   \
                             ta1_z_xxx_xy_0,   \
                             ta1_z_xxx_xy_1,   \
                             ta1_z_xxxz_xx_0,  \
                             ta1_z_xxxz_xx_1,  \
                             ta1_z_xxxz_xy_0,  \
                             ta1_z_xxxz_xy_1,  \
                             ta1_z_xxxzz_xx_0, \
                             ta1_z_xxxzz_xy_0, \
                             ta1_z_xxxzz_xz_0, \
                             ta1_z_xxxzz_yy_0, \
                             ta1_z_xxxzz_yz_0, \
                             ta1_z_xxxzz_zz_0, \
                             ta1_z_xxzz_xz_0,  \
                             ta1_z_xxzz_xz_1,  \
                             ta1_z_xxzz_yy_0,  \
                             ta1_z_xxzz_yy_1,  \
                             ta1_z_xxzz_yz_0,  \
                             ta1_z_xxzz_yz_1,  \
                             ta1_z_xxzz_z_0,   \
                             ta1_z_xxzz_z_1,   \
                             ta1_z_xxzz_zz_0,  \
                             ta1_z_xxzz_zz_1,  \
                             ta1_z_xzz_xz_0,   \
                             ta1_z_xzz_xz_1,   \
                             ta1_z_xzz_yy_0,   \
                             ta1_z_xzz_yy_1,   \
                             ta1_z_xzz_yz_0,   \
                             ta1_z_xzz_yz_1,   \
                             ta1_z_xzz_zz_0,   \
                             ta1_z_xzz_zz_1,   \
                             ta_xxxz_xx_1,     \
                             ta_xxxz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzz_xx_0[i] =
            ta1_z_xxx_xx_0[i] * fe_0 - ta1_z_xxx_xx_1[i] * fe_0 + ta_xxxz_xx_1[i] + ta1_z_xxxz_xx_0[i] * pa_z[i] - ta1_z_xxxz_xx_1[i] * pc_z[i];

        ta1_z_xxxzz_xy_0[i] =
            ta1_z_xxx_xy_0[i] * fe_0 - ta1_z_xxx_xy_1[i] * fe_0 + ta_xxxz_xy_1[i] + ta1_z_xxxz_xy_0[i] * pa_z[i] - ta1_z_xxxz_xy_1[i] * pc_z[i];

        ta1_z_xxxzz_xz_0[i] = 2.0 * ta1_z_xzz_xz_0[i] * fe_0 - 2.0 * ta1_z_xzz_xz_1[i] * fe_0 + ta1_z_xxzz_z_0[i] * fe_0 - ta1_z_xxzz_z_1[i] * fe_0 +
                              ta1_z_xxzz_xz_0[i] * pa_x[i] - ta1_z_xxzz_xz_1[i] * pc_x[i];

        ta1_z_xxxzz_yy_0[i] =
            2.0 * ta1_z_xzz_yy_0[i] * fe_0 - 2.0 * ta1_z_xzz_yy_1[i] * fe_0 + ta1_z_xxzz_yy_0[i] * pa_x[i] - ta1_z_xxzz_yy_1[i] * pc_x[i];

        ta1_z_xxxzz_yz_0[i] =
            2.0 * ta1_z_xzz_yz_0[i] * fe_0 - 2.0 * ta1_z_xzz_yz_1[i] * fe_0 + ta1_z_xxzz_yz_0[i] * pa_x[i] - ta1_z_xxzz_yz_1[i] * pc_x[i];

        ta1_z_xxxzz_zz_0[i] =
            2.0 * ta1_z_xzz_zz_0[i] * fe_0 - 2.0 * ta1_z_xzz_zz_1[i] * fe_0 + ta1_z_xxzz_zz_0[i] * pa_x[i] - ta1_z_xxzz_zz_1[i] * pc_x[i];
    }

    // Set up 288-294 components of targeted buffer : HD

    auto ta1_z_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 288);

    auto ta1_z_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 289);

    auto ta1_z_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 290);

    auto ta1_z_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 291);

    auto ta1_z_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 292);

    auto ta1_z_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 293);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxy_xx_0,   \
                             ta1_z_xxy_xx_1,   \
                             ta1_z_xxy_xz_0,   \
                             ta1_z_xxy_xz_1,   \
                             ta1_z_xxyy_xx_0,  \
                             ta1_z_xxyy_xx_1,  \
                             ta1_z_xxyy_xz_0,  \
                             ta1_z_xxyy_xz_1,  \
                             ta1_z_xxyyy_xx_0, \
                             ta1_z_xxyyy_xy_0, \
                             ta1_z_xxyyy_xz_0, \
                             ta1_z_xxyyy_yy_0, \
                             ta1_z_xxyyy_yz_0, \
                             ta1_z_xxyyy_zz_0, \
                             ta1_z_xyyy_xy_0,  \
                             ta1_z_xyyy_xy_1,  \
                             ta1_z_xyyy_y_0,   \
                             ta1_z_xyyy_y_1,   \
                             ta1_z_xyyy_yy_0,  \
                             ta1_z_xyyy_yy_1,  \
                             ta1_z_xyyy_yz_0,  \
                             ta1_z_xyyy_yz_1,  \
                             ta1_z_xyyy_zz_0,  \
                             ta1_z_xyyy_zz_1,  \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yz_0,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_zz_0,   \
                             ta1_z_yyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyy_xx_0[i] =
            2.0 * ta1_z_xxy_xx_0[i] * fe_0 - 2.0 * ta1_z_xxy_xx_1[i] * fe_0 + ta1_z_xxyy_xx_0[i] * pa_y[i] - ta1_z_xxyy_xx_1[i] * pc_y[i];

        ta1_z_xxyyy_xy_0[i] = ta1_z_yyy_xy_0[i] * fe_0 - ta1_z_yyy_xy_1[i] * fe_0 + ta1_z_xyyy_y_0[i] * fe_0 - ta1_z_xyyy_y_1[i] * fe_0 +
                              ta1_z_xyyy_xy_0[i] * pa_x[i] - ta1_z_xyyy_xy_1[i] * pc_x[i];

        ta1_z_xxyyy_xz_0[i] =
            2.0 * ta1_z_xxy_xz_0[i] * fe_0 - 2.0 * ta1_z_xxy_xz_1[i] * fe_0 + ta1_z_xxyy_xz_0[i] * pa_y[i] - ta1_z_xxyy_xz_1[i] * pc_y[i];

        ta1_z_xxyyy_yy_0[i] = ta1_z_yyy_yy_0[i] * fe_0 - ta1_z_yyy_yy_1[i] * fe_0 + ta1_z_xyyy_yy_0[i] * pa_x[i] - ta1_z_xyyy_yy_1[i] * pc_x[i];

        ta1_z_xxyyy_yz_0[i] = ta1_z_yyy_yz_0[i] * fe_0 - ta1_z_yyy_yz_1[i] * fe_0 + ta1_z_xyyy_yz_0[i] * pa_x[i] - ta1_z_xyyy_yz_1[i] * pc_x[i];

        ta1_z_xxyyy_zz_0[i] = ta1_z_yyy_zz_0[i] * fe_0 - ta1_z_yyy_zz_1[i] * fe_0 + ta1_z_xyyy_zz_0[i] * pa_x[i] - ta1_z_xyyy_zz_1[i] * pc_x[i];
    }

    // Set up 294-300 components of targeted buffer : HD

    auto ta1_z_xxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 294);

    auto ta1_z_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 295);

    auto ta1_z_xxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 296);

    auto ta1_z_xxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 297);

    auto ta1_z_xxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 298);

    auto ta1_z_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 299);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pa_z,             \
                             pc_x,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_xxyy_xx_0,  \
                             ta1_z_xxyy_xx_1,  \
                             ta1_z_xxyy_xy_0,  \
                             ta1_z_xxyy_xy_1,  \
                             ta1_z_xxyy_yy_0,  \
                             ta1_z_xxyy_yy_1,  \
                             ta1_z_xxyyz_xx_0, \
                             ta1_z_xxyyz_xy_0, \
                             ta1_z_xxyyz_xz_0, \
                             ta1_z_xxyyz_yy_0, \
                             ta1_z_xxyyz_yz_0, \
                             ta1_z_xxyyz_zz_0, \
                             ta1_z_xxyz_xz_0,  \
                             ta1_z_xxyz_xz_1,  \
                             ta1_z_xxz_xz_0,   \
                             ta1_z_xxz_xz_1,   \
                             ta1_z_xyyz_yz_0,  \
                             ta1_z_xyyz_yz_1,  \
                             ta1_z_xyyz_zz_0,  \
                             ta1_z_xyyz_zz_1,  \
                             ta1_z_yyz_yz_0,   \
                             ta1_z_yyz_yz_1,   \
                             ta1_z_yyz_zz_0,   \
                             ta1_z_yyz_zz_1,   \
                             ta_xxyy_xx_1,     \
                             ta_xxyy_xy_1,     \
                             ta_xxyy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyz_xx_0[i] = ta_xxyy_xx_1[i] + ta1_z_xxyy_xx_0[i] * pa_z[i] - ta1_z_xxyy_xx_1[i] * pc_z[i];

        ta1_z_xxyyz_xy_0[i] = ta_xxyy_xy_1[i] + ta1_z_xxyy_xy_0[i] * pa_z[i] - ta1_z_xxyy_xy_1[i] * pc_z[i];

        ta1_z_xxyyz_xz_0[i] = ta1_z_xxz_xz_0[i] * fe_0 - ta1_z_xxz_xz_1[i] * fe_0 + ta1_z_xxyz_xz_0[i] * pa_y[i] - ta1_z_xxyz_xz_1[i] * pc_y[i];

        ta1_z_xxyyz_yy_0[i] = ta_xxyy_yy_1[i] + ta1_z_xxyy_yy_0[i] * pa_z[i] - ta1_z_xxyy_yy_1[i] * pc_z[i];

        ta1_z_xxyyz_yz_0[i] = ta1_z_yyz_yz_0[i] * fe_0 - ta1_z_yyz_yz_1[i] * fe_0 + ta1_z_xyyz_yz_0[i] * pa_x[i] - ta1_z_xyyz_yz_1[i] * pc_x[i];

        ta1_z_xxyyz_zz_0[i] = ta1_z_yyz_zz_0[i] * fe_0 - ta1_z_yyz_zz_1[i] * fe_0 + ta1_z_xyyz_zz_0[i] * pa_x[i] - ta1_z_xyyz_zz_1[i] * pc_x[i];
    }

    // Set up 300-306 components of targeted buffer : HD

    auto ta1_z_xxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 300);

    auto ta1_z_xxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 301);

    auto ta1_z_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 302);

    auto ta1_z_xxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 303);

    auto ta1_z_xxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 304);

    auto ta1_z_xxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 305);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xxyzz_xx_0, \
                             ta1_z_xxyzz_xy_0, \
                             ta1_z_xxyzz_xz_0, \
                             ta1_z_xxyzz_yy_0, \
                             ta1_z_xxyzz_yz_0, \
                             ta1_z_xxyzz_zz_0, \
                             ta1_z_xxzz_x_0,   \
                             ta1_z_xxzz_x_1,   \
                             ta1_z_xxzz_xx_0,  \
                             ta1_z_xxzz_xx_1,  \
                             ta1_z_xxzz_xy_0,  \
                             ta1_z_xxzz_xy_1,  \
                             ta1_z_xxzz_xz_0,  \
                             ta1_z_xxzz_xz_1,  \
                             ta1_z_xxzz_zz_0,  \
                             ta1_z_xxzz_zz_1,  \
                             ta1_z_xyzz_yy_0,  \
                             ta1_z_xyzz_yy_1,  \
                             ta1_z_xyzz_yz_0,  \
                             ta1_z_xyzz_yz_1,  \
                             ta1_z_yzz_yy_0,   \
                             ta1_z_yzz_yy_1,   \
                             ta1_z_yzz_yz_0,   \
                             ta1_z_yzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzz_xx_0[i] = ta1_z_xxzz_xx_0[i] * pa_y[i] - ta1_z_xxzz_xx_1[i] * pc_y[i];

        ta1_z_xxyzz_xy_0[i] = ta1_z_xxzz_x_0[i] * fe_0 - ta1_z_xxzz_x_1[i] * fe_0 + ta1_z_xxzz_xy_0[i] * pa_y[i] - ta1_z_xxzz_xy_1[i] * pc_y[i];

        ta1_z_xxyzz_xz_0[i] = ta1_z_xxzz_xz_0[i] * pa_y[i] - ta1_z_xxzz_xz_1[i] * pc_y[i];

        ta1_z_xxyzz_yy_0[i] = ta1_z_yzz_yy_0[i] * fe_0 - ta1_z_yzz_yy_1[i] * fe_0 + ta1_z_xyzz_yy_0[i] * pa_x[i] - ta1_z_xyzz_yy_1[i] * pc_x[i];

        ta1_z_xxyzz_yz_0[i] = ta1_z_yzz_yz_0[i] * fe_0 - ta1_z_yzz_yz_1[i] * fe_0 + ta1_z_xyzz_yz_0[i] * pa_x[i] - ta1_z_xyzz_yz_1[i] * pc_x[i];

        ta1_z_xxyzz_zz_0[i] = ta1_z_xxzz_zz_0[i] * pa_y[i] - ta1_z_xxzz_zz_1[i] * pc_y[i];
    }

    // Set up 306-312 components of targeted buffer : HD

    auto ta1_z_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 306);

    auto ta1_z_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 307);

    auto ta1_z_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 308);

    auto ta1_z_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 309);

    auto ta1_z_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 310);

    auto ta1_z_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 311);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xxz_xx_0,   \
                             ta1_z_xxz_xx_1,   \
                             ta1_z_xxz_xy_0,   \
                             ta1_z_xxz_xy_1,   \
                             ta1_z_xxzz_xx_0,  \
                             ta1_z_xxzz_xx_1,  \
                             ta1_z_xxzz_xy_0,  \
                             ta1_z_xxzz_xy_1,  \
                             ta1_z_xxzzz_xx_0, \
                             ta1_z_xxzzz_xy_0, \
                             ta1_z_xxzzz_xz_0, \
                             ta1_z_xxzzz_yy_0, \
                             ta1_z_xxzzz_yz_0, \
                             ta1_z_xxzzz_zz_0, \
                             ta1_z_xzzz_xz_0,  \
                             ta1_z_xzzz_xz_1,  \
                             ta1_z_xzzz_yy_0,  \
                             ta1_z_xzzz_yy_1,  \
                             ta1_z_xzzz_yz_0,  \
                             ta1_z_xzzz_yz_1,  \
                             ta1_z_xzzz_z_0,   \
                             ta1_z_xzzz_z_1,   \
                             ta1_z_xzzz_zz_0,  \
                             ta1_z_xzzz_zz_1,  \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yy_0,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta_xxzz_xx_1,     \
                             ta_xxzz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzz_xx_0[i] = 2.0 * ta1_z_xxz_xx_0[i] * fe_0 - 2.0 * ta1_z_xxz_xx_1[i] * fe_0 + ta_xxzz_xx_1[i] + ta1_z_xxzz_xx_0[i] * pa_z[i] -
                              ta1_z_xxzz_xx_1[i] * pc_z[i];

        ta1_z_xxzzz_xy_0[i] = 2.0 * ta1_z_xxz_xy_0[i] * fe_0 - 2.0 * ta1_z_xxz_xy_1[i] * fe_0 + ta_xxzz_xy_1[i] + ta1_z_xxzz_xy_0[i] * pa_z[i] -
                              ta1_z_xxzz_xy_1[i] * pc_z[i];

        ta1_z_xxzzz_xz_0[i] = ta1_z_zzz_xz_0[i] * fe_0 - ta1_z_zzz_xz_1[i] * fe_0 + ta1_z_xzzz_z_0[i] * fe_0 - ta1_z_xzzz_z_1[i] * fe_0 +
                              ta1_z_xzzz_xz_0[i] * pa_x[i] - ta1_z_xzzz_xz_1[i] * pc_x[i];

        ta1_z_xxzzz_yy_0[i] = ta1_z_zzz_yy_0[i] * fe_0 - ta1_z_zzz_yy_1[i] * fe_0 + ta1_z_xzzz_yy_0[i] * pa_x[i] - ta1_z_xzzz_yy_1[i] * pc_x[i];

        ta1_z_xxzzz_yz_0[i] = ta1_z_zzz_yz_0[i] * fe_0 - ta1_z_zzz_yz_1[i] * fe_0 + ta1_z_xzzz_yz_0[i] * pa_x[i] - ta1_z_xzzz_yz_1[i] * pc_x[i];

        ta1_z_xxzzz_zz_0[i] = ta1_z_zzz_zz_0[i] * fe_0 - ta1_z_zzz_zz_1[i] * fe_0 + ta1_z_xzzz_zz_0[i] * pa_x[i] - ta1_z_xzzz_zz_1[i] * pc_x[i];
    }

    // Set up 312-318 components of targeted buffer : HD

    auto ta1_z_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 312);

    auto ta1_z_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 313);

    auto ta1_z_xyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 314);

    auto ta1_z_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 315);

    auto ta1_z_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 316);

    auto ta1_z_xyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 317);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xyyyy_xx_0, \
                             ta1_z_xyyyy_xy_0, \
                             ta1_z_xyyyy_xz_0, \
                             ta1_z_xyyyy_yy_0, \
                             ta1_z_xyyyy_yz_0, \
                             ta1_z_xyyyy_zz_0, \
                             ta1_z_yyyy_x_0,   \
                             ta1_z_yyyy_x_1,   \
                             ta1_z_yyyy_xx_0,  \
                             ta1_z_yyyy_xx_1,  \
                             ta1_z_yyyy_xy_0,  \
                             ta1_z_yyyy_xy_1,  \
                             ta1_z_yyyy_xz_0,  \
                             ta1_z_yyyy_xz_1,  \
                             ta1_z_yyyy_y_0,   \
                             ta1_z_yyyy_y_1,   \
                             ta1_z_yyyy_yy_0,  \
                             ta1_z_yyyy_yy_1,  \
                             ta1_z_yyyy_yz_0,  \
                             ta1_z_yyyy_yz_1,  \
                             ta1_z_yyyy_z_0,   \
                             ta1_z_yyyy_z_1,   \
                             ta1_z_yyyy_zz_0,  \
                             ta1_z_yyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyy_xx_0[i] =
            2.0 * ta1_z_yyyy_x_0[i] * fe_0 - 2.0 * ta1_z_yyyy_x_1[i] * fe_0 + ta1_z_yyyy_xx_0[i] * pa_x[i] - ta1_z_yyyy_xx_1[i] * pc_x[i];

        ta1_z_xyyyy_xy_0[i] = ta1_z_yyyy_y_0[i] * fe_0 - ta1_z_yyyy_y_1[i] * fe_0 + ta1_z_yyyy_xy_0[i] * pa_x[i] - ta1_z_yyyy_xy_1[i] * pc_x[i];

        ta1_z_xyyyy_xz_0[i] = ta1_z_yyyy_z_0[i] * fe_0 - ta1_z_yyyy_z_1[i] * fe_0 + ta1_z_yyyy_xz_0[i] * pa_x[i] - ta1_z_yyyy_xz_1[i] * pc_x[i];

        ta1_z_xyyyy_yy_0[i] = ta1_z_yyyy_yy_0[i] * pa_x[i] - ta1_z_yyyy_yy_1[i] * pc_x[i];

        ta1_z_xyyyy_yz_0[i] = ta1_z_yyyy_yz_0[i] * pa_x[i] - ta1_z_yyyy_yz_1[i] * pc_x[i];

        ta1_z_xyyyy_zz_0[i] = ta1_z_yyyy_zz_0[i] * pa_x[i] - ta1_z_yyyy_zz_1[i] * pc_x[i];
    }

    // Set up 318-324 components of targeted buffer : HD

    auto ta1_z_xyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 318);

    auto ta1_z_xyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 319);

    auto ta1_z_xyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 320);

    auto ta1_z_xyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 321);

    auto ta1_z_xyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 322);

    auto ta1_z_xyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 323);

#pragma omp simd aligned(pa_x,                 \
                             pa_z,             \
                             pc_x,             \
                             pc_z,             \
                             ta1_z_xyyy_xx_0,  \
                             ta1_z_xyyy_xx_1,  \
                             ta1_z_xyyy_xy_0,  \
                             ta1_z_xyyy_xy_1,  \
                             ta1_z_xyyyz_xx_0, \
                             ta1_z_xyyyz_xy_0, \
                             ta1_z_xyyyz_xz_0, \
                             ta1_z_xyyyz_yy_0, \
                             ta1_z_xyyyz_yz_0, \
                             ta1_z_xyyyz_zz_0, \
                             ta1_z_yyyz_xz_0,  \
                             ta1_z_yyyz_xz_1,  \
                             ta1_z_yyyz_yy_0,  \
                             ta1_z_yyyz_yy_1,  \
                             ta1_z_yyyz_yz_0,  \
                             ta1_z_yyyz_yz_1,  \
                             ta1_z_yyyz_z_0,   \
                             ta1_z_yyyz_z_1,   \
                             ta1_z_yyyz_zz_0,  \
                             ta1_z_yyyz_zz_1,  \
                             ta_xyyy_xx_1,     \
                             ta_xyyy_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyz_xx_0[i] = ta_xyyy_xx_1[i] + ta1_z_xyyy_xx_0[i] * pa_z[i] - ta1_z_xyyy_xx_1[i] * pc_z[i];

        ta1_z_xyyyz_xy_0[i] = ta_xyyy_xy_1[i] + ta1_z_xyyy_xy_0[i] * pa_z[i] - ta1_z_xyyy_xy_1[i] * pc_z[i];

        ta1_z_xyyyz_xz_0[i] = ta1_z_yyyz_z_0[i] * fe_0 - ta1_z_yyyz_z_1[i] * fe_0 + ta1_z_yyyz_xz_0[i] * pa_x[i] - ta1_z_yyyz_xz_1[i] * pc_x[i];

        ta1_z_xyyyz_yy_0[i] = ta1_z_yyyz_yy_0[i] * pa_x[i] - ta1_z_yyyz_yy_1[i] * pc_x[i];

        ta1_z_xyyyz_yz_0[i] = ta1_z_yyyz_yz_0[i] * pa_x[i] - ta1_z_yyyz_yz_1[i] * pc_x[i];

        ta1_z_xyyyz_zz_0[i] = ta1_z_yyyz_zz_0[i] * pa_x[i] - ta1_z_yyyz_zz_1[i] * pc_x[i];
    }

    // Set up 324-330 components of targeted buffer : HD

    auto ta1_z_xyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 324);

    auto ta1_z_xyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 325);

    auto ta1_z_xyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 326);

    auto ta1_z_xyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 327);

    auto ta1_z_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 328);

    auto ta1_z_xyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 329);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xyyzz_xx_0, \
                             ta1_z_xyyzz_xy_0, \
                             ta1_z_xyyzz_xz_0, \
                             ta1_z_xyyzz_yy_0, \
                             ta1_z_xyyzz_yz_0, \
                             ta1_z_xyyzz_zz_0, \
                             ta1_z_yyzz_x_0,   \
                             ta1_z_yyzz_x_1,   \
                             ta1_z_yyzz_xx_0,  \
                             ta1_z_yyzz_xx_1,  \
                             ta1_z_yyzz_xy_0,  \
                             ta1_z_yyzz_xy_1,  \
                             ta1_z_yyzz_xz_0,  \
                             ta1_z_yyzz_xz_1,  \
                             ta1_z_yyzz_y_0,   \
                             ta1_z_yyzz_y_1,   \
                             ta1_z_yyzz_yy_0,  \
                             ta1_z_yyzz_yy_1,  \
                             ta1_z_yyzz_yz_0,  \
                             ta1_z_yyzz_yz_1,  \
                             ta1_z_yyzz_z_0,   \
                             ta1_z_yyzz_z_1,   \
                             ta1_z_yyzz_zz_0,  \
                             ta1_z_yyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzz_xx_0[i] =
            2.0 * ta1_z_yyzz_x_0[i] * fe_0 - 2.0 * ta1_z_yyzz_x_1[i] * fe_0 + ta1_z_yyzz_xx_0[i] * pa_x[i] - ta1_z_yyzz_xx_1[i] * pc_x[i];

        ta1_z_xyyzz_xy_0[i] = ta1_z_yyzz_y_0[i] * fe_0 - ta1_z_yyzz_y_1[i] * fe_0 + ta1_z_yyzz_xy_0[i] * pa_x[i] - ta1_z_yyzz_xy_1[i] * pc_x[i];

        ta1_z_xyyzz_xz_0[i] = ta1_z_yyzz_z_0[i] * fe_0 - ta1_z_yyzz_z_1[i] * fe_0 + ta1_z_yyzz_xz_0[i] * pa_x[i] - ta1_z_yyzz_xz_1[i] * pc_x[i];

        ta1_z_xyyzz_yy_0[i] = ta1_z_yyzz_yy_0[i] * pa_x[i] - ta1_z_yyzz_yy_1[i] * pc_x[i];

        ta1_z_xyyzz_yz_0[i] = ta1_z_yyzz_yz_0[i] * pa_x[i] - ta1_z_yyzz_yz_1[i] * pc_x[i];

        ta1_z_xyyzz_zz_0[i] = ta1_z_yyzz_zz_0[i] * pa_x[i] - ta1_z_yyzz_zz_1[i] * pc_x[i];
    }

    // Set up 330-336 components of targeted buffer : HD

    auto ta1_z_xyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 330);

    auto ta1_z_xyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 331);

    auto ta1_z_xyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 332);

    auto ta1_z_xyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 333);

    auto ta1_z_xyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 334);

    auto ta1_z_xyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 335);

#pragma omp simd aligned(pa_x,                 \
                             pa_y,             \
                             pc_x,             \
                             pc_y,             \
                             ta1_z_xyzzz_xx_0, \
                             ta1_z_xyzzz_xy_0, \
                             ta1_z_xyzzz_xz_0, \
                             ta1_z_xyzzz_yy_0, \
                             ta1_z_xyzzz_yz_0, \
                             ta1_z_xyzzz_zz_0, \
                             ta1_z_xzzz_xx_0,  \
                             ta1_z_xzzz_xx_1,  \
                             ta1_z_xzzz_xz_0,  \
                             ta1_z_xzzz_xz_1,  \
                             ta1_z_yzzz_xy_0,  \
                             ta1_z_yzzz_xy_1,  \
                             ta1_z_yzzz_y_0,   \
                             ta1_z_yzzz_y_1,   \
                             ta1_z_yzzz_yy_0,  \
                             ta1_z_yzzz_yy_1,  \
                             ta1_z_yzzz_yz_0,  \
                             ta1_z_yzzz_yz_1,  \
                             ta1_z_yzzz_zz_0,  \
                             ta1_z_yzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzz_xx_0[i] = ta1_z_xzzz_xx_0[i] * pa_y[i] - ta1_z_xzzz_xx_1[i] * pc_y[i];

        ta1_z_xyzzz_xy_0[i] = ta1_z_yzzz_y_0[i] * fe_0 - ta1_z_yzzz_y_1[i] * fe_0 + ta1_z_yzzz_xy_0[i] * pa_x[i] - ta1_z_yzzz_xy_1[i] * pc_x[i];

        ta1_z_xyzzz_xz_0[i] = ta1_z_xzzz_xz_0[i] * pa_y[i] - ta1_z_xzzz_xz_1[i] * pc_y[i];

        ta1_z_xyzzz_yy_0[i] = ta1_z_yzzz_yy_0[i] * pa_x[i] - ta1_z_yzzz_yy_1[i] * pc_x[i];

        ta1_z_xyzzz_yz_0[i] = ta1_z_yzzz_yz_0[i] * pa_x[i] - ta1_z_yzzz_yz_1[i] * pc_x[i];

        ta1_z_xyzzz_zz_0[i] = ta1_z_yzzz_zz_0[i] * pa_x[i] - ta1_z_yzzz_zz_1[i] * pc_x[i];
    }

    // Set up 336-342 components of targeted buffer : HD

    auto ta1_z_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 336);

    auto ta1_z_xzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 337);

    auto ta1_z_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 338);

    auto ta1_z_xzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 339);

    auto ta1_z_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 340);

    auto ta1_z_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 341);

#pragma omp simd aligned(pa_x,                 \
                             pc_x,             \
                             ta1_z_xzzzz_xx_0, \
                             ta1_z_xzzzz_xy_0, \
                             ta1_z_xzzzz_xz_0, \
                             ta1_z_xzzzz_yy_0, \
                             ta1_z_xzzzz_yz_0, \
                             ta1_z_xzzzz_zz_0, \
                             ta1_z_zzzz_x_0,   \
                             ta1_z_zzzz_x_1,   \
                             ta1_z_zzzz_xx_0,  \
                             ta1_z_zzzz_xx_1,  \
                             ta1_z_zzzz_xy_0,  \
                             ta1_z_zzzz_xy_1,  \
                             ta1_z_zzzz_xz_0,  \
                             ta1_z_zzzz_xz_1,  \
                             ta1_z_zzzz_y_0,   \
                             ta1_z_zzzz_y_1,   \
                             ta1_z_zzzz_yy_0,  \
                             ta1_z_zzzz_yy_1,  \
                             ta1_z_zzzz_yz_0,  \
                             ta1_z_zzzz_yz_1,  \
                             ta1_z_zzzz_z_0,   \
                             ta1_z_zzzz_z_1,   \
                             ta1_z_zzzz_zz_0,  \
                             ta1_z_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzz_xx_0[i] =
            2.0 * ta1_z_zzzz_x_0[i] * fe_0 - 2.0 * ta1_z_zzzz_x_1[i] * fe_0 + ta1_z_zzzz_xx_0[i] * pa_x[i] - ta1_z_zzzz_xx_1[i] * pc_x[i];

        ta1_z_xzzzz_xy_0[i] = ta1_z_zzzz_y_0[i] * fe_0 - ta1_z_zzzz_y_1[i] * fe_0 + ta1_z_zzzz_xy_0[i] * pa_x[i] - ta1_z_zzzz_xy_1[i] * pc_x[i];

        ta1_z_xzzzz_xz_0[i] = ta1_z_zzzz_z_0[i] * fe_0 - ta1_z_zzzz_z_1[i] * fe_0 + ta1_z_zzzz_xz_0[i] * pa_x[i] - ta1_z_zzzz_xz_1[i] * pc_x[i];

        ta1_z_xzzzz_yy_0[i] = ta1_z_zzzz_yy_0[i] * pa_x[i] - ta1_z_zzzz_yy_1[i] * pc_x[i];

        ta1_z_xzzzz_yz_0[i] = ta1_z_zzzz_yz_0[i] * pa_x[i] - ta1_z_zzzz_yz_1[i] * pc_x[i];

        ta1_z_xzzzz_zz_0[i] = ta1_z_zzzz_zz_0[i] * pa_x[i] - ta1_z_zzzz_zz_1[i] * pc_x[i];
    }

    // Set up 342-348 components of targeted buffer : HD

    auto ta1_z_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 342);

    auto ta1_z_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 343);

    auto ta1_z_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 344);

    auto ta1_z_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 345);

    auto ta1_z_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 346);

    auto ta1_z_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 347);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yyy_xx_0,   \
                             ta1_z_yyy_xx_1,   \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_xz_0,   \
                             ta1_z_yyy_xz_1,   \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyy_yz_0,   \
                             ta1_z_yyy_yz_1,   \
                             ta1_z_yyy_zz_0,   \
                             ta1_z_yyy_zz_1,   \
                             ta1_z_yyyy_x_0,   \
                             ta1_z_yyyy_x_1,   \
                             ta1_z_yyyy_xx_0,  \
                             ta1_z_yyyy_xx_1,  \
                             ta1_z_yyyy_xy_0,  \
                             ta1_z_yyyy_xy_1,  \
                             ta1_z_yyyy_xz_0,  \
                             ta1_z_yyyy_xz_1,  \
                             ta1_z_yyyy_y_0,   \
                             ta1_z_yyyy_y_1,   \
                             ta1_z_yyyy_yy_0,  \
                             ta1_z_yyyy_yy_1,  \
                             ta1_z_yyyy_yz_0,  \
                             ta1_z_yyyy_yz_1,  \
                             ta1_z_yyyy_z_0,   \
                             ta1_z_yyyy_z_1,   \
                             ta1_z_yyyy_zz_0,  \
                             ta1_z_yyyy_zz_1,  \
                             ta1_z_yyyyy_xx_0, \
                             ta1_z_yyyyy_xy_0, \
                             ta1_z_yyyyy_xz_0, \
                             ta1_z_yyyyy_yy_0, \
                             ta1_z_yyyyy_yz_0, \
                             ta1_z_yyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyy_xx_0[i] =
            4.0 * ta1_z_yyy_xx_0[i] * fe_0 - 4.0 * ta1_z_yyy_xx_1[i] * fe_0 + ta1_z_yyyy_xx_0[i] * pa_y[i] - ta1_z_yyyy_xx_1[i] * pc_y[i];

        ta1_z_yyyyy_xy_0[i] = 4.0 * ta1_z_yyy_xy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xy_1[i] * fe_0 + ta1_z_yyyy_x_0[i] * fe_0 - ta1_z_yyyy_x_1[i] * fe_0 +
                              ta1_z_yyyy_xy_0[i] * pa_y[i] - ta1_z_yyyy_xy_1[i] * pc_y[i];

        ta1_z_yyyyy_xz_0[i] =
            4.0 * ta1_z_yyy_xz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xz_1[i] * fe_0 + ta1_z_yyyy_xz_0[i] * pa_y[i] - ta1_z_yyyy_xz_1[i] * pc_y[i];

        ta1_z_yyyyy_yy_0[i] = 4.0 * ta1_z_yyy_yy_0[i] * fe_0 - 4.0 * ta1_z_yyy_yy_1[i] * fe_0 + 2.0 * ta1_z_yyyy_y_0[i] * fe_0 -
                              2.0 * ta1_z_yyyy_y_1[i] * fe_0 + ta1_z_yyyy_yy_0[i] * pa_y[i] - ta1_z_yyyy_yy_1[i] * pc_y[i];

        ta1_z_yyyyy_yz_0[i] = 4.0 * ta1_z_yyy_yz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yz_1[i] * fe_0 + ta1_z_yyyy_z_0[i] * fe_0 - ta1_z_yyyy_z_1[i] * fe_0 +
                              ta1_z_yyyy_yz_0[i] * pa_y[i] - ta1_z_yyyy_yz_1[i] * pc_y[i];

        ta1_z_yyyyy_zz_0[i] =
            4.0 * ta1_z_yyy_zz_0[i] * fe_0 - 4.0 * ta1_z_yyy_zz_1[i] * fe_0 + ta1_z_yyyy_zz_0[i] * pa_y[i] - ta1_z_yyyy_zz_1[i] * pc_y[i];
    }

    // Set up 348-354 components of targeted buffer : HD

    auto ta1_z_yyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 348);

    auto ta1_z_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 349);

    auto ta1_z_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 350);

    auto ta1_z_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 351);

    auto ta1_z_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 352);

    auto ta1_z_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 353);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyyy_xx_0,  \
                             ta1_z_yyyy_xx_1,  \
                             ta1_z_yyyy_xy_0,  \
                             ta1_z_yyyy_xy_1,  \
                             ta1_z_yyyy_y_0,   \
                             ta1_z_yyyy_y_1,   \
                             ta1_z_yyyy_yy_0,  \
                             ta1_z_yyyy_yy_1,  \
                             ta1_z_yyyy_yz_0,  \
                             ta1_z_yyyy_yz_1,  \
                             ta1_z_yyyyz_xx_0, \
                             ta1_z_yyyyz_xy_0, \
                             ta1_z_yyyyz_xz_0, \
                             ta1_z_yyyyz_yy_0, \
                             ta1_z_yyyyz_yz_0, \
                             ta1_z_yyyyz_zz_0, \
                             ta1_z_yyyz_xz_0,  \
                             ta1_z_yyyz_xz_1,  \
                             ta1_z_yyyz_zz_0,  \
                             ta1_z_yyyz_zz_1,  \
                             ta1_z_yyz_xz_0,   \
                             ta1_z_yyz_xz_1,   \
                             ta1_z_yyz_zz_0,   \
                             ta1_z_yyz_zz_1,   \
                             ta_yyyy_xx_1,     \
                             ta_yyyy_xy_1,     \
                             ta_yyyy_yy_1,     \
                             ta_yyyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyz_xx_0[i] = ta_yyyy_xx_1[i] + ta1_z_yyyy_xx_0[i] * pa_z[i] - ta1_z_yyyy_xx_1[i] * pc_z[i];

        ta1_z_yyyyz_xy_0[i] = ta_yyyy_xy_1[i] + ta1_z_yyyy_xy_0[i] * pa_z[i] - ta1_z_yyyy_xy_1[i] * pc_z[i];

        ta1_z_yyyyz_xz_0[i] =
            3.0 * ta1_z_yyz_xz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xz_1[i] * fe_0 + ta1_z_yyyz_xz_0[i] * pa_y[i] - ta1_z_yyyz_xz_1[i] * pc_y[i];

        ta1_z_yyyyz_yy_0[i] = ta_yyyy_yy_1[i] + ta1_z_yyyy_yy_0[i] * pa_z[i] - ta1_z_yyyy_yy_1[i] * pc_z[i];

        ta1_z_yyyyz_yz_0[i] =
            ta1_z_yyyy_y_0[i] * fe_0 - ta1_z_yyyy_y_1[i] * fe_0 + ta_yyyy_yz_1[i] + ta1_z_yyyy_yz_0[i] * pa_z[i] - ta1_z_yyyy_yz_1[i] * pc_z[i];

        ta1_z_yyyyz_zz_0[i] =
            3.0 * ta1_z_yyz_zz_0[i] * fe_0 - 3.0 * ta1_z_yyz_zz_1[i] * fe_0 + ta1_z_yyyz_zz_0[i] * pa_y[i] - ta1_z_yyyz_zz_1[i] * pc_y[i];
    }

    // Set up 354-360 components of targeted buffer : HD

    auto ta1_z_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 354);

    auto ta1_z_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 355);

    auto ta1_z_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 356);

    auto ta1_z_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 357);

    auto ta1_z_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 358);

    auto ta1_z_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 359);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyy_xy_0,   \
                             ta1_z_yyy_xy_1,   \
                             ta1_z_yyy_yy_0,   \
                             ta1_z_yyy_yy_1,   \
                             ta1_z_yyyz_xy_0,  \
                             ta1_z_yyyz_xy_1,  \
                             ta1_z_yyyz_yy_0,  \
                             ta1_z_yyyz_yy_1,  \
                             ta1_z_yyyzz_xx_0, \
                             ta1_z_yyyzz_xy_0, \
                             ta1_z_yyyzz_xz_0, \
                             ta1_z_yyyzz_yy_0, \
                             ta1_z_yyyzz_yz_0, \
                             ta1_z_yyyzz_zz_0, \
                             ta1_z_yyzz_xx_0,  \
                             ta1_z_yyzz_xx_1,  \
                             ta1_z_yyzz_xz_0,  \
                             ta1_z_yyzz_xz_1,  \
                             ta1_z_yyzz_yz_0,  \
                             ta1_z_yyzz_yz_1,  \
                             ta1_z_yyzz_z_0,   \
                             ta1_z_yyzz_z_1,   \
                             ta1_z_yyzz_zz_0,  \
                             ta1_z_yyzz_zz_1,  \
                             ta1_z_yzz_xx_0,   \
                             ta1_z_yzz_xx_1,   \
                             ta1_z_yzz_xz_0,   \
                             ta1_z_yzz_xz_1,   \
                             ta1_z_yzz_yz_0,   \
                             ta1_z_yzz_yz_1,   \
                             ta1_z_yzz_zz_0,   \
                             ta1_z_yzz_zz_1,   \
                             ta_yyyz_xy_1,     \
                             ta_yyyz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzz_xx_0[i] =
            2.0 * ta1_z_yzz_xx_0[i] * fe_0 - 2.0 * ta1_z_yzz_xx_1[i] * fe_0 + ta1_z_yyzz_xx_0[i] * pa_y[i] - ta1_z_yyzz_xx_1[i] * pc_y[i];

        ta1_z_yyyzz_xy_0[i] =
            ta1_z_yyy_xy_0[i] * fe_0 - ta1_z_yyy_xy_1[i] * fe_0 + ta_yyyz_xy_1[i] + ta1_z_yyyz_xy_0[i] * pa_z[i] - ta1_z_yyyz_xy_1[i] * pc_z[i];

        ta1_z_yyyzz_xz_0[i] =
            2.0 * ta1_z_yzz_xz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xz_1[i] * fe_0 + ta1_z_yyzz_xz_0[i] * pa_y[i] - ta1_z_yyzz_xz_1[i] * pc_y[i];

        ta1_z_yyyzz_yy_0[i] =
            ta1_z_yyy_yy_0[i] * fe_0 - ta1_z_yyy_yy_1[i] * fe_0 + ta_yyyz_yy_1[i] + ta1_z_yyyz_yy_0[i] * pa_z[i] - ta1_z_yyyz_yy_1[i] * pc_z[i];

        ta1_z_yyyzz_yz_0[i] = 2.0 * ta1_z_yzz_yz_0[i] * fe_0 - 2.0 * ta1_z_yzz_yz_1[i] * fe_0 + ta1_z_yyzz_z_0[i] * fe_0 - ta1_z_yyzz_z_1[i] * fe_0 +
                              ta1_z_yyzz_yz_0[i] * pa_y[i] - ta1_z_yyzz_yz_1[i] * pc_y[i];

        ta1_z_yyyzz_zz_0[i] =
            2.0 * ta1_z_yzz_zz_0[i] * fe_0 - 2.0 * ta1_z_yzz_zz_1[i] * fe_0 + ta1_z_yyzz_zz_0[i] * pa_y[i] - ta1_z_yyzz_zz_1[i] * pc_y[i];
    }

    // Set up 360-366 components of targeted buffer : HD

    auto ta1_z_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 360);

    auto ta1_z_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 361);

    auto ta1_z_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 362);

    auto ta1_z_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 363);

    auto ta1_z_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 364);

    auto ta1_z_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 365);

#pragma omp simd aligned(pa_y,                 \
                             pa_z,             \
                             pc_y,             \
                             pc_z,             \
                             ta1_z_yyz_xy_0,   \
                             ta1_z_yyz_xy_1,   \
                             ta1_z_yyz_yy_0,   \
                             ta1_z_yyz_yy_1,   \
                             ta1_z_yyzz_xy_0,  \
                             ta1_z_yyzz_xy_1,  \
                             ta1_z_yyzz_yy_0,  \
                             ta1_z_yyzz_yy_1,  \
                             ta1_z_yyzzz_xx_0, \
                             ta1_z_yyzzz_xy_0, \
                             ta1_z_yyzzz_xz_0, \
                             ta1_z_yyzzz_yy_0, \
                             ta1_z_yyzzz_yz_0, \
                             ta1_z_yyzzz_zz_0, \
                             ta1_z_yzzz_xx_0,  \
                             ta1_z_yzzz_xx_1,  \
                             ta1_z_yzzz_xz_0,  \
                             ta1_z_yzzz_xz_1,  \
                             ta1_z_yzzz_yz_0,  \
                             ta1_z_yzzz_yz_1,  \
                             ta1_z_yzzz_z_0,   \
                             ta1_z_yzzz_z_1,   \
                             ta1_z_yzzz_zz_0,  \
                             ta1_z_yzzz_zz_1,  \
                             ta1_z_zzz_xx_0,   \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta_yyzz_xy_1,     \
                             ta_yyzz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzz_xx_0[i] = ta1_z_zzz_xx_0[i] * fe_0 - ta1_z_zzz_xx_1[i] * fe_0 + ta1_z_yzzz_xx_0[i] * pa_y[i] - ta1_z_yzzz_xx_1[i] * pc_y[i];

        ta1_z_yyzzz_xy_0[i] = 2.0 * ta1_z_yyz_xy_0[i] * fe_0 - 2.0 * ta1_z_yyz_xy_1[i] * fe_0 + ta_yyzz_xy_1[i] + ta1_z_yyzz_xy_0[i] * pa_z[i] -
                              ta1_z_yyzz_xy_1[i] * pc_z[i];

        ta1_z_yyzzz_xz_0[i] = ta1_z_zzz_xz_0[i] * fe_0 - ta1_z_zzz_xz_1[i] * fe_0 + ta1_z_yzzz_xz_0[i] * pa_y[i] - ta1_z_yzzz_xz_1[i] * pc_y[i];

        ta1_z_yyzzz_yy_0[i] = 2.0 * ta1_z_yyz_yy_0[i] * fe_0 - 2.0 * ta1_z_yyz_yy_1[i] * fe_0 + ta_yyzz_yy_1[i] + ta1_z_yyzz_yy_0[i] * pa_z[i] -
                              ta1_z_yyzz_yy_1[i] * pc_z[i];

        ta1_z_yyzzz_yz_0[i] = ta1_z_zzz_yz_0[i] * fe_0 - ta1_z_zzz_yz_1[i] * fe_0 + ta1_z_yzzz_z_0[i] * fe_0 - ta1_z_yzzz_z_1[i] * fe_0 +
                              ta1_z_yzzz_yz_0[i] * pa_y[i] - ta1_z_yzzz_yz_1[i] * pc_y[i];

        ta1_z_yyzzz_zz_0[i] = ta1_z_zzz_zz_0[i] * fe_0 - ta1_z_zzz_zz_1[i] * fe_0 + ta1_z_yzzz_zz_0[i] * pa_y[i] - ta1_z_yzzz_zz_1[i] * pc_y[i];
    }

    // Set up 366-372 components of targeted buffer : HD

    auto ta1_z_yzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 366);

    auto ta1_z_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 367);

    auto ta1_z_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 368);

    auto ta1_z_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 369);

    auto ta1_z_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 370);

    auto ta1_z_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 371);

#pragma omp simd aligned(pa_y,                 \
                             pc_y,             \
                             ta1_z_yzzzz_xx_0, \
                             ta1_z_yzzzz_xy_0, \
                             ta1_z_yzzzz_xz_0, \
                             ta1_z_yzzzz_yy_0, \
                             ta1_z_yzzzz_yz_0, \
                             ta1_z_yzzzz_zz_0, \
                             ta1_z_zzzz_x_0,   \
                             ta1_z_zzzz_x_1,   \
                             ta1_z_zzzz_xx_0,  \
                             ta1_z_zzzz_xx_1,  \
                             ta1_z_zzzz_xy_0,  \
                             ta1_z_zzzz_xy_1,  \
                             ta1_z_zzzz_xz_0,  \
                             ta1_z_zzzz_xz_1,  \
                             ta1_z_zzzz_y_0,   \
                             ta1_z_zzzz_y_1,   \
                             ta1_z_zzzz_yy_0,  \
                             ta1_z_zzzz_yy_1,  \
                             ta1_z_zzzz_yz_0,  \
                             ta1_z_zzzz_yz_1,  \
                             ta1_z_zzzz_z_0,   \
                             ta1_z_zzzz_z_1,   \
                             ta1_z_zzzz_zz_0,  \
                             ta1_z_zzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzz_xx_0[i] = ta1_z_zzzz_xx_0[i] * pa_y[i] - ta1_z_zzzz_xx_1[i] * pc_y[i];

        ta1_z_yzzzz_xy_0[i] = ta1_z_zzzz_x_0[i] * fe_0 - ta1_z_zzzz_x_1[i] * fe_0 + ta1_z_zzzz_xy_0[i] * pa_y[i] - ta1_z_zzzz_xy_1[i] * pc_y[i];

        ta1_z_yzzzz_xz_0[i] = ta1_z_zzzz_xz_0[i] * pa_y[i] - ta1_z_zzzz_xz_1[i] * pc_y[i];

        ta1_z_yzzzz_yy_0[i] =
            2.0 * ta1_z_zzzz_y_0[i] * fe_0 - 2.0 * ta1_z_zzzz_y_1[i] * fe_0 + ta1_z_zzzz_yy_0[i] * pa_y[i] - ta1_z_zzzz_yy_1[i] * pc_y[i];

        ta1_z_yzzzz_yz_0[i] = ta1_z_zzzz_z_0[i] * fe_0 - ta1_z_zzzz_z_1[i] * fe_0 + ta1_z_zzzz_yz_0[i] * pa_y[i] - ta1_z_zzzz_yz_1[i] * pc_y[i];

        ta1_z_yzzzz_zz_0[i] = ta1_z_zzzz_zz_0[i] * pa_y[i] - ta1_z_zzzz_zz_1[i] * pc_y[i];
    }

    // Set up 372-378 components of targeted buffer : HD

    auto ta1_z_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 372);

    auto ta1_z_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 373);

    auto ta1_z_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 374);

    auto ta1_z_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 375);

    auto ta1_z_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 376);

    auto ta1_z_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 377);

#pragma omp simd aligned(pa_z,                 \
                             pc_z,             \
                             ta1_z_zzz_xx_0,   \
                             ta1_z_zzz_xx_1,   \
                             ta1_z_zzz_xy_0,   \
                             ta1_z_zzz_xy_1,   \
                             ta1_z_zzz_xz_0,   \
                             ta1_z_zzz_xz_1,   \
                             ta1_z_zzz_yy_0,   \
                             ta1_z_zzz_yy_1,   \
                             ta1_z_zzz_yz_0,   \
                             ta1_z_zzz_yz_1,   \
                             ta1_z_zzz_zz_0,   \
                             ta1_z_zzz_zz_1,   \
                             ta1_z_zzzz_x_0,   \
                             ta1_z_zzzz_x_1,   \
                             ta1_z_zzzz_xx_0,  \
                             ta1_z_zzzz_xx_1,  \
                             ta1_z_zzzz_xy_0,  \
                             ta1_z_zzzz_xy_1,  \
                             ta1_z_zzzz_xz_0,  \
                             ta1_z_zzzz_xz_1,  \
                             ta1_z_zzzz_y_0,   \
                             ta1_z_zzzz_y_1,   \
                             ta1_z_zzzz_yy_0,  \
                             ta1_z_zzzz_yy_1,  \
                             ta1_z_zzzz_yz_0,  \
                             ta1_z_zzzz_yz_1,  \
                             ta1_z_zzzz_z_0,   \
                             ta1_z_zzzz_z_1,   \
                             ta1_z_zzzz_zz_0,  \
                             ta1_z_zzzz_zz_1,  \
                             ta1_z_zzzzz_xx_0, \
                             ta1_z_zzzzz_xy_0, \
                             ta1_z_zzzzz_xz_0, \
                             ta1_z_zzzzz_yy_0, \
                             ta1_z_zzzzz_yz_0, \
                             ta1_z_zzzzz_zz_0, \
                             ta_zzzz_xx_1,     \
                             ta_zzzz_xy_1,     \
                             ta_zzzz_xz_1,     \
                             ta_zzzz_yy_1,     \
                             ta_zzzz_yz_1,     \
                             ta_zzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzz_xx_0[i] = 4.0 * ta1_z_zzz_xx_0[i] * fe_0 - 4.0 * ta1_z_zzz_xx_1[i] * fe_0 + ta_zzzz_xx_1[i] + ta1_z_zzzz_xx_0[i] * pa_z[i] -
                              ta1_z_zzzz_xx_1[i] * pc_z[i];

        ta1_z_zzzzz_xy_0[i] = 4.0 * ta1_z_zzz_xy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xy_1[i] * fe_0 + ta_zzzz_xy_1[i] + ta1_z_zzzz_xy_0[i] * pa_z[i] -
                              ta1_z_zzzz_xy_1[i] * pc_z[i];

        ta1_z_zzzzz_xz_0[i] = 4.0 * ta1_z_zzz_xz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xz_1[i] * fe_0 + ta1_z_zzzz_x_0[i] * fe_0 - ta1_z_zzzz_x_1[i] * fe_0 +
                              ta_zzzz_xz_1[i] + ta1_z_zzzz_xz_0[i] * pa_z[i] - ta1_z_zzzz_xz_1[i] * pc_z[i];

        ta1_z_zzzzz_yy_0[i] = 4.0 * ta1_z_zzz_yy_0[i] * fe_0 - 4.0 * ta1_z_zzz_yy_1[i] * fe_0 + ta_zzzz_yy_1[i] + ta1_z_zzzz_yy_0[i] * pa_z[i] -
                              ta1_z_zzzz_yy_1[i] * pc_z[i];

        ta1_z_zzzzz_yz_0[i] = 4.0 * ta1_z_zzz_yz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yz_1[i] * fe_0 + ta1_z_zzzz_y_0[i] * fe_0 - ta1_z_zzzz_y_1[i] * fe_0 +
                              ta_zzzz_yz_1[i] + ta1_z_zzzz_yz_0[i] * pa_z[i] - ta1_z_zzzz_yz_1[i] * pc_z[i];

        ta1_z_zzzzz_zz_0[i] = 4.0 * ta1_z_zzz_zz_0[i] * fe_0 - 4.0 * ta1_z_zzz_zz_1[i] * fe_0 + 2.0 * ta1_z_zzzz_z_0[i] * fe_0 -
                              2.0 * ta1_z_zzzz_z_1[i] * fe_0 + ta_zzzz_zz_1[i] + ta1_z_zzzz_zz_0[i] * pa_z[i] - ta1_z_zzzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
