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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXGG.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xgg(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xgg,
                                        const size_t idx_xfg,
                                        const size_t idx_geom_10_xfg,
                                        const size_t idx_geom_10_xfh,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_cd,
                                        const int a_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});

    // Set up R(CD) distances

    auto cd_x = factors.data(idx_cd);

    auto cd_y = factors.data(idx_cd + 1);

    auto cd_z = factors.data(idx_cd + 2);

    for (int i = 0; i < acomps; i++)
    {
        /// Set up components of auxilary buffer : SFG

        const auto fg_off = idx_xfg + i * 150;

        auto g_xxx_xxxx = cbuffer.data(fg_off + 0);

        auto g_xxx_xxxy = cbuffer.data(fg_off + 1);

        auto g_xxx_xxxz = cbuffer.data(fg_off + 2);

        auto g_xxx_xxyy = cbuffer.data(fg_off + 3);

        auto g_xxx_xxyz = cbuffer.data(fg_off + 4);

        auto g_xxx_xxzz = cbuffer.data(fg_off + 5);

        auto g_xxx_xyyy = cbuffer.data(fg_off + 6);

        auto g_xxx_xyyz = cbuffer.data(fg_off + 7);

        auto g_xxx_xyzz = cbuffer.data(fg_off + 8);

        auto g_xxx_xzzz = cbuffer.data(fg_off + 9);

        auto g_xxx_yyyy = cbuffer.data(fg_off + 10);

        auto g_xxx_yyyz = cbuffer.data(fg_off + 11);

        auto g_xxx_yyzz = cbuffer.data(fg_off + 12);

        auto g_xxx_yzzz = cbuffer.data(fg_off + 13);

        auto g_xxx_zzzz = cbuffer.data(fg_off + 14);

        auto g_yyy_xxxx = cbuffer.data(fg_off + 90);

        auto g_yyy_xxxy = cbuffer.data(fg_off + 91);

        auto g_yyy_xxxz = cbuffer.data(fg_off + 92);

        auto g_yyy_xxyy = cbuffer.data(fg_off + 93);

        auto g_yyy_xxyz = cbuffer.data(fg_off + 94);

        auto g_yyy_xxzz = cbuffer.data(fg_off + 95);

        auto g_yyy_xyyy = cbuffer.data(fg_off + 96);

        auto g_yyy_xyyz = cbuffer.data(fg_off + 97);

        auto g_yyy_xyzz = cbuffer.data(fg_off + 98);

        auto g_yyy_xzzz = cbuffer.data(fg_off + 99);

        auto g_yyy_yyyy = cbuffer.data(fg_off + 100);

        auto g_yyy_yyyz = cbuffer.data(fg_off + 101);

        auto g_yyy_yyzz = cbuffer.data(fg_off + 102);

        auto g_yyy_yzzz = cbuffer.data(fg_off + 103);

        auto g_yyy_zzzz = cbuffer.data(fg_off + 104);

        auto g_zzz_xxxx = cbuffer.data(fg_off + 135);

        auto g_zzz_xxxy = cbuffer.data(fg_off + 136);

        auto g_zzz_xxxz = cbuffer.data(fg_off + 137);

        auto g_zzz_xxyy = cbuffer.data(fg_off + 138);

        auto g_zzz_xxyz = cbuffer.data(fg_off + 139);

        auto g_zzz_xxzz = cbuffer.data(fg_off + 140);

        auto g_zzz_xyyy = cbuffer.data(fg_off + 141);

        auto g_zzz_xyyz = cbuffer.data(fg_off + 142);

        auto g_zzz_xyzz = cbuffer.data(fg_off + 143);

        auto g_zzz_xzzz = cbuffer.data(fg_off + 144);

        auto g_zzz_yyyy = cbuffer.data(fg_off + 145);

        auto g_zzz_yyyz = cbuffer.data(fg_off + 146);

        auto g_zzz_yyzz = cbuffer.data(fg_off + 147);

        auto g_zzz_yzzz = cbuffer.data(fg_off + 148);

        auto g_zzz_zzzz = cbuffer.data(fg_off + 149);

        /// Set up components of auxilary buffer : SFG

        const auto fg_geom_10_off = idx_geom_10_xfg + i * 150;

        auto g_x_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 30);

        auto g_x_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 35);

        auto g_x_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 41);

        auto g_x_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 42);

        auto g_x_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 48);

        auto g_x_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 55);

        auto g_x_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 56);

        auto g_x_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 59);

        auto g_x_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 60);

        auto g_x_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 61);

        auto g_x_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 62);

        auto g_x_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 63);

        auto g_x_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 64);

        auto g_x_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 65);

        auto g_x_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 66);

        auto g_x_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 67);

        auto g_x_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 68);

        auto g_x_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 69);

        auto g_x_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 70);

        auto g_x_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 71);

        auto g_x_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 72);

        auto g_x_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 73);

        auto g_x_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 74);

        auto g_x_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 75);

        auto g_x_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 76);

        auto g_x_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 77);

        auto g_x_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 78);

        auto g_x_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 79);

        auto g_x_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 80);

        auto g_x_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 81);

        auto g_x_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 82);

        auto g_x_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 83);

        auto g_x_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 84);

        auto g_x_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 85);

        auto g_x_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 86);

        auto g_x_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 87);

        auto g_x_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 88);

        auto g_x_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 89);

        auto g_x_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 90);

        auto g_x_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 91);

        auto g_x_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 92);

        auto g_x_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 93);

        auto g_x_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 94);

        auto g_x_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 95);

        auto g_x_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 96);

        auto g_x_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 97);

        auto g_x_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 98);

        auto g_x_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 99);

        auto g_x_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 100);

        auto g_x_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 101);

        auto g_x_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 102);

        auto g_x_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 103);

        auto g_x_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 104);

        auto g_x_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 105);

        auto g_x_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 106);

        auto g_x_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 107);

        auto g_x_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 108);

        auto g_x_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 109);

        auto g_x_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 110);

        auto g_x_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 111);

        auto g_x_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 112);

        auto g_x_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 113);

        auto g_x_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 114);

        auto g_x_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 115);

        auto g_x_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 116);

        auto g_x_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 117);

        auto g_x_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 118);

        auto g_x_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 119);

        auto g_x_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 120);

        auto g_x_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 121);

        auto g_x_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 122);

        auto g_x_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 123);

        auto g_x_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 124);

        auto g_x_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 125);

        auto g_x_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 126);

        auto g_x_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 127);

        auto g_x_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 128);

        auto g_x_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 129);

        auto g_x_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 130);

        auto g_x_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 131);

        auto g_x_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 132);

        auto g_x_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 133);

        auto g_x_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 134);

        auto g_x_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 0 * acomps + 135);

        auto g_x_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 0 * acomps + 136);

        auto g_x_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 0 * acomps + 137);

        auto g_x_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 138);

        auto g_x_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 139);

        auto g_x_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 140);

        auto g_x_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 141);

        auto g_x_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 142);

        auto g_x_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 143);

        auto g_x_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 144);

        auto g_x_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 0 * acomps + 145);

        auto g_x_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 0 * acomps + 146);

        auto g_x_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 147);

        auto g_x_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 148);

        auto g_x_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 0 * acomps + 149);

        auto g_y_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 0);

        auto g_y_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 1);

        auto g_y_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 2);

        auto g_y_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 3);

        auto g_y_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 4);

        auto g_y_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 5);

        auto g_y_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 6);

        auto g_y_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 7);

        auto g_y_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 8);

        auto g_y_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 9);

        auto g_y_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 10);

        auto g_y_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 11);

        auto g_y_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 12);

        auto g_y_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 13);

        auto g_y_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 14);

        auto g_y_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 15);

        auto g_y_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 16);

        auto g_y_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 17);

        auto g_y_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 18);

        auto g_y_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 19);

        auto g_y_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 20);

        auto g_y_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 21);

        auto g_y_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 22);

        auto g_y_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 23);

        auto g_y_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 24);

        auto g_y_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 25);

        auto g_y_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 26);

        auto g_y_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 27);

        auto g_y_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 28);

        auto g_y_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 29);

        auto g_y_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 30);

        auto g_y_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 31);

        auto g_y_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 32);

        auto g_y_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 33);

        auto g_y_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 34);

        auto g_y_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 35);

        auto g_y_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 36);

        auto g_y_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 37);

        auto g_y_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 38);

        auto g_y_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 39);

        auto g_y_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 40);

        auto g_y_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 41);

        auto g_y_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 42);

        auto g_y_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 43);

        auto g_y_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 44);

        auto g_y_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 45);

        auto g_y_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 46);

        auto g_y_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 47);

        auto g_y_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 48);

        auto g_y_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 49);

        auto g_y_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 50);

        auto g_y_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 51);

        auto g_y_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 52);

        auto g_y_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 53);

        auto g_y_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 54);

        auto g_y_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 55);

        auto g_y_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 56);

        auto g_y_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 57);

        auto g_y_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 58);

        auto g_y_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 59);

        auto g_y_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 60);

        auto g_y_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 61);

        auto g_y_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 62);

        auto g_y_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 63);

        auto g_y_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 64);

        auto g_y_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 65);

        auto g_y_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 66);

        auto g_y_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 67);

        auto g_y_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 68);

        auto g_y_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 69);

        auto g_y_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 70);

        auto g_y_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 71);

        auto g_y_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 72);

        auto g_y_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 73);

        auto g_y_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 74);

        auto g_y_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 75);

        auto g_y_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 76);

        auto g_y_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 77);

        auto g_y_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 78);

        auto g_y_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 79);

        auto g_y_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 80);

        auto g_y_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 81);

        auto g_y_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 82);

        auto g_y_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 83);

        auto g_y_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 84);

        auto g_y_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 85);

        auto g_y_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 86);

        auto g_y_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 87);

        auto g_y_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 88);

        auto g_y_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 89);

        auto g_y_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 90);

        auto g_y_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 91);

        auto g_y_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 92);

        auto g_y_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 93);

        auto g_y_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 94);

        auto g_y_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 95);

        auto g_y_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 96);

        auto g_y_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 97);

        auto g_y_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 98);

        auto g_y_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 99);

        auto g_y_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 100);

        auto g_y_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 101);

        auto g_y_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 102);

        auto g_y_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 103);

        auto g_y_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 104);

        auto g_y_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 105);

        auto g_y_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 106);

        auto g_y_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 107);

        auto g_y_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 108);

        auto g_y_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 109);

        auto g_y_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 110);

        auto g_y_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 111);

        auto g_y_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 112);

        auto g_y_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 113);

        auto g_y_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 114);

        auto g_y_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 115);

        auto g_y_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 116);

        auto g_y_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 117);

        auto g_y_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 118);

        auto g_y_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 119);

        auto g_y_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 120);

        auto g_y_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 121);

        auto g_y_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 122);

        auto g_y_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 123);

        auto g_y_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 124);

        auto g_y_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 125);

        auto g_y_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 126);

        auto g_y_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 127);

        auto g_y_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 128);

        auto g_y_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 129);

        auto g_y_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 130);

        auto g_y_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 131);

        auto g_y_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 132);

        auto g_y_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 133);

        auto g_y_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 134);

        auto g_y_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 150 * acomps + 135);

        auto g_y_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 150 * acomps + 136);

        auto g_y_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 150 * acomps + 137);

        auto g_y_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 138);

        auto g_y_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 139);

        auto g_y_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 140);

        auto g_y_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 141);

        auto g_y_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 142);

        auto g_y_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 143);

        auto g_y_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 144);

        auto g_y_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 150 * acomps + 145);

        auto g_y_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 150 * acomps + 146);

        auto g_y_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 147);

        auto g_y_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 148);

        auto g_y_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 150 * acomps + 149);

        auto g_z_0_xxx_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 0);

        auto g_z_0_xxx_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 1);

        auto g_z_0_xxx_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 2);

        auto g_z_0_xxx_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 3);

        auto g_z_0_xxx_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 4);

        auto g_z_0_xxx_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 5);

        auto g_z_0_xxx_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 6);

        auto g_z_0_xxx_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 7);

        auto g_z_0_xxx_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 8);

        auto g_z_0_xxx_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 9);

        auto g_z_0_xxx_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 10);

        auto g_z_0_xxx_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 11);

        auto g_z_0_xxx_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 12);

        auto g_z_0_xxx_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 13);

        auto g_z_0_xxx_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 14);

        auto g_z_0_xxy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 15);

        auto g_z_0_xxy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 16);

        auto g_z_0_xxy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 17);

        auto g_z_0_xxy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 18);

        auto g_z_0_xxy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 19);

        auto g_z_0_xxy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 20);

        auto g_z_0_xxy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 21);

        auto g_z_0_xxy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 22);

        auto g_z_0_xxy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 23);

        auto g_z_0_xxy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 24);

        auto g_z_0_xxy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 25);

        auto g_z_0_xxy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 26);

        auto g_z_0_xxy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 27);

        auto g_z_0_xxy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 28);

        auto g_z_0_xxy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 29);

        auto g_z_0_xxz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 30);

        auto g_z_0_xxz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 31);

        auto g_z_0_xxz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 32);

        auto g_z_0_xxz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 33);

        auto g_z_0_xxz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 34);

        auto g_z_0_xxz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 35);

        auto g_z_0_xxz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 36);

        auto g_z_0_xxz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 37);

        auto g_z_0_xxz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 38);

        auto g_z_0_xxz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 39);

        auto g_z_0_xxz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 40);

        auto g_z_0_xxz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 41);

        auto g_z_0_xxz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 42);

        auto g_z_0_xxz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 43);

        auto g_z_0_xxz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 44);

        auto g_z_0_xyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 45);

        auto g_z_0_xyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 46);

        auto g_z_0_xyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 47);

        auto g_z_0_xyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 48);

        auto g_z_0_xyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 49);

        auto g_z_0_xyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 50);

        auto g_z_0_xyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 51);

        auto g_z_0_xyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 52);

        auto g_z_0_xyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 53);

        auto g_z_0_xyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 54);

        auto g_z_0_xyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 55);

        auto g_z_0_xyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 56);

        auto g_z_0_xyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 57);

        auto g_z_0_xyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 58);

        auto g_z_0_xyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 59);

        auto g_z_0_xyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 60);

        auto g_z_0_xyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 61);

        auto g_z_0_xyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 62);

        auto g_z_0_xyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 63);

        auto g_z_0_xyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 64);

        auto g_z_0_xyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 65);

        auto g_z_0_xyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 66);

        auto g_z_0_xyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 67);

        auto g_z_0_xyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 68);

        auto g_z_0_xyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 69);

        auto g_z_0_xyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 70);

        auto g_z_0_xyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 71);

        auto g_z_0_xyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 72);

        auto g_z_0_xyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 73);

        auto g_z_0_xyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 74);

        auto g_z_0_xzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 75);

        auto g_z_0_xzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 76);

        auto g_z_0_xzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 77);

        auto g_z_0_xzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 78);

        auto g_z_0_xzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 79);

        auto g_z_0_xzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 80);

        auto g_z_0_xzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 81);

        auto g_z_0_xzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 82);

        auto g_z_0_xzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 83);

        auto g_z_0_xzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 84);

        auto g_z_0_xzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 85);

        auto g_z_0_xzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 86);

        auto g_z_0_xzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 87);

        auto g_z_0_xzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 88);

        auto g_z_0_xzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 89);

        auto g_z_0_yyy_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 90);

        auto g_z_0_yyy_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 91);

        auto g_z_0_yyy_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 92);

        auto g_z_0_yyy_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 93);

        auto g_z_0_yyy_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 94);

        auto g_z_0_yyy_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 95);

        auto g_z_0_yyy_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 96);

        auto g_z_0_yyy_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 97);

        auto g_z_0_yyy_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 98);

        auto g_z_0_yyy_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 99);

        auto g_z_0_yyy_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 100);

        auto g_z_0_yyy_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 101);

        auto g_z_0_yyy_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 102);

        auto g_z_0_yyy_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 103);

        auto g_z_0_yyy_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 104);

        auto g_z_0_yyz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 105);

        auto g_z_0_yyz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 106);

        auto g_z_0_yyz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 107);

        auto g_z_0_yyz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 108);

        auto g_z_0_yyz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 109);

        auto g_z_0_yyz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 110);

        auto g_z_0_yyz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 111);

        auto g_z_0_yyz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 112);

        auto g_z_0_yyz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 113);

        auto g_z_0_yyz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 114);

        auto g_z_0_yyz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 115);

        auto g_z_0_yyz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 116);

        auto g_z_0_yyz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 117);

        auto g_z_0_yyz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 118);

        auto g_z_0_yyz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 119);

        auto g_z_0_yzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 120);

        auto g_z_0_yzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 121);

        auto g_z_0_yzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 122);

        auto g_z_0_yzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 123);

        auto g_z_0_yzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 124);

        auto g_z_0_yzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 125);

        auto g_z_0_yzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 126);

        auto g_z_0_yzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 127);

        auto g_z_0_yzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 128);

        auto g_z_0_yzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 129);

        auto g_z_0_yzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 130);

        auto g_z_0_yzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 131);

        auto g_z_0_yzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 132);

        auto g_z_0_yzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 133);

        auto g_z_0_yzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 134);

        auto g_z_0_zzz_xxxx = cbuffer.data(fg_geom_10_off + 300 * acomps + 135);

        auto g_z_0_zzz_xxxy = cbuffer.data(fg_geom_10_off + 300 * acomps + 136);

        auto g_z_0_zzz_xxxz = cbuffer.data(fg_geom_10_off + 300 * acomps + 137);

        auto g_z_0_zzz_xxyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 138);

        auto g_z_0_zzz_xxyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 139);

        auto g_z_0_zzz_xxzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 140);

        auto g_z_0_zzz_xyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 141);

        auto g_z_0_zzz_xyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 142);

        auto g_z_0_zzz_xyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 143);

        auto g_z_0_zzz_xzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 144);

        auto g_z_0_zzz_yyyy = cbuffer.data(fg_geom_10_off + 300 * acomps + 145);

        auto g_z_0_zzz_yyyz = cbuffer.data(fg_geom_10_off + 300 * acomps + 146);

        auto g_z_0_zzz_yyzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 147);

        auto g_z_0_zzz_yzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 148);

        auto g_z_0_zzz_zzzz = cbuffer.data(fg_geom_10_off + 300 * acomps + 149);

        /// Set up components of auxilary buffer : SFH

        const auto fh_geom_10_off = idx_geom_10_xfh + i * 210;

        auto g_x_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xxx_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xxx_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xxx_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xxx_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xxx_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 19);

        auto g_x_0_xxx_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 20);

        auto g_x_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 29);

        auto g_x_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 31);

        auto g_x_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 32);

        auto g_x_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 33);

        auto g_x_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 34);

        auto g_x_0_xxy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 36);

        auto g_x_0_xxy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 37);

        auto g_x_0_xxy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 38);

        auto g_x_0_xxy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 39);

        auto g_x_0_xxy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 40);

        auto g_x_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 43);

        auto g_x_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps + 44);

        auto g_x_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 45);

        auto g_x_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 46);

        auto g_x_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 47);

        auto g_x_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 48);

        auto g_x_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 49);

        auto g_x_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 50);

        auto g_x_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 51);

        auto g_x_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 52);

        auto g_x_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 53);

        auto g_x_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 54);

        auto g_x_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 55);

        auto g_x_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 56);

        auto g_x_0_xxz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 57);

        auto g_x_0_xxz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 58);

        auto g_x_0_xxz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 59);

        auto g_x_0_xxz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 60);

        auto g_x_0_xxz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 61);

        auto g_x_0_xxz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 62);

        auto g_x_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 64);

        auto g_x_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 66);

        auto g_x_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 67);

        auto g_x_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 69);

        auto g_x_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 70);

        auto g_x_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 71);

        auto g_x_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 73);

        auto g_x_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 74);

        auto g_x_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 75);

        auto g_x_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 76);

        auto g_x_0_xyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 78);

        auto g_x_0_xyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 79);

        auto g_x_0_xyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 80);

        auto g_x_0_xyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 81);

        auto g_x_0_xyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 82);

        auto g_x_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 85);

        auto g_x_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 87);

        auto g_x_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 88);

        auto g_x_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 90);

        auto g_x_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 91);

        auto g_x_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 92);

        auto g_x_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 94);

        auto g_x_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 95);

        auto g_x_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 96);

        auto g_x_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 97);

        auto g_x_0_xyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 99);

        auto g_x_0_xyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 100);

        auto g_x_0_xyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 101);

        auto g_x_0_xyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 102);

        auto g_x_0_xyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 103);

        auto g_x_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 106);

        auto g_x_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps + 107);

        auto g_x_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 108);

        auto g_x_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 109);

        auto g_x_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 110);

        auto g_x_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 111);

        auto g_x_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 112);

        auto g_x_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 113);

        auto g_x_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 114);

        auto g_x_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 115);

        auto g_x_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 116);

        auto g_x_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 117);

        auto g_x_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 118);

        auto g_x_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 119);

        auto g_x_0_xzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 120);

        auto g_x_0_xzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 121);

        auto g_x_0_xzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 122);

        auto g_x_0_xzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 123);

        auto g_x_0_xzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 124);

        auto g_x_0_xzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 125);

        auto g_x_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 127);

        auto g_x_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 129);

        auto g_x_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 130);

        auto g_x_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 132);

        auto g_x_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 133);

        auto g_x_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 134);

        auto g_x_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 136);

        auto g_x_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 137);

        auto g_x_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 138);

        auto g_x_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 139);

        auto g_x_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 141);

        auto g_x_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 142);

        auto g_x_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 143);

        auto g_x_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 144);

        auto g_x_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 145);

        auto g_x_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 148);

        auto g_x_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 150);

        auto g_x_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 151);

        auto g_x_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 153);

        auto g_x_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 154);

        auto g_x_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 155);

        auto g_x_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 157);

        auto g_x_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 158);

        auto g_x_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 159);

        auto g_x_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 160);

        auto g_x_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 162);

        auto g_x_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 163);

        auto g_x_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 164);

        auto g_x_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 165);

        auto g_x_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 166);

        auto g_x_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 169);

        auto g_x_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 171);

        auto g_x_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 172);

        auto g_x_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 174);

        auto g_x_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 175);

        auto g_x_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 176);

        auto g_x_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 178);

        auto g_x_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 179);

        auto g_x_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 180);

        auto g_x_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 181);

        auto g_x_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 183);

        auto g_x_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 184);

        auto g_x_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 185);

        auto g_x_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 186);

        auto g_x_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 187);

        auto g_x_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 0 * acomps + 190);

        auto g_x_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 0 * acomps + 191);

        auto g_x_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 192);

        auto g_x_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 193);

        auto g_x_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 194);

        auto g_x_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 195);

        auto g_x_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 196);

        auto g_x_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 197);

        auto g_x_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 198);

        auto g_x_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 199);

        auto g_x_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 200);

        auto g_x_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 201);

        auto g_x_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 202);

        auto g_x_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 203);

        auto g_x_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 0 * acomps + 204);

        auto g_x_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 0 * acomps + 205);

        auto g_x_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 206);

        auto g_x_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 207);

        auto g_x_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 208);

        auto g_x_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 0 * acomps + 209);

        auto g_y_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 0);

        auto g_y_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 1);

        auto g_y_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 2);

        auto g_y_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 3);

        auto g_y_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 4);

        auto g_y_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 5);

        auto g_y_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 6);

        auto g_y_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 7);

        auto g_y_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 8);

        auto g_y_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 9);

        auto g_y_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 10);

        auto g_y_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 11);

        auto g_y_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 12);

        auto g_y_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 13);

        auto g_y_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 14);

        auto g_y_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 21);

        auto g_y_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 22);

        auto g_y_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 23);

        auto g_y_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 24);

        auto g_y_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 25);

        auto g_y_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 26);

        auto g_y_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 27);

        auto g_y_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 28);

        auto g_y_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 29);

        auto g_y_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 30);

        auto g_y_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 31);

        auto g_y_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 32);

        auto g_y_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 33);

        auto g_y_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 34);

        auto g_y_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 35);

        auto g_y_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 42);

        auto g_y_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 43);

        auto g_y_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 44);

        auto g_y_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 45);

        auto g_y_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 46);

        auto g_y_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 47);

        auto g_y_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 48);

        auto g_y_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 49);

        auto g_y_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 50);

        auto g_y_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 51);

        auto g_y_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 52);

        auto g_y_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 53);

        auto g_y_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 54);

        auto g_y_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 55);

        auto g_y_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 56);

        auto g_y_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 63);

        auto g_y_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 64);

        auto g_y_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 65);

        auto g_y_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 66);

        auto g_y_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 67);

        auto g_y_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 68);

        auto g_y_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 69);

        auto g_y_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 70);

        auto g_y_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 71);

        auto g_y_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 72);

        auto g_y_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 73);

        auto g_y_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 74);

        auto g_y_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 75);

        auto g_y_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 76);

        auto g_y_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 77);

        auto g_y_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 84);

        auto g_y_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 85);

        auto g_y_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 86);

        auto g_y_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 87);

        auto g_y_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 88);

        auto g_y_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 89);

        auto g_y_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 90);

        auto g_y_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 91);

        auto g_y_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 92);

        auto g_y_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 93);

        auto g_y_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 94);

        auto g_y_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 95);

        auto g_y_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 96);

        auto g_y_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 97);

        auto g_y_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 98);

        auto g_y_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 105);

        auto g_y_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 106);

        auto g_y_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 107);

        auto g_y_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 108);

        auto g_y_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 109);

        auto g_y_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 110);

        auto g_y_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 111);

        auto g_y_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 112);

        auto g_y_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 113);

        auto g_y_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 114);

        auto g_y_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 115);

        auto g_y_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 116);

        auto g_y_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 117);

        auto g_y_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 118);

        auto g_y_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 119);

        auto g_y_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 126);

        auto g_y_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 127);

        auto g_y_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 128);

        auto g_y_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 129);

        auto g_y_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 130);

        auto g_y_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 131);

        auto g_y_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 132);

        auto g_y_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 133);

        auto g_y_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 134);

        auto g_y_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 135);

        auto g_y_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 136);

        auto g_y_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 137);

        auto g_y_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 138);

        auto g_y_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 139);

        auto g_y_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 140);

        auto g_y_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 141);

        auto g_y_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 142);

        auto g_y_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 143);

        auto g_y_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 144);

        auto g_y_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 145);

        auto g_y_0_yyy_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 146);

        auto g_y_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 147);

        auto g_y_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 148);

        auto g_y_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 149);

        auto g_y_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 150);

        auto g_y_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 151);

        auto g_y_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 152);

        auto g_y_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 153);

        auto g_y_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 154);

        auto g_y_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 155);

        auto g_y_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 156);

        auto g_y_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 157);

        auto g_y_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 158);

        auto g_y_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 159);

        auto g_y_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 160);

        auto g_y_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 161);

        auto g_y_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 163);

        auto g_y_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 164);

        auto g_y_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 165);

        auto g_y_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 166);

        auto g_y_0_yyz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 167);

        auto g_y_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 168);

        auto g_y_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 169);

        auto g_y_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 170);

        auto g_y_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 171);

        auto g_y_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 172);

        auto g_y_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 173);

        auto g_y_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 174);

        auto g_y_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 175);

        auto g_y_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 176);

        auto g_y_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 177);

        auto g_y_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 178);

        auto g_y_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 179);

        auto g_y_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 180);

        auto g_y_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 181);

        auto g_y_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 182);

        auto g_y_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 184);

        auto g_y_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 185);

        auto g_y_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 186);

        auto g_y_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 187);

        auto g_y_0_yzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 188);

        auto g_y_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 210 * acomps + 189);

        auto g_y_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 210 * acomps + 190);

        auto g_y_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 210 * acomps + 191);

        auto g_y_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 192);

        auto g_y_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 193);

        auto g_y_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 194);

        auto g_y_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 195);

        auto g_y_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 196);

        auto g_y_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 197);

        auto g_y_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 198);

        auto g_y_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 210 * acomps + 199);

        auto g_y_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 200);

        auto g_y_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 201);

        auto g_y_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 202);

        auto g_y_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 203);

        auto g_y_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 210 * acomps + 205);

        auto g_y_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 206);

        auto g_y_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 207);

        auto g_y_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 208);

        auto g_y_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 210 * acomps + 209);

        auto g_z_0_xxx_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 0);

        auto g_z_0_xxx_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 1);

        auto g_z_0_xxx_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 2);

        auto g_z_0_xxx_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 3);

        auto g_z_0_xxx_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 4);

        auto g_z_0_xxx_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 5);

        auto g_z_0_xxx_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 6);

        auto g_z_0_xxx_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 7);

        auto g_z_0_xxx_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 8);

        auto g_z_0_xxx_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 9);

        auto g_z_0_xxx_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 10);

        auto g_z_0_xxx_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 11);

        auto g_z_0_xxx_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 12);

        auto g_z_0_xxx_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 13);

        auto g_z_0_xxx_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 14);

        auto g_z_0_xxy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 21);

        auto g_z_0_xxy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 22);

        auto g_z_0_xxy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 23);

        auto g_z_0_xxy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 24);

        auto g_z_0_xxy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 25);

        auto g_z_0_xxy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 26);

        auto g_z_0_xxy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 27);

        auto g_z_0_xxy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 28);

        auto g_z_0_xxy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 29);

        auto g_z_0_xxy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 30);

        auto g_z_0_xxy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 31);

        auto g_z_0_xxy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 32);

        auto g_z_0_xxy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 33);

        auto g_z_0_xxy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 34);

        auto g_z_0_xxy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 35);

        auto g_z_0_xxz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 42);

        auto g_z_0_xxz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 43);

        auto g_z_0_xxz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 44);

        auto g_z_0_xxz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 45);

        auto g_z_0_xxz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 46);

        auto g_z_0_xxz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 47);

        auto g_z_0_xxz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 48);

        auto g_z_0_xxz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 49);

        auto g_z_0_xxz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 50);

        auto g_z_0_xxz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 51);

        auto g_z_0_xxz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 52);

        auto g_z_0_xxz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 53);

        auto g_z_0_xxz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 54);

        auto g_z_0_xxz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 55);

        auto g_z_0_xxz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 56);

        auto g_z_0_xyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 63);

        auto g_z_0_xyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 64);

        auto g_z_0_xyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 65);

        auto g_z_0_xyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 66);

        auto g_z_0_xyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 67);

        auto g_z_0_xyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 68);

        auto g_z_0_xyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 69);

        auto g_z_0_xyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 70);

        auto g_z_0_xyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 71);

        auto g_z_0_xyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 72);

        auto g_z_0_xyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 73);

        auto g_z_0_xyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 74);

        auto g_z_0_xyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 75);

        auto g_z_0_xyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 76);

        auto g_z_0_xyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 77);

        auto g_z_0_xyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 84);

        auto g_z_0_xyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 85);

        auto g_z_0_xyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 86);

        auto g_z_0_xyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 87);

        auto g_z_0_xyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 88);

        auto g_z_0_xyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 89);

        auto g_z_0_xyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 90);

        auto g_z_0_xyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 91);

        auto g_z_0_xyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 92);

        auto g_z_0_xyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 93);

        auto g_z_0_xyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 94);

        auto g_z_0_xyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 95);

        auto g_z_0_xyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 96);

        auto g_z_0_xyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 97);

        auto g_z_0_xyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 98);

        auto g_z_0_xzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 105);

        auto g_z_0_xzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 106);

        auto g_z_0_xzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 107);

        auto g_z_0_xzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 108);

        auto g_z_0_xzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 109);

        auto g_z_0_xzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 110);

        auto g_z_0_xzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 111);

        auto g_z_0_xzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 112);

        auto g_z_0_xzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 113);

        auto g_z_0_xzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 114);

        auto g_z_0_xzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 115);

        auto g_z_0_xzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 116);

        auto g_z_0_xzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 117);

        auto g_z_0_xzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 118);

        auto g_z_0_xzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 119);

        auto g_z_0_yyy_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 126);

        auto g_z_0_yyy_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 127);

        auto g_z_0_yyy_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 128);

        auto g_z_0_yyy_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 129);

        auto g_z_0_yyy_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 130);

        auto g_z_0_yyy_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 131);

        auto g_z_0_yyy_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 132);

        auto g_z_0_yyy_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 133);

        auto g_z_0_yyy_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 134);

        auto g_z_0_yyy_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 135);

        auto g_z_0_yyy_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 136);

        auto g_z_0_yyy_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 137);

        auto g_z_0_yyy_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 138);

        auto g_z_0_yyy_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 139);

        auto g_z_0_yyy_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 140);

        auto g_z_0_yyy_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 141);

        auto g_z_0_yyy_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 142);

        auto g_z_0_yyy_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 143);

        auto g_z_0_yyy_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 144);

        auto g_z_0_yyy_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 145);

        auto g_z_0_yyz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 147);

        auto g_z_0_yyz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 148);

        auto g_z_0_yyz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 149);

        auto g_z_0_yyz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 150);

        auto g_z_0_yyz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 151);

        auto g_z_0_yyz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 152);

        auto g_z_0_yyz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 153);

        auto g_z_0_yyz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 154);

        auto g_z_0_yyz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 155);

        auto g_z_0_yyz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 156);

        auto g_z_0_yyz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 157);

        auto g_z_0_yyz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 158);

        auto g_z_0_yyz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 159);

        auto g_z_0_yyz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 160);

        auto g_z_0_yyz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 161);

        auto g_z_0_yyz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 162);

        auto g_z_0_yyz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 163);

        auto g_z_0_yyz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 164);

        auto g_z_0_yyz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 165);

        auto g_z_0_yyz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 166);

        auto g_z_0_yzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 168);

        auto g_z_0_yzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 169);

        auto g_z_0_yzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 170);

        auto g_z_0_yzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 171);

        auto g_z_0_yzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 172);

        auto g_z_0_yzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 173);

        auto g_z_0_yzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 174);

        auto g_z_0_yzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 175);

        auto g_z_0_yzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 176);

        auto g_z_0_yzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 177);

        auto g_z_0_yzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 178);

        auto g_z_0_yzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 179);

        auto g_z_0_yzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 180);

        auto g_z_0_yzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 181);

        auto g_z_0_yzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 182);

        auto g_z_0_yzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 183);

        auto g_z_0_yzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 184);

        auto g_z_0_yzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 185);

        auto g_z_0_yzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 186);

        auto g_z_0_yzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 187);

        auto g_z_0_zzz_xxxxx = cbuffer.data(fh_geom_10_off + 420 * acomps + 189);

        auto g_z_0_zzz_xxxxy = cbuffer.data(fh_geom_10_off + 420 * acomps + 190);

        auto g_z_0_zzz_xxxxz = cbuffer.data(fh_geom_10_off + 420 * acomps + 191);

        auto g_z_0_zzz_xxxyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 192);

        auto g_z_0_zzz_xxxyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 193);

        auto g_z_0_zzz_xxxzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 194);

        auto g_z_0_zzz_xxyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 195);

        auto g_z_0_zzz_xxyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 196);

        auto g_z_0_zzz_xxyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 197);

        auto g_z_0_zzz_xxzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 198);

        auto g_z_0_zzz_xyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 199);

        auto g_z_0_zzz_xyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 200);

        auto g_z_0_zzz_xyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 201);

        auto g_z_0_zzz_xyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 202);

        auto g_z_0_zzz_xzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 203);

        auto g_z_0_zzz_yyyyy = cbuffer.data(fh_geom_10_off + 420 * acomps + 204);

        auto g_z_0_zzz_yyyyz = cbuffer.data(fh_geom_10_off + 420 * acomps + 205);

        auto g_z_0_zzz_yyyzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 206);

        auto g_z_0_zzz_yyzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 207);

        auto g_z_0_zzz_yzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 208);

        auto g_z_0_zzz_zzzzz = cbuffer.data(fh_geom_10_off + 420 * acomps + 209);

        /// set up bra offset for contr_buffer_xxgg

        const auto gg_geom_10_off = idx_geom_10_xgg + i * 225;

        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 5);

        auto g_x_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 11);

        auto g_x_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_x_0_xxx_xxxx, g_x_0_xxx_xxxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyz, g_x_0_xxx_yyzz, g_x_0_xxx_yzzz, g_x_0_xxx_zzzz, g_x_0_xxxx_xxxx, g_x_0_xxxx_xxxy, g_x_0_xxxx_xxxz, g_x_0_xxxx_xxyy, g_x_0_xxxx_xxyz, g_x_0_xxxx_xxzz, g_x_0_xxxx_xyyy, g_x_0_xxxx_xyyz, g_x_0_xxxx_xyzz, g_x_0_xxxx_xzzz, g_x_0_xxxx_yyyy, g_x_0_xxxx_yyyz, g_x_0_xxxx_yyzz, g_x_0_xxxx_yzzz, g_x_0_xxxx_zzzz, g_xxx_xxxx, g_xxx_xxxy, g_xxx_xxxz, g_xxx_xxyy, g_xxx_xxyz, g_xxx_xxzz, g_xxx_xyyy, g_xxx_xyyz, g_xxx_xyzz, g_xxx_xzzz, g_xxx_yyyy, g_xxx_yyyz, g_xxx_yyzz, g_xxx_yzzz, g_xxx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxx_xxxx[k] = -g_xxx_xxxx[k] - g_x_0_xxx_xxxx[k] * cd_x[k] + g_x_0_xxx_xxxxx[k];

            g_x_0_xxxx_xxxy[k] = -g_xxx_xxxy[k] - g_x_0_xxx_xxxy[k] * cd_x[k] + g_x_0_xxx_xxxxy[k];

            g_x_0_xxxx_xxxz[k] = -g_xxx_xxxz[k] - g_x_0_xxx_xxxz[k] * cd_x[k] + g_x_0_xxx_xxxxz[k];

            g_x_0_xxxx_xxyy[k] = -g_xxx_xxyy[k] - g_x_0_xxx_xxyy[k] * cd_x[k] + g_x_0_xxx_xxxyy[k];

            g_x_0_xxxx_xxyz[k] = -g_xxx_xxyz[k] - g_x_0_xxx_xxyz[k] * cd_x[k] + g_x_0_xxx_xxxyz[k];

            g_x_0_xxxx_xxzz[k] = -g_xxx_xxzz[k] - g_x_0_xxx_xxzz[k] * cd_x[k] + g_x_0_xxx_xxxzz[k];

            g_x_0_xxxx_xyyy[k] = -g_xxx_xyyy[k] - g_x_0_xxx_xyyy[k] * cd_x[k] + g_x_0_xxx_xxyyy[k];

            g_x_0_xxxx_xyyz[k] = -g_xxx_xyyz[k] - g_x_0_xxx_xyyz[k] * cd_x[k] + g_x_0_xxx_xxyyz[k];

            g_x_0_xxxx_xyzz[k] = -g_xxx_xyzz[k] - g_x_0_xxx_xyzz[k] * cd_x[k] + g_x_0_xxx_xxyzz[k];

            g_x_0_xxxx_xzzz[k] = -g_xxx_xzzz[k] - g_x_0_xxx_xzzz[k] * cd_x[k] + g_x_0_xxx_xxzzz[k];

            g_x_0_xxxx_yyyy[k] = -g_xxx_yyyy[k] - g_x_0_xxx_yyyy[k] * cd_x[k] + g_x_0_xxx_xyyyy[k];

            g_x_0_xxxx_yyyz[k] = -g_xxx_yyyz[k] - g_x_0_xxx_yyyz[k] * cd_x[k] + g_x_0_xxx_xyyyz[k];

            g_x_0_xxxx_yyzz[k] = -g_xxx_yyzz[k] - g_x_0_xxx_yyzz[k] * cd_x[k] + g_x_0_xxx_xyyzz[k];

            g_x_0_xxxx_yzzz[k] = -g_xxx_yzzz[k] - g_x_0_xxx_yzzz[k] * cd_x[k] + g_x_0_xxx_xyzzz[k];

            g_x_0_xxxx_zzzz[k] = -g_xxx_zzzz[k] - g_x_0_xxx_zzzz[k] * cd_x[k] + g_x_0_xxx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 17);

        auto g_x_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 23);

        auto g_x_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_xxx_xxxx, g_x_0_xxx_xxxxy, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzz, g_x_0_xxxy_xxxx, g_x_0_xxxy_xxxy, g_x_0_xxxy_xxxz, g_x_0_xxxy_xxyy, g_x_0_xxxy_xxyz, g_x_0_xxxy_xxzz, g_x_0_xxxy_xyyy, g_x_0_xxxy_xyyz, g_x_0_xxxy_xyzz, g_x_0_xxxy_xzzz, g_x_0_xxxy_yyyy, g_x_0_xxxy_yyyz, g_x_0_xxxy_yyzz, g_x_0_xxxy_yzzz, g_x_0_xxxy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxy_xxxx[k] = -g_x_0_xxx_xxxx[k] * cd_y[k] + g_x_0_xxx_xxxxy[k];

            g_x_0_xxxy_xxxy[k] = -g_x_0_xxx_xxxy[k] * cd_y[k] + g_x_0_xxx_xxxyy[k];

            g_x_0_xxxy_xxxz[k] = -g_x_0_xxx_xxxz[k] * cd_y[k] + g_x_0_xxx_xxxyz[k];

            g_x_0_xxxy_xxyy[k] = -g_x_0_xxx_xxyy[k] * cd_y[k] + g_x_0_xxx_xxyyy[k];

            g_x_0_xxxy_xxyz[k] = -g_x_0_xxx_xxyz[k] * cd_y[k] + g_x_0_xxx_xxyyz[k];

            g_x_0_xxxy_xxzz[k] = -g_x_0_xxx_xxzz[k] * cd_y[k] + g_x_0_xxx_xxyzz[k];

            g_x_0_xxxy_xyyy[k] = -g_x_0_xxx_xyyy[k] * cd_y[k] + g_x_0_xxx_xyyyy[k];

            g_x_0_xxxy_xyyz[k] = -g_x_0_xxx_xyyz[k] * cd_y[k] + g_x_0_xxx_xyyyz[k];

            g_x_0_xxxy_xyzz[k] = -g_x_0_xxx_xyzz[k] * cd_y[k] + g_x_0_xxx_xyyzz[k];

            g_x_0_xxxy_xzzz[k] = -g_x_0_xxx_xzzz[k] * cd_y[k] + g_x_0_xxx_xyzzz[k];

            g_x_0_xxxy_yyyy[k] = -g_x_0_xxx_yyyy[k] * cd_y[k] + g_x_0_xxx_yyyyy[k];

            g_x_0_xxxy_yyyz[k] = -g_x_0_xxx_yyyz[k] * cd_y[k] + g_x_0_xxx_yyyyz[k];

            g_x_0_xxxy_yyzz[k] = -g_x_0_xxx_yyzz[k] * cd_y[k] + g_x_0_xxx_yyyzz[k];

            g_x_0_xxxy_yzzz[k] = -g_x_0_xxx_yzzz[k] * cd_y[k] + g_x_0_xxx_yyzzz[k];

            g_x_0_xxxy_zzzz[k] = -g_x_0_xxx_zzzz[k] * cd_y[k] + g_x_0_xxx_yzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 35);

        auto g_x_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 41);

        auto g_x_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 44);

        #pragma omp simd aligned(cd_z, g_x_0_xxx_xxxx, g_x_0_xxx_xxxxz, g_x_0_xxx_xxxy, g_x_0_xxx_xxxyz, g_x_0_xxx_xxxz, g_x_0_xxx_xxxzz, g_x_0_xxx_xxyy, g_x_0_xxx_xxyyz, g_x_0_xxx_xxyz, g_x_0_xxx_xxyzz, g_x_0_xxx_xxzz, g_x_0_xxx_xxzzz, g_x_0_xxx_xyyy, g_x_0_xxx_xyyyz, g_x_0_xxx_xyyz, g_x_0_xxx_xyyzz, g_x_0_xxx_xyzz, g_x_0_xxx_xyzzz, g_x_0_xxx_xzzz, g_x_0_xxx_xzzzz, g_x_0_xxx_yyyy, g_x_0_xxx_yyyyz, g_x_0_xxx_yyyz, g_x_0_xxx_yyyzz, g_x_0_xxx_yyzz, g_x_0_xxx_yyzzz, g_x_0_xxx_yzzz, g_x_0_xxx_yzzzz, g_x_0_xxx_zzzz, g_x_0_xxx_zzzzz, g_x_0_xxxz_xxxx, g_x_0_xxxz_xxxy, g_x_0_xxxz_xxxz, g_x_0_xxxz_xxyy, g_x_0_xxxz_xxyz, g_x_0_xxxz_xxzz, g_x_0_xxxz_xyyy, g_x_0_xxxz_xyyz, g_x_0_xxxz_xyzz, g_x_0_xxxz_xzzz, g_x_0_xxxz_yyyy, g_x_0_xxxz_yyyz, g_x_0_xxxz_yyzz, g_x_0_xxxz_yzzz, g_x_0_xxxz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxxz_xxxx[k] = -g_x_0_xxx_xxxx[k] * cd_z[k] + g_x_0_xxx_xxxxz[k];

            g_x_0_xxxz_xxxy[k] = -g_x_0_xxx_xxxy[k] * cd_z[k] + g_x_0_xxx_xxxyz[k];

            g_x_0_xxxz_xxxz[k] = -g_x_0_xxx_xxxz[k] * cd_z[k] + g_x_0_xxx_xxxzz[k];

            g_x_0_xxxz_xxyy[k] = -g_x_0_xxx_xxyy[k] * cd_z[k] + g_x_0_xxx_xxyyz[k];

            g_x_0_xxxz_xxyz[k] = -g_x_0_xxx_xxyz[k] * cd_z[k] + g_x_0_xxx_xxyzz[k];

            g_x_0_xxxz_xxzz[k] = -g_x_0_xxx_xxzz[k] * cd_z[k] + g_x_0_xxx_xxzzz[k];

            g_x_0_xxxz_xyyy[k] = -g_x_0_xxx_xyyy[k] * cd_z[k] + g_x_0_xxx_xyyyz[k];

            g_x_0_xxxz_xyyz[k] = -g_x_0_xxx_xyyz[k] * cd_z[k] + g_x_0_xxx_xyyzz[k];

            g_x_0_xxxz_xyzz[k] = -g_x_0_xxx_xyzz[k] * cd_z[k] + g_x_0_xxx_xyzzz[k];

            g_x_0_xxxz_xzzz[k] = -g_x_0_xxx_xzzz[k] * cd_z[k] + g_x_0_xxx_xzzzz[k];

            g_x_0_xxxz_yyyy[k] = -g_x_0_xxx_yyyy[k] * cd_z[k] + g_x_0_xxx_yyyyz[k];

            g_x_0_xxxz_yyyz[k] = -g_x_0_xxx_yyyz[k] * cd_z[k] + g_x_0_xxx_yyyzz[k];

            g_x_0_xxxz_yyzz[k] = -g_x_0_xxx_yyzz[k] * cd_z[k] + g_x_0_xxx_yyzzz[k];

            g_x_0_xxxz_yzzz[k] = -g_x_0_xxx_yzzz[k] * cd_z[k] + g_x_0_xxx_yzzzz[k];

            g_x_0_xxxz_zzzz[k] = -g_x_0_xxx_zzzz[k] * cd_z[k] + g_x_0_xxx_zzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 47);

        auto g_x_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 53);

        auto g_x_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_y, g_x_0_xxy_xxxx, g_x_0_xxy_xxxxy, g_x_0_xxy_xxxy, g_x_0_xxy_xxxyy, g_x_0_xxy_xxxyz, g_x_0_xxy_xxxz, g_x_0_xxy_xxyy, g_x_0_xxy_xxyyy, g_x_0_xxy_xxyyz, g_x_0_xxy_xxyz, g_x_0_xxy_xxyzz, g_x_0_xxy_xxzz, g_x_0_xxy_xyyy, g_x_0_xxy_xyyyy, g_x_0_xxy_xyyyz, g_x_0_xxy_xyyz, g_x_0_xxy_xyyzz, g_x_0_xxy_xyzz, g_x_0_xxy_xyzzz, g_x_0_xxy_xzzz, g_x_0_xxy_yyyy, g_x_0_xxy_yyyyy, g_x_0_xxy_yyyyz, g_x_0_xxy_yyyz, g_x_0_xxy_yyyzz, g_x_0_xxy_yyzz, g_x_0_xxy_yyzzz, g_x_0_xxy_yzzz, g_x_0_xxy_yzzzz, g_x_0_xxy_zzzz, g_x_0_xxyy_xxxx, g_x_0_xxyy_xxxy, g_x_0_xxyy_xxxz, g_x_0_xxyy_xxyy, g_x_0_xxyy_xxyz, g_x_0_xxyy_xxzz, g_x_0_xxyy_xyyy, g_x_0_xxyy_xyyz, g_x_0_xxyy_xyzz, g_x_0_xxyy_xzzz, g_x_0_xxyy_yyyy, g_x_0_xxyy_yyyz, g_x_0_xxyy_yyzz, g_x_0_xxyy_yzzz, g_x_0_xxyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyy_xxxx[k] = -g_x_0_xxy_xxxx[k] * cd_y[k] + g_x_0_xxy_xxxxy[k];

            g_x_0_xxyy_xxxy[k] = -g_x_0_xxy_xxxy[k] * cd_y[k] + g_x_0_xxy_xxxyy[k];

            g_x_0_xxyy_xxxz[k] = -g_x_0_xxy_xxxz[k] * cd_y[k] + g_x_0_xxy_xxxyz[k];

            g_x_0_xxyy_xxyy[k] = -g_x_0_xxy_xxyy[k] * cd_y[k] + g_x_0_xxy_xxyyy[k];

            g_x_0_xxyy_xxyz[k] = -g_x_0_xxy_xxyz[k] * cd_y[k] + g_x_0_xxy_xxyyz[k];

            g_x_0_xxyy_xxzz[k] = -g_x_0_xxy_xxzz[k] * cd_y[k] + g_x_0_xxy_xxyzz[k];

            g_x_0_xxyy_xyyy[k] = -g_x_0_xxy_xyyy[k] * cd_y[k] + g_x_0_xxy_xyyyy[k];

            g_x_0_xxyy_xyyz[k] = -g_x_0_xxy_xyyz[k] * cd_y[k] + g_x_0_xxy_xyyyz[k];

            g_x_0_xxyy_xyzz[k] = -g_x_0_xxy_xyzz[k] * cd_y[k] + g_x_0_xxy_xyyzz[k];

            g_x_0_xxyy_xzzz[k] = -g_x_0_xxy_xzzz[k] * cd_y[k] + g_x_0_xxy_xyzzz[k];

            g_x_0_xxyy_yyyy[k] = -g_x_0_xxy_yyyy[k] * cd_y[k] + g_x_0_xxy_yyyyy[k];

            g_x_0_xxyy_yyyz[k] = -g_x_0_xxy_yyyz[k] * cd_y[k] + g_x_0_xxy_yyyyz[k];

            g_x_0_xxyy_yyzz[k] = -g_x_0_xxy_yyzz[k] * cd_y[k] + g_x_0_xxy_yyyzz[k];

            g_x_0_xxyy_yzzz[k] = -g_x_0_xxy_yzzz[k] * cd_y[k] + g_x_0_xxy_yyzzz[k];

            g_x_0_xxyy_zzzz[k] = -g_x_0_xxy_zzzz[k] * cd_y[k] + g_x_0_xxy_yzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 60);

        auto g_x_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 61);

        auto g_x_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 62);

        auto g_x_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 63);

        auto g_x_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 64);

        auto g_x_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 65);

        auto g_x_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 66);

        auto g_x_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 67);

        auto g_x_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 68);

        auto g_x_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 69);

        auto g_x_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 70);

        auto g_x_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 71);

        auto g_x_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 72);

        auto g_x_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 73);

        auto g_x_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 74);

        #pragma omp simd aligned(cd_y, g_x_0_xxyz_xxxx, g_x_0_xxyz_xxxy, g_x_0_xxyz_xxxz, g_x_0_xxyz_xxyy, g_x_0_xxyz_xxyz, g_x_0_xxyz_xxzz, g_x_0_xxyz_xyyy, g_x_0_xxyz_xyyz, g_x_0_xxyz_xyzz, g_x_0_xxyz_xzzz, g_x_0_xxyz_yyyy, g_x_0_xxyz_yyyz, g_x_0_xxyz_yyzz, g_x_0_xxyz_yzzz, g_x_0_xxyz_zzzz, g_x_0_xxz_xxxx, g_x_0_xxz_xxxxy, g_x_0_xxz_xxxy, g_x_0_xxz_xxxyy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxz, g_x_0_xxz_xxyy, g_x_0_xxz_xxyyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzz, g_x_0_xxz_xyyy, g_x_0_xxz_xyyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzz, g_x_0_xxz_yyyy, g_x_0_xxz_yyyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxyz_xxxx[k] = -g_x_0_xxz_xxxx[k] * cd_y[k] + g_x_0_xxz_xxxxy[k];

            g_x_0_xxyz_xxxy[k] = -g_x_0_xxz_xxxy[k] * cd_y[k] + g_x_0_xxz_xxxyy[k];

            g_x_0_xxyz_xxxz[k] = -g_x_0_xxz_xxxz[k] * cd_y[k] + g_x_0_xxz_xxxyz[k];

            g_x_0_xxyz_xxyy[k] = -g_x_0_xxz_xxyy[k] * cd_y[k] + g_x_0_xxz_xxyyy[k];

            g_x_0_xxyz_xxyz[k] = -g_x_0_xxz_xxyz[k] * cd_y[k] + g_x_0_xxz_xxyyz[k];

            g_x_0_xxyz_xxzz[k] = -g_x_0_xxz_xxzz[k] * cd_y[k] + g_x_0_xxz_xxyzz[k];

            g_x_0_xxyz_xyyy[k] = -g_x_0_xxz_xyyy[k] * cd_y[k] + g_x_0_xxz_xyyyy[k];

            g_x_0_xxyz_xyyz[k] = -g_x_0_xxz_xyyz[k] * cd_y[k] + g_x_0_xxz_xyyyz[k];

            g_x_0_xxyz_xyzz[k] = -g_x_0_xxz_xyzz[k] * cd_y[k] + g_x_0_xxz_xyyzz[k];

            g_x_0_xxyz_xzzz[k] = -g_x_0_xxz_xzzz[k] * cd_y[k] + g_x_0_xxz_xyzzz[k];

            g_x_0_xxyz_yyyy[k] = -g_x_0_xxz_yyyy[k] * cd_y[k] + g_x_0_xxz_yyyyy[k];

            g_x_0_xxyz_yyyz[k] = -g_x_0_xxz_yyyz[k] * cd_y[k] + g_x_0_xxz_yyyyz[k];

            g_x_0_xxyz_yyzz[k] = -g_x_0_xxz_yyzz[k] * cd_y[k] + g_x_0_xxz_yyyzz[k];

            g_x_0_xxyz_yzzz[k] = -g_x_0_xxz_yzzz[k] * cd_y[k] + g_x_0_xxz_yyzzz[k];

            g_x_0_xxyz_zzzz[k] = -g_x_0_xxz_zzzz[k] * cd_y[k] + g_x_0_xxz_yzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 75);

        auto g_x_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 76);

        auto g_x_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 77);

        auto g_x_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 78);

        auto g_x_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 79);

        auto g_x_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 80);

        auto g_x_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 81);

        auto g_x_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 82);

        auto g_x_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 83);

        auto g_x_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 84);

        auto g_x_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 85);

        auto g_x_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 86);

        auto g_x_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 87);

        auto g_x_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 88);

        auto g_x_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 89);

        #pragma omp simd aligned(cd_z, g_x_0_xxz_xxxx, g_x_0_xxz_xxxxz, g_x_0_xxz_xxxy, g_x_0_xxz_xxxyz, g_x_0_xxz_xxxz, g_x_0_xxz_xxxzz, g_x_0_xxz_xxyy, g_x_0_xxz_xxyyz, g_x_0_xxz_xxyz, g_x_0_xxz_xxyzz, g_x_0_xxz_xxzz, g_x_0_xxz_xxzzz, g_x_0_xxz_xyyy, g_x_0_xxz_xyyyz, g_x_0_xxz_xyyz, g_x_0_xxz_xyyzz, g_x_0_xxz_xyzz, g_x_0_xxz_xyzzz, g_x_0_xxz_xzzz, g_x_0_xxz_xzzzz, g_x_0_xxz_yyyy, g_x_0_xxz_yyyyz, g_x_0_xxz_yyyz, g_x_0_xxz_yyyzz, g_x_0_xxz_yyzz, g_x_0_xxz_yyzzz, g_x_0_xxz_yzzz, g_x_0_xxz_yzzzz, g_x_0_xxz_zzzz, g_x_0_xxz_zzzzz, g_x_0_xxzz_xxxx, g_x_0_xxzz_xxxy, g_x_0_xxzz_xxxz, g_x_0_xxzz_xxyy, g_x_0_xxzz_xxyz, g_x_0_xxzz_xxzz, g_x_0_xxzz_xyyy, g_x_0_xxzz_xyyz, g_x_0_xxzz_xyzz, g_x_0_xxzz_xzzz, g_x_0_xxzz_yyyy, g_x_0_xxzz_yyyz, g_x_0_xxzz_yyzz, g_x_0_xxzz_yzzz, g_x_0_xxzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxzz_xxxx[k] = -g_x_0_xxz_xxxx[k] * cd_z[k] + g_x_0_xxz_xxxxz[k];

            g_x_0_xxzz_xxxy[k] = -g_x_0_xxz_xxxy[k] * cd_z[k] + g_x_0_xxz_xxxyz[k];

            g_x_0_xxzz_xxxz[k] = -g_x_0_xxz_xxxz[k] * cd_z[k] + g_x_0_xxz_xxxzz[k];

            g_x_0_xxzz_xxyy[k] = -g_x_0_xxz_xxyy[k] * cd_z[k] + g_x_0_xxz_xxyyz[k];

            g_x_0_xxzz_xxyz[k] = -g_x_0_xxz_xxyz[k] * cd_z[k] + g_x_0_xxz_xxyzz[k];

            g_x_0_xxzz_xxzz[k] = -g_x_0_xxz_xxzz[k] * cd_z[k] + g_x_0_xxz_xxzzz[k];

            g_x_0_xxzz_xyyy[k] = -g_x_0_xxz_xyyy[k] * cd_z[k] + g_x_0_xxz_xyyyz[k];

            g_x_0_xxzz_xyyz[k] = -g_x_0_xxz_xyyz[k] * cd_z[k] + g_x_0_xxz_xyyzz[k];

            g_x_0_xxzz_xyzz[k] = -g_x_0_xxz_xyzz[k] * cd_z[k] + g_x_0_xxz_xyzzz[k];

            g_x_0_xxzz_xzzz[k] = -g_x_0_xxz_xzzz[k] * cd_z[k] + g_x_0_xxz_xzzzz[k];

            g_x_0_xxzz_yyyy[k] = -g_x_0_xxz_yyyy[k] * cd_z[k] + g_x_0_xxz_yyyyz[k];

            g_x_0_xxzz_yyyz[k] = -g_x_0_xxz_yyyz[k] * cd_z[k] + g_x_0_xxz_yyyzz[k];

            g_x_0_xxzz_yyzz[k] = -g_x_0_xxz_yyzz[k] * cd_z[k] + g_x_0_xxz_yyzzz[k];

            g_x_0_xxzz_yzzz[k] = -g_x_0_xxz_yzzz[k] * cd_z[k] + g_x_0_xxz_yzzzz[k];

            g_x_0_xxzz_zzzz[k] = -g_x_0_xxz_zzzz[k] * cd_z[k] + g_x_0_xxz_zzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 90);

        auto g_x_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 91);

        auto g_x_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 92);

        auto g_x_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 93);

        auto g_x_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 94);

        auto g_x_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 95);

        auto g_x_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 96);

        auto g_x_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 97);

        auto g_x_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 98);

        auto g_x_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 99);

        auto g_x_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 100);

        auto g_x_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 101);

        auto g_x_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 102);

        auto g_x_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 103);

        auto g_x_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 104);

        #pragma omp simd aligned(cd_y, g_x_0_xyy_xxxx, g_x_0_xyy_xxxxy, g_x_0_xyy_xxxy, g_x_0_xyy_xxxyy, g_x_0_xyy_xxxyz, g_x_0_xyy_xxxz, g_x_0_xyy_xxyy, g_x_0_xyy_xxyyy, g_x_0_xyy_xxyyz, g_x_0_xyy_xxyz, g_x_0_xyy_xxyzz, g_x_0_xyy_xxzz, g_x_0_xyy_xyyy, g_x_0_xyy_xyyyy, g_x_0_xyy_xyyyz, g_x_0_xyy_xyyz, g_x_0_xyy_xyyzz, g_x_0_xyy_xyzz, g_x_0_xyy_xyzzz, g_x_0_xyy_xzzz, g_x_0_xyy_yyyy, g_x_0_xyy_yyyyy, g_x_0_xyy_yyyyz, g_x_0_xyy_yyyz, g_x_0_xyy_yyyzz, g_x_0_xyy_yyzz, g_x_0_xyy_yyzzz, g_x_0_xyy_yzzz, g_x_0_xyy_yzzzz, g_x_0_xyy_zzzz, g_x_0_xyyy_xxxx, g_x_0_xyyy_xxxy, g_x_0_xyyy_xxxz, g_x_0_xyyy_xxyy, g_x_0_xyyy_xxyz, g_x_0_xyyy_xxzz, g_x_0_xyyy_xyyy, g_x_0_xyyy_xyyz, g_x_0_xyyy_xyzz, g_x_0_xyyy_xzzz, g_x_0_xyyy_yyyy, g_x_0_xyyy_yyyz, g_x_0_xyyy_yyzz, g_x_0_xyyy_yzzz, g_x_0_xyyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyy_xxxx[k] = -g_x_0_xyy_xxxx[k] * cd_y[k] + g_x_0_xyy_xxxxy[k];

            g_x_0_xyyy_xxxy[k] = -g_x_0_xyy_xxxy[k] * cd_y[k] + g_x_0_xyy_xxxyy[k];

            g_x_0_xyyy_xxxz[k] = -g_x_0_xyy_xxxz[k] * cd_y[k] + g_x_0_xyy_xxxyz[k];

            g_x_0_xyyy_xxyy[k] = -g_x_0_xyy_xxyy[k] * cd_y[k] + g_x_0_xyy_xxyyy[k];

            g_x_0_xyyy_xxyz[k] = -g_x_0_xyy_xxyz[k] * cd_y[k] + g_x_0_xyy_xxyyz[k];

            g_x_0_xyyy_xxzz[k] = -g_x_0_xyy_xxzz[k] * cd_y[k] + g_x_0_xyy_xxyzz[k];

            g_x_0_xyyy_xyyy[k] = -g_x_0_xyy_xyyy[k] * cd_y[k] + g_x_0_xyy_xyyyy[k];

            g_x_0_xyyy_xyyz[k] = -g_x_0_xyy_xyyz[k] * cd_y[k] + g_x_0_xyy_xyyyz[k];

            g_x_0_xyyy_xyzz[k] = -g_x_0_xyy_xyzz[k] * cd_y[k] + g_x_0_xyy_xyyzz[k];

            g_x_0_xyyy_xzzz[k] = -g_x_0_xyy_xzzz[k] * cd_y[k] + g_x_0_xyy_xyzzz[k];

            g_x_0_xyyy_yyyy[k] = -g_x_0_xyy_yyyy[k] * cd_y[k] + g_x_0_xyy_yyyyy[k];

            g_x_0_xyyy_yyyz[k] = -g_x_0_xyy_yyyz[k] * cd_y[k] + g_x_0_xyy_yyyyz[k];

            g_x_0_xyyy_yyzz[k] = -g_x_0_xyy_yyzz[k] * cd_y[k] + g_x_0_xyy_yyyzz[k];

            g_x_0_xyyy_yzzz[k] = -g_x_0_xyy_yzzz[k] * cd_y[k] + g_x_0_xyy_yyzzz[k];

            g_x_0_xyyy_zzzz[k] = -g_x_0_xyy_zzzz[k] * cd_y[k] + g_x_0_xyy_yzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 105);

        auto g_x_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 106);

        auto g_x_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 107);

        auto g_x_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 108);

        auto g_x_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 109);

        auto g_x_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 110);

        auto g_x_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 111);

        auto g_x_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 112);

        auto g_x_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 113);

        auto g_x_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 114);

        auto g_x_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 115);

        auto g_x_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 116);

        auto g_x_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 117);

        auto g_x_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 118);

        auto g_x_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 119);

        #pragma omp simd aligned(cd_y, g_x_0_xyyz_xxxx, g_x_0_xyyz_xxxy, g_x_0_xyyz_xxxz, g_x_0_xyyz_xxyy, g_x_0_xyyz_xxyz, g_x_0_xyyz_xxzz, g_x_0_xyyz_xyyy, g_x_0_xyyz_xyyz, g_x_0_xyyz_xyzz, g_x_0_xyyz_xzzz, g_x_0_xyyz_yyyy, g_x_0_xyyz_yyyz, g_x_0_xyyz_yyzz, g_x_0_xyyz_yzzz, g_x_0_xyyz_zzzz, g_x_0_xyz_xxxx, g_x_0_xyz_xxxxy, g_x_0_xyz_xxxy, g_x_0_xyz_xxxyy, g_x_0_xyz_xxxyz, g_x_0_xyz_xxxz, g_x_0_xyz_xxyy, g_x_0_xyz_xxyyy, g_x_0_xyz_xxyyz, g_x_0_xyz_xxyz, g_x_0_xyz_xxyzz, g_x_0_xyz_xxzz, g_x_0_xyz_xyyy, g_x_0_xyz_xyyyy, g_x_0_xyz_xyyyz, g_x_0_xyz_xyyz, g_x_0_xyz_xyyzz, g_x_0_xyz_xyzz, g_x_0_xyz_xyzzz, g_x_0_xyz_xzzz, g_x_0_xyz_yyyy, g_x_0_xyz_yyyyy, g_x_0_xyz_yyyyz, g_x_0_xyz_yyyz, g_x_0_xyz_yyyzz, g_x_0_xyz_yyzz, g_x_0_xyz_yyzzz, g_x_0_xyz_yzzz, g_x_0_xyz_yzzzz, g_x_0_xyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyyz_xxxx[k] = -g_x_0_xyz_xxxx[k] * cd_y[k] + g_x_0_xyz_xxxxy[k];

            g_x_0_xyyz_xxxy[k] = -g_x_0_xyz_xxxy[k] * cd_y[k] + g_x_0_xyz_xxxyy[k];

            g_x_0_xyyz_xxxz[k] = -g_x_0_xyz_xxxz[k] * cd_y[k] + g_x_0_xyz_xxxyz[k];

            g_x_0_xyyz_xxyy[k] = -g_x_0_xyz_xxyy[k] * cd_y[k] + g_x_0_xyz_xxyyy[k];

            g_x_0_xyyz_xxyz[k] = -g_x_0_xyz_xxyz[k] * cd_y[k] + g_x_0_xyz_xxyyz[k];

            g_x_0_xyyz_xxzz[k] = -g_x_0_xyz_xxzz[k] * cd_y[k] + g_x_0_xyz_xxyzz[k];

            g_x_0_xyyz_xyyy[k] = -g_x_0_xyz_xyyy[k] * cd_y[k] + g_x_0_xyz_xyyyy[k];

            g_x_0_xyyz_xyyz[k] = -g_x_0_xyz_xyyz[k] * cd_y[k] + g_x_0_xyz_xyyyz[k];

            g_x_0_xyyz_xyzz[k] = -g_x_0_xyz_xyzz[k] * cd_y[k] + g_x_0_xyz_xyyzz[k];

            g_x_0_xyyz_xzzz[k] = -g_x_0_xyz_xzzz[k] * cd_y[k] + g_x_0_xyz_xyzzz[k];

            g_x_0_xyyz_yyyy[k] = -g_x_0_xyz_yyyy[k] * cd_y[k] + g_x_0_xyz_yyyyy[k];

            g_x_0_xyyz_yyyz[k] = -g_x_0_xyz_yyyz[k] * cd_y[k] + g_x_0_xyz_yyyyz[k];

            g_x_0_xyyz_yyzz[k] = -g_x_0_xyz_yyzz[k] * cd_y[k] + g_x_0_xyz_yyyzz[k];

            g_x_0_xyyz_yzzz[k] = -g_x_0_xyz_yzzz[k] * cd_y[k] + g_x_0_xyz_yyzzz[k];

            g_x_0_xyyz_zzzz[k] = -g_x_0_xyz_zzzz[k] * cd_y[k] + g_x_0_xyz_yzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 120);

        auto g_x_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 121);

        auto g_x_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 122);

        auto g_x_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 123);

        auto g_x_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 124);

        auto g_x_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 125);

        auto g_x_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 126);

        auto g_x_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 127);

        auto g_x_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 128);

        auto g_x_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 129);

        auto g_x_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 130);

        auto g_x_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 131);

        auto g_x_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 132);

        auto g_x_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 133);

        auto g_x_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 134);

        #pragma omp simd aligned(cd_y, g_x_0_xyzz_xxxx, g_x_0_xyzz_xxxy, g_x_0_xyzz_xxxz, g_x_0_xyzz_xxyy, g_x_0_xyzz_xxyz, g_x_0_xyzz_xxzz, g_x_0_xyzz_xyyy, g_x_0_xyzz_xyyz, g_x_0_xyzz_xyzz, g_x_0_xyzz_xzzz, g_x_0_xyzz_yyyy, g_x_0_xyzz_yyyz, g_x_0_xyzz_yyzz, g_x_0_xyzz_yzzz, g_x_0_xyzz_zzzz, g_x_0_xzz_xxxx, g_x_0_xzz_xxxxy, g_x_0_xzz_xxxy, g_x_0_xzz_xxxyy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxz, g_x_0_xzz_xxyy, g_x_0_xzz_xxyyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzz, g_x_0_xzz_xyyy, g_x_0_xzz_xyyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzz, g_x_0_xzz_yyyy, g_x_0_xzz_yyyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyzz_xxxx[k] = -g_x_0_xzz_xxxx[k] * cd_y[k] + g_x_0_xzz_xxxxy[k];

            g_x_0_xyzz_xxxy[k] = -g_x_0_xzz_xxxy[k] * cd_y[k] + g_x_0_xzz_xxxyy[k];

            g_x_0_xyzz_xxxz[k] = -g_x_0_xzz_xxxz[k] * cd_y[k] + g_x_0_xzz_xxxyz[k];

            g_x_0_xyzz_xxyy[k] = -g_x_0_xzz_xxyy[k] * cd_y[k] + g_x_0_xzz_xxyyy[k];

            g_x_0_xyzz_xxyz[k] = -g_x_0_xzz_xxyz[k] * cd_y[k] + g_x_0_xzz_xxyyz[k];

            g_x_0_xyzz_xxzz[k] = -g_x_0_xzz_xxzz[k] * cd_y[k] + g_x_0_xzz_xxyzz[k];

            g_x_0_xyzz_xyyy[k] = -g_x_0_xzz_xyyy[k] * cd_y[k] + g_x_0_xzz_xyyyy[k];

            g_x_0_xyzz_xyyz[k] = -g_x_0_xzz_xyyz[k] * cd_y[k] + g_x_0_xzz_xyyyz[k];

            g_x_0_xyzz_xyzz[k] = -g_x_0_xzz_xyzz[k] * cd_y[k] + g_x_0_xzz_xyyzz[k];

            g_x_0_xyzz_xzzz[k] = -g_x_0_xzz_xzzz[k] * cd_y[k] + g_x_0_xzz_xyzzz[k];

            g_x_0_xyzz_yyyy[k] = -g_x_0_xzz_yyyy[k] * cd_y[k] + g_x_0_xzz_yyyyy[k];

            g_x_0_xyzz_yyyz[k] = -g_x_0_xzz_yyyz[k] * cd_y[k] + g_x_0_xzz_yyyyz[k];

            g_x_0_xyzz_yyzz[k] = -g_x_0_xzz_yyzz[k] * cd_y[k] + g_x_0_xzz_yyyzz[k];

            g_x_0_xyzz_yzzz[k] = -g_x_0_xzz_yzzz[k] * cd_y[k] + g_x_0_xzz_yyzzz[k];

            g_x_0_xyzz_zzzz[k] = -g_x_0_xzz_zzzz[k] * cd_y[k] + g_x_0_xzz_yzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 135);

        auto g_x_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 136);

        auto g_x_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 137);

        auto g_x_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 138);

        auto g_x_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 139);

        auto g_x_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 140);

        auto g_x_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 141);

        auto g_x_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 142);

        auto g_x_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 143);

        auto g_x_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 144);

        auto g_x_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 145);

        auto g_x_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 146);

        auto g_x_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 147);

        auto g_x_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 148);

        auto g_x_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 149);

        #pragma omp simd aligned(cd_z, g_x_0_xzz_xxxx, g_x_0_xzz_xxxxz, g_x_0_xzz_xxxy, g_x_0_xzz_xxxyz, g_x_0_xzz_xxxz, g_x_0_xzz_xxxzz, g_x_0_xzz_xxyy, g_x_0_xzz_xxyyz, g_x_0_xzz_xxyz, g_x_0_xzz_xxyzz, g_x_0_xzz_xxzz, g_x_0_xzz_xxzzz, g_x_0_xzz_xyyy, g_x_0_xzz_xyyyz, g_x_0_xzz_xyyz, g_x_0_xzz_xyyzz, g_x_0_xzz_xyzz, g_x_0_xzz_xyzzz, g_x_0_xzz_xzzz, g_x_0_xzz_xzzzz, g_x_0_xzz_yyyy, g_x_0_xzz_yyyyz, g_x_0_xzz_yyyz, g_x_0_xzz_yyyzz, g_x_0_xzz_yyzz, g_x_0_xzz_yyzzz, g_x_0_xzz_yzzz, g_x_0_xzz_yzzzz, g_x_0_xzz_zzzz, g_x_0_xzz_zzzzz, g_x_0_xzzz_xxxx, g_x_0_xzzz_xxxy, g_x_0_xzzz_xxxz, g_x_0_xzzz_xxyy, g_x_0_xzzz_xxyz, g_x_0_xzzz_xxzz, g_x_0_xzzz_xyyy, g_x_0_xzzz_xyyz, g_x_0_xzzz_xyzz, g_x_0_xzzz_xzzz, g_x_0_xzzz_yyyy, g_x_0_xzzz_yyyz, g_x_0_xzzz_yyzz, g_x_0_xzzz_yzzz, g_x_0_xzzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzzz_xxxx[k] = -g_x_0_xzz_xxxx[k] * cd_z[k] + g_x_0_xzz_xxxxz[k];

            g_x_0_xzzz_xxxy[k] = -g_x_0_xzz_xxxy[k] * cd_z[k] + g_x_0_xzz_xxxyz[k];

            g_x_0_xzzz_xxxz[k] = -g_x_0_xzz_xxxz[k] * cd_z[k] + g_x_0_xzz_xxxzz[k];

            g_x_0_xzzz_xxyy[k] = -g_x_0_xzz_xxyy[k] * cd_z[k] + g_x_0_xzz_xxyyz[k];

            g_x_0_xzzz_xxyz[k] = -g_x_0_xzz_xxyz[k] * cd_z[k] + g_x_0_xzz_xxyzz[k];

            g_x_0_xzzz_xxzz[k] = -g_x_0_xzz_xxzz[k] * cd_z[k] + g_x_0_xzz_xxzzz[k];

            g_x_0_xzzz_xyyy[k] = -g_x_0_xzz_xyyy[k] * cd_z[k] + g_x_0_xzz_xyyyz[k];

            g_x_0_xzzz_xyyz[k] = -g_x_0_xzz_xyyz[k] * cd_z[k] + g_x_0_xzz_xyyzz[k];

            g_x_0_xzzz_xyzz[k] = -g_x_0_xzz_xyzz[k] * cd_z[k] + g_x_0_xzz_xyzzz[k];

            g_x_0_xzzz_xzzz[k] = -g_x_0_xzz_xzzz[k] * cd_z[k] + g_x_0_xzz_xzzzz[k];

            g_x_0_xzzz_yyyy[k] = -g_x_0_xzz_yyyy[k] * cd_z[k] + g_x_0_xzz_yyyyz[k];

            g_x_0_xzzz_yyyz[k] = -g_x_0_xzz_yyyz[k] * cd_z[k] + g_x_0_xzz_yyyzz[k];

            g_x_0_xzzz_yyzz[k] = -g_x_0_xzz_yyzz[k] * cd_z[k] + g_x_0_xzz_yyzzz[k];

            g_x_0_xzzz_yzzz[k] = -g_x_0_xzz_yzzz[k] * cd_z[k] + g_x_0_xzz_yzzzz[k];

            g_x_0_xzzz_zzzz[k] = -g_x_0_xzz_zzzz[k] * cd_z[k] + g_x_0_xzz_zzzzz[k];
        }

        /// Set up 150-165 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 150);

        auto g_x_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 151);

        auto g_x_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 152);

        auto g_x_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 153);

        auto g_x_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 154);

        auto g_x_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 155);

        auto g_x_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 156);

        auto g_x_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 157);

        auto g_x_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 158);

        auto g_x_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 159);

        auto g_x_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 160);

        auto g_x_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 161);

        auto g_x_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 162);

        auto g_x_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 163);

        auto g_x_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 164);

        #pragma omp simd aligned(cd_y, g_x_0_yyy_xxxx, g_x_0_yyy_xxxxy, g_x_0_yyy_xxxy, g_x_0_yyy_xxxyy, g_x_0_yyy_xxxyz, g_x_0_yyy_xxxz, g_x_0_yyy_xxyy, g_x_0_yyy_xxyyy, g_x_0_yyy_xxyyz, g_x_0_yyy_xxyz, g_x_0_yyy_xxyzz, g_x_0_yyy_xxzz, g_x_0_yyy_xyyy, g_x_0_yyy_xyyyy, g_x_0_yyy_xyyyz, g_x_0_yyy_xyyz, g_x_0_yyy_xyyzz, g_x_0_yyy_xyzz, g_x_0_yyy_xyzzz, g_x_0_yyy_xzzz, g_x_0_yyy_yyyy, g_x_0_yyy_yyyyy, g_x_0_yyy_yyyyz, g_x_0_yyy_yyyz, g_x_0_yyy_yyyzz, g_x_0_yyy_yyzz, g_x_0_yyy_yyzzz, g_x_0_yyy_yzzz, g_x_0_yyy_yzzzz, g_x_0_yyy_zzzz, g_x_0_yyyy_xxxx, g_x_0_yyyy_xxxy, g_x_0_yyyy_xxxz, g_x_0_yyyy_xxyy, g_x_0_yyyy_xxyz, g_x_0_yyyy_xxzz, g_x_0_yyyy_xyyy, g_x_0_yyyy_xyyz, g_x_0_yyyy_xyzz, g_x_0_yyyy_xzzz, g_x_0_yyyy_yyyy, g_x_0_yyyy_yyyz, g_x_0_yyyy_yyzz, g_x_0_yyyy_yzzz, g_x_0_yyyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyy_xxxx[k] = -g_x_0_yyy_xxxx[k] * cd_y[k] + g_x_0_yyy_xxxxy[k];

            g_x_0_yyyy_xxxy[k] = -g_x_0_yyy_xxxy[k] * cd_y[k] + g_x_0_yyy_xxxyy[k];

            g_x_0_yyyy_xxxz[k] = -g_x_0_yyy_xxxz[k] * cd_y[k] + g_x_0_yyy_xxxyz[k];

            g_x_0_yyyy_xxyy[k] = -g_x_0_yyy_xxyy[k] * cd_y[k] + g_x_0_yyy_xxyyy[k];

            g_x_0_yyyy_xxyz[k] = -g_x_0_yyy_xxyz[k] * cd_y[k] + g_x_0_yyy_xxyyz[k];

            g_x_0_yyyy_xxzz[k] = -g_x_0_yyy_xxzz[k] * cd_y[k] + g_x_0_yyy_xxyzz[k];

            g_x_0_yyyy_xyyy[k] = -g_x_0_yyy_xyyy[k] * cd_y[k] + g_x_0_yyy_xyyyy[k];

            g_x_0_yyyy_xyyz[k] = -g_x_0_yyy_xyyz[k] * cd_y[k] + g_x_0_yyy_xyyyz[k];

            g_x_0_yyyy_xyzz[k] = -g_x_0_yyy_xyzz[k] * cd_y[k] + g_x_0_yyy_xyyzz[k];

            g_x_0_yyyy_xzzz[k] = -g_x_0_yyy_xzzz[k] * cd_y[k] + g_x_0_yyy_xyzzz[k];

            g_x_0_yyyy_yyyy[k] = -g_x_0_yyy_yyyy[k] * cd_y[k] + g_x_0_yyy_yyyyy[k];

            g_x_0_yyyy_yyyz[k] = -g_x_0_yyy_yyyz[k] * cd_y[k] + g_x_0_yyy_yyyyz[k];

            g_x_0_yyyy_yyzz[k] = -g_x_0_yyy_yyzz[k] * cd_y[k] + g_x_0_yyy_yyyzz[k];

            g_x_0_yyyy_yzzz[k] = -g_x_0_yyy_yzzz[k] * cd_y[k] + g_x_0_yyy_yyzzz[k];

            g_x_0_yyyy_zzzz[k] = -g_x_0_yyy_zzzz[k] * cd_y[k] + g_x_0_yyy_yzzzz[k];
        }

        /// Set up 165-180 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 165);

        auto g_x_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 166);

        auto g_x_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 167);

        auto g_x_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 168);

        auto g_x_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 169);

        auto g_x_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 170);

        auto g_x_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 171);

        auto g_x_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 172);

        auto g_x_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 173);

        auto g_x_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 174);

        auto g_x_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 175);

        auto g_x_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 176);

        auto g_x_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 177);

        auto g_x_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 178);

        auto g_x_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 179);

        #pragma omp simd aligned(cd_y, g_x_0_yyyz_xxxx, g_x_0_yyyz_xxxy, g_x_0_yyyz_xxxz, g_x_0_yyyz_xxyy, g_x_0_yyyz_xxyz, g_x_0_yyyz_xxzz, g_x_0_yyyz_xyyy, g_x_0_yyyz_xyyz, g_x_0_yyyz_xyzz, g_x_0_yyyz_xzzz, g_x_0_yyyz_yyyy, g_x_0_yyyz_yyyz, g_x_0_yyyz_yyzz, g_x_0_yyyz_yzzz, g_x_0_yyyz_zzzz, g_x_0_yyz_xxxx, g_x_0_yyz_xxxxy, g_x_0_yyz_xxxy, g_x_0_yyz_xxxyy, g_x_0_yyz_xxxyz, g_x_0_yyz_xxxz, g_x_0_yyz_xxyy, g_x_0_yyz_xxyyy, g_x_0_yyz_xxyyz, g_x_0_yyz_xxyz, g_x_0_yyz_xxyzz, g_x_0_yyz_xxzz, g_x_0_yyz_xyyy, g_x_0_yyz_xyyyy, g_x_0_yyz_xyyyz, g_x_0_yyz_xyyz, g_x_0_yyz_xyyzz, g_x_0_yyz_xyzz, g_x_0_yyz_xyzzz, g_x_0_yyz_xzzz, g_x_0_yyz_yyyy, g_x_0_yyz_yyyyy, g_x_0_yyz_yyyyz, g_x_0_yyz_yyyz, g_x_0_yyz_yyyzz, g_x_0_yyz_yyzz, g_x_0_yyz_yyzzz, g_x_0_yyz_yzzz, g_x_0_yyz_yzzzz, g_x_0_yyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyyz_xxxx[k] = -g_x_0_yyz_xxxx[k] * cd_y[k] + g_x_0_yyz_xxxxy[k];

            g_x_0_yyyz_xxxy[k] = -g_x_0_yyz_xxxy[k] * cd_y[k] + g_x_0_yyz_xxxyy[k];

            g_x_0_yyyz_xxxz[k] = -g_x_0_yyz_xxxz[k] * cd_y[k] + g_x_0_yyz_xxxyz[k];

            g_x_0_yyyz_xxyy[k] = -g_x_0_yyz_xxyy[k] * cd_y[k] + g_x_0_yyz_xxyyy[k];

            g_x_0_yyyz_xxyz[k] = -g_x_0_yyz_xxyz[k] * cd_y[k] + g_x_0_yyz_xxyyz[k];

            g_x_0_yyyz_xxzz[k] = -g_x_0_yyz_xxzz[k] * cd_y[k] + g_x_0_yyz_xxyzz[k];

            g_x_0_yyyz_xyyy[k] = -g_x_0_yyz_xyyy[k] * cd_y[k] + g_x_0_yyz_xyyyy[k];

            g_x_0_yyyz_xyyz[k] = -g_x_0_yyz_xyyz[k] * cd_y[k] + g_x_0_yyz_xyyyz[k];

            g_x_0_yyyz_xyzz[k] = -g_x_0_yyz_xyzz[k] * cd_y[k] + g_x_0_yyz_xyyzz[k];

            g_x_0_yyyz_xzzz[k] = -g_x_0_yyz_xzzz[k] * cd_y[k] + g_x_0_yyz_xyzzz[k];

            g_x_0_yyyz_yyyy[k] = -g_x_0_yyz_yyyy[k] * cd_y[k] + g_x_0_yyz_yyyyy[k];

            g_x_0_yyyz_yyyz[k] = -g_x_0_yyz_yyyz[k] * cd_y[k] + g_x_0_yyz_yyyyz[k];

            g_x_0_yyyz_yyzz[k] = -g_x_0_yyz_yyzz[k] * cd_y[k] + g_x_0_yyz_yyyzz[k];

            g_x_0_yyyz_yzzz[k] = -g_x_0_yyz_yzzz[k] * cd_y[k] + g_x_0_yyz_yyzzz[k];

            g_x_0_yyyz_zzzz[k] = -g_x_0_yyz_zzzz[k] * cd_y[k] + g_x_0_yyz_yzzzz[k];
        }

        /// Set up 180-195 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 180);

        auto g_x_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 181);

        auto g_x_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 182);

        auto g_x_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 183);

        auto g_x_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 184);

        auto g_x_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 185);

        auto g_x_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 186);

        auto g_x_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 187);

        auto g_x_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 188);

        auto g_x_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 189);

        auto g_x_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 190);

        auto g_x_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 191);

        auto g_x_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 192);

        auto g_x_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 193);

        auto g_x_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 194);

        #pragma omp simd aligned(cd_y, g_x_0_yyzz_xxxx, g_x_0_yyzz_xxxy, g_x_0_yyzz_xxxz, g_x_0_yyzz_xxyy, g_x_0_yyzz_xxyz, g_x_0_yyzz_xxzz, g_x_0_yyzz_xyyy, g_x_0_yyzz_xyyz, g_x_0_yyzz_xyzz, g_x_0_yyzz_xzzz, g_x_0_yyzz_yyyy, g_x_0_yyzz_yyyz, g_x_0_yyzz_yyzz, g_x_0_yyzz_yzzz, g_x_0_yyzz_zzzz, g_x_0_yzz_xxxx, g_x_0_yzz_xxxxy, g_x_0_yzz_xxxy, g_x_0_yzz_xxxyy, g_x_0_yzz_xxxyz, g_x_0_yzz_xxxz, g_x_0_yzz_xxyy, g_x_0_yzz_xxyyy, g_x_0_yzz_xxyyz, g_x_0_yzz_xxyz, g_x_0_yzz_xxyzz, g_x_0_yzz_xxzz, g_x_0_yzz_xyyy, g_x_0_yzz_xyyyy, g_x_0_yzz_xyyyz, g_x_0_yzz_xyyz, g_x_0_yzz_xyyzz, g_x_0_yzz_xyzz, g_x_0_yzz_xyzzz, g_x_0_yzz_xzzz, g_x_0_yzz_yyyy, g_x_0_yzz_yyyyy, g_x_0_yzz_yyyyz, g_x_0_yzz_yyyz, g_x_0_yzz_yyyzz, g_x_0_yzz_yyzz, g_x_0_yzz_yyzzz, g_x_0_yzz_yzzz, g_x_0_yzz_yzzzz, g_x_0_yzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyzz_xxxx[k] = -g_x_0_yzz_xxxx[k] * cd_y[k] + g_x_0_yzz_xxxxy[k];

            g_x_0_yyzz_xxxy[k] = -g_x_0_yzz_xxxy[k] * cd_y[k] + g_x_0_yzz_xxxyy[k];

            g_x_0_yyzz_xxxz[k] = -g_x_0_yzz_xxxz[k] * cd_y[k] + g_x_0_yzz_xxxyz[k];

            g_x_0_yyzz_xxyy[k] = -g_x_0_yzz_xxyy[k] * cd_y[k] + g_x_0_yzz_xxyyy[k];

            g_x_0_yyzz_xxyz[k] = -g_x_0_yzz_xxyz[k] * cd_y[k] + g_x_0_yzz_xxyyz[k];

            g_x_0_yyzz_xxzz[k] = -g_x_0_yzz_xxzz[k] * cd_y[k] + g_x_0_yzz_xxyzz[k];

            g_x_0_yyzz_xyyy[k] = -g_x_0_yzz_xyyy[k] * cd_y[k] + g_x_0_yzz_xyyyy[k];

            g_x_0_yyzz_xyyz[k] = -g_x_0_yzz_xyyz[k] * cd_y[k] + g_x_0_yzz_xyyyz[k];

            g_x_0_yyzz_xyzz[k] = -g_x_0_yzz_xyzz[k] * cd_y[k] + g_x_0_yzz_xyyzz[k];

            g_x_0_yyzz_xzzz[k] = -g_x_0_yzz_xzzz[k] * cd_y[k] + g_x_0_yzz_xyzzz[k];

            g_x_0_yyzz_yyyy[k] = -g_x_0_yzz_yyyy[k] * cd_y[k] + g_x_0_yzz_yyyyy[k];

            g_x_0_yyzz_yyyz[k] = -g_x_0_yzz_yyyz[k] * cd_y[k] + g_x_0_yzz_yyyyz[k];

            g_x_0_yyzz_yyzz[k] = -g_x_0_yzz_yyzz[k] * cd_y[k] + g_x_0_yzz_yyyzz[k];

            g_x_0_yyzz_yzzz[k] = -g_x_0_yzz_yzzz[k] * cd_y[k] + g_x_0_yzz_yyzzz[k];

            g_x_0_yyzz_zzzz[k] = -g_x_0_yzz_zzzz[k] * cd_y[k] + g_x_0_yzz_yzzzz[k];
        }

        /// Set up 195-210 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 195);

        auto g_x_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 196);

        auto g_x_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 197);

        auto g_x_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 198);

        auto g_x_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 199);

        auto g_x_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 200);

        auto g_x_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 201);

        auto g_x_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 202);

        auto g_x_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 203);

        auto g_x_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 204);

        auto g_x_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 205);

        auto g_x_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 206);

        auto g_x_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 207);

        auto g_x_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 208);

        auto g_x_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 209);

        #pragma omp simd aligned(cd_y, g_x_0_yzzz_xxxx, g_x_0_yzzz_xxxy, g_x_0_yzzz_xxxz, g_x_0_yzzz_xxyy, g_x_0_yzzz_xxyz, g_x_0_yzzz_xxzz, g_x_0_yzzz_xyyy, g_x_0_yzzz_xyyz, g_x_0_yzzz_xyzz, g_x_0_yzzz_xzzz, g_x_0_yzzz_yyyy, g_x_0_yzzz_yyyz, g_x_0_yzzz_yyzz, g_x_0_yzzz_yzzz, g_x_0_yzzz_zzzz, g_x_0_zzz_xxxx, g_x_0_zzz_xxxxy, g_x_0_zzz_xxxy, g_x_0_zzz_xxxyy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxz, g_x_0_zzz_xxyy, g_x_0_zzz_xxyyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzz, g_x_0_zzz_xyyy, g_x_0_zzz_xyyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzz, g_x_0_zzz_yyyy, g_x_0_zzz_yyyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzzz_xxxx[k] = -g_x_0_zzz_xxxx[k] * cd_y[k] + g_x_0_zzz_xxxxy[k];

            g_x_0_yzzz_xxxy[k] = -g_x_0_zzz_xxxy[k] * cd_y[k] + g_x_0_zzz_xxxyy[k];

            g_x_0_yzzz_xxxz[k] = -g_x_0_zzz_xxxz[k] * cd_y[k] + g_x_0_zzz_xxxyz[k];

            g_x_0_yzzz_xxyy[k] = -g_x_0_zzz_xxyy[k] * cd_y[k] + g_x_0_zzz_xxyyy[k];

            g_x_0_yzzz_xxyz[k] = -g_x_0_zzz_xxyz[k] * cd_y[k] + g_x_0_zzz_xxyyz[k];

            g_x_0_yzzz_xxzz[k] = -g_x_0_zzz_xxzz[k] * cd_y[k] + g_x_0_zzz_xxyzz[k];

            g_x_0_yzzz_xyyy[k] = -g_x_0_zzz_xyyy[k] * cd_y[k] + g_x_0_zzz_xyyyy[k];

            g_x_0_yzzz_xyyz[k] = -g_x_0_zzz_xyyz[k] * cd_y[k] + g_x_0_zzz_xyyyz[k];

            g_x_0_yzzz_xyzz[k] = -g_x_0_zzz_xyzz[k] * cd_y[k] + g_x_0_zzz_xyyzz[k];

            g_x_0_yzzz_xzzz[k] = -g_x_0_zzz_xzzz[k] * cd_y[k] + g_x_0_zzz_xyzzz[k];

            g_x_0_yzzz_yyyy[k] = -g_x_0_zzz_yyyy[k] * cd_y[k] + g_x_0_zzz_yyyyy[k];

            g_x_0_yzzz_yyyz[k] = -g_x_0_zzz_yyyz[k] * cd_y[k] + g_x_0_zzz_yyyyz[k];

            g_x_0_yzzz_yyzz[k] = -g_x_0_zzz_yyzz[k] * cd_y[k] + g_x_0_zzz_yyyzz[k];

            g_x_0_yzzz_yzzz[k] = -g_x_0_zzz_yzzz[k] * cd_y[k] + g_x_0_zzz_yyzzz[k];

            g_x_0_yzzz_zzzz[k] = -g_x_0_zzz_zzzz[k] * cd_y[k] + g_x_0_zzz_yzzzz[k];
        }

        /// Set up 210-225 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 0 * acomps  + 210);

        auto g_x_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 211);

        auto g_x_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 212);

        auto g_x_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 213);

        auto g_x_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 214);

        auto g_x_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 215);

        auto g_x_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 216);

        auto g_x_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 217);

        auto g_x_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 218);

        auto g_x_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 219);

        auto g_x_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 0 * acomps  + 220);

        auto g_x_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 221);

        auto g_x_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 222);

        auto g_x_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 223);

        auto g_x_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 0 * acomps  + 224);

        #pragma omp simd aligned(cd_z, g_x_0_zzz_xxxx, g_x_0_zzz_xxxxz, g_x_0_zzz_xxxy, g_x_0_zzz_xxxyz, g_x_0_zzz_xxxz, g_x_0_zzz_xxxzz, g_x_0_zzz_xxyy, g_x_0_zzz_xxyyz, g_x_0_zzz_xxyz, g_x_0_zzz_xxyzz, g_x_0_zzz_xxzz, g_x_0_zzz_xxzzz, g_x_0_zzz_xyyy, g_x_0_zzz_xyyyz, g_x_0_zzz_xyyz, g_x_0_zzz_xyyzz, g_x_0_zzz_xyzz, g_x_0_zzz_xyzzz, g_x_0_zzz_xzzz, g_x_0_zzz_xzzzz, g_x_0_zzz_yyyy, g_x_0_zzz_yyyyz, g_x_0_zzz_yyyz, g_x_0_zzz_yyyzz, g_x_0_zzz_yyzz, g_x_0_zzz_yyzzz, g_x_0_zzz_yzzz, g_x_0_zzz_yzzzz, g_x_0_zzz_zzzz, g_x_0_zzz_zzzzz, g_x_0_zzzz_xxxx, g_x_0_zzzz_xxxy, g_x_0_zzzz_xxxz, g_x_0_zzzz_xxyy, g_x_0_zzzz_xxyz, g_x_0_zzzz_xxzz, g_x_0_zzzz_xyyy, g_x_0_zzzz_xyyz, g_x_0_zzzz_xyzz, g_x_0_zzzz_xzzz, g_x_0_zzzz_yyyy, g_x_0_zzzz_yyyz, g_x_0_zzzz_yyzz, g_x_0_zzzz_yzzz, g_x_0_zzzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzzz_xxxx[k] = -g_x_0_zzz_xxxx[k] * cd_z[k] + g_x_0_zzz_xxxxz[k];

            g_x_0_zzzz_xxxy[k] = -g_x_0_zzz_xxxy[k] * cd_z[k] + g_x_0_zzz_xxxyz[k];

            g_x_0_zzzz_xxxz[k] = -g_x_0_zzz_xxxz[k] * cd_z[k] + g_x_0_zzz_xxxzz[k];

            g_x_0_zzzz_xxyy[k] = -g_x_0_zzz_xxyy[k] * cd_z[k] + g_x_0_zzz_xxyyz[k];

            g_x_0_zzzz_xxyz[k] = -g_x_0_zzz_xxyz[k] * cd_z[k] + g_x_0_zzz_xxyzz[k];

            g_x_0_zzzz_xxzz[k] = -g_x_0_zzz_xxzz[k] * cd_z[k] + g_x_0_zzz_xxzzz[k];

            g_x_0_zzzz_xyyy[k] = -g_x_0_zzz_xyyy[k] * cd_z[k] + g_x_0_zzz_xyyyz[k];

            g_x_0_zzzz_xyyz[k] = -g_x_0_zzz_xyyz[k] * cd_z[k] + g_x_0_zzz_xyyzz[k];

            g_x_0_zzzz_xyzz[k] = -g_x_0_zzz_xyzz[k] * cd_z[k] + g_x_0_zzz_xyzzz[k];

            g_x_0_zzzz_xzzz[k] = -g_x_0_zzz_xzzz[k] * cd_z[k] + g_x_0_zzz_xzzzz[k];

            g_x_0_zzzz_yyyy[k] = -g_x_0_zzz_yyyy[k] * cd_z[k] + g_x_0_zzz_yyyyz[k];

            g_x_0_zzzz_yyyz[k] = -g_x_0_zzz_yyyz[k] * cd_z[k] + g_x_0_zzz_yyyzz[k];

            g_x_0_zzzz_yyzz[k] = -g_x_0_zzz_yyzz[k] * cd_z[k] + g_x_0_zzz_yyzzz[k];

            g_x_0_zzzz_yzzz[k] = -g_x_0_zzz_yzzz[k] * cd_z[k] + g_x_0_zzz_yzzzz[k];

            g_x_0_zzzz_zzzz[k] = -g_x_0_zzz_zzzz[k] * cd_z[k] + g_x_0_zzz_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 0);

        auto g_y_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 1);

        auto g_y_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 2);

        auto g_y_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 3);

        auto g_y_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 4);

        auto g_y_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 5);

        auto g_y_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 6);

        auto g_y_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 7);

        auto g_y_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 8);

        auto g_y_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 9);

        auto g_y_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 10);

        auto g_y_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 11);

        auto g_y_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 12);

        auto g_y_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 13);

        auto g_y_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_y_0_xxx_xxxx, g_y_0_xxx_xxxxx, g_y_0_xxx_xxxxy, g_y_0_xxx_xxxxz, g_y_0_xxx_xxxy, g_y_0_xxx_xxxyy, g_y_0_xxx_xxxyz, g_y_0_xxx_xxxz, g_y_0_xxx_xxxzz, g_y_0_xxx_xxyy, g_y_0_xxx_xxyyy, g_y_0_xxx_xxyyz, g_y_0_xxx_xxyz, g_y_0_xxx_xxyzz, g_y_0_xxx_xxzz, g_y_0_xxx_xxzzz, g_y_0_xxx_xyyy, g_y_0_xxx_xyyyy, g_y_0_xxx_xyyyz, g_y_0_xxx_xyyz, g_y_0_xxx_xyyzz, g_y_0_xxx_xyzz, g_y_0_xxx_xyzzz, g_y_0_xxx_xzzz, g_y_0_xxx_xzzzz, g_y_0_xxx_yyyy, g_y_0_xxx_yyyz, g_y_0_xxx_yyzz, g_y_0_xxx_yzzz, g_y_0_xxx_zzzz, g_y_0_xxxx_xxxx, g_y_0_xxxx_xxxy, g_y_0_xxxx_xxxz, g_y_0_xxxx_xxyy, g_y_0_xxxx_xxyz, g_y_0_xxxx_xxzz, g_y_0_xxxx_xyyy, g_y_0_xxxx_xyyz, g_y_0_xxxx_xyzz, g_y_0_xxxx_xzzz, g_y_0_xxxx_yyyy, g_y_0_xxxx_yyyz, g_y_0_xxxx_yyzz, g_y_0_xxxx_yzzz, g_y_0_xxxx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxx_xxxx[k] = -g_y_0_xxx_xxxx[k] * cd_x[k] + g_y_0_xxx_xxxxx[k];

            g_y_0_xxxx_xxxy[k] = -g_y_0_xxx_xxxy[k] * cd_x[k] + g_y_0_xxx_xxxxy[k];

            g_y_0_xxxx_xxxz[k] = -g_y_0_xxx_xxxz[k] * cd_x[k] + g_y_0_xxx_xxxxz[k];

            g_y_0_xxxx_xxyy[k] = -g_y_0_xxx_xxyy[k] * cd_x[k] + g_y_0_xxx_xxxyy[k];

            g_y_0_xxxx_xxyz[k] = -g_y_0_xxx_xxyz[k] * cd_x[k] + g_y_0_xxx_xxxyz[k];

            g_y_0_xxxx_xxzz[k] = -g_y_0_xxx_xxzz[k] * cd_x[k] + g_y_0_xxx_xxxzz[k];

            g_y_0_xxxx_xyyy[k] = -g_y_0_xxx_xyyy[k] * cd_x[k] + g_y_0_xxx_xxyyy[k];

            g_y_0_xxxx_xyyz[k] = -g_y_0_xxx_xyyz[k] * cd_x[k] + g_y_0_xxx_xxyyz[k];

            g_y_0_xxxx_xyzz[k] = -g_y_0_xxx_xyzz[k] * cd_x[k] + g_y_0_xxx_xxyzz[k];

            g_y_0_xxxx_xzzz[k] = -g_y_0_xxx_xzzz[k] * cd_x[k] + g_y_0_xxx_xxzzz[k];

            g_y_0_xxxx_yyyy[k] = -g_y_0_xxx_yyyy[k] * cd_x[k] + g_y_0_xxx_xyyyy[k];

            g_y_0_xxxx_yyyz[k] = -g_y_0_xxx_yyyz[k] * cd_x[k] + g_y_0_xxx_xyyyz[k];

            g_y_0_xxxx_yyzz[k] = -g_y_0_xxx_yyzz[k] * cd_x[k] + g_y_0_xxx_xyyzz[k];

            g_y_0_xxxx_yzzz[k] = -g_y_0_xxx_yzzz[k] * cd_x[k] + g_y_0_xxx_xyzzz[k];

            g_y_0_xxxx_zzzz[k] = -g_y_0_xxx_zzzz[k] * cd_x[k] + g_y_0_xxx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 15);

        auto g_y_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 16);

        auto g_y_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 17);

        auto g_y_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 18);

        auto g_y_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 19);

        auto g_y_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 20);

        auto g_y_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 21);

        auto g_y_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 22);

        auto g_y_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 23);

        auto g_y_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 24);

        auto g_y_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 25);

        auto g_y_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 26);

        auto g_y_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 27);

        auto g_y_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 28);

        auto g_y_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xxxy_xxxx, g_y_0_xxxy_xxxy, g_y_0_xxxy_xxxz, g_y_0_xxxy_xxyy, g_y_0_xxxy_xxyz, g_y_0_xxxy_xxzz, g_y_0_xxxy_xyyy, g_y_0_xxxy_xyyz, g_y_0_xxxy_xyzz, g_y_0_xxxy_xzzz, g_y_0_xxxy_yyyy, g_y_0_xxxy_yyyz, g_y_0_xxxy_yyzz, g_y_0_xxxy_yzzz, g_y_0_xxxy_zzzz, g_y_0_xxy_xxxx, g_y_0_xxy_xxxxx, g_y_0_xxy_xxxxy, g_y_0_xxy_xxxxz, g_y_0_xxy_xxxy, g_y_0_xxy_xxxyy, g_y_0_xxy_xxxyz, g_y_0_xxy_xxxz, g_y_0_xxy_xxxzz, g_y_0_xxy_xxyy, g_y_0_xxy_xxyyy, g_y_0_xxy_xxyyz, g_y_0_xxy_xxyz, g_y_0_xxy_xxyzz, g_y_0_xxy_xxzz, g_y_0_xxy_xxzzz, g_y_0_xxy_xyyy, g_y_0_xxy_xyyyy, g_y_0_xxy_xyyyz, g_y_0_xxy_xyyz, g_y_0_xxy_xyyzz, g_y_0_xxy_xyzz, g_y_0_xxy_xyzzz, g_y_0_xxy_xzzz, g_y_0_xxy_xzzzz, g_y_0_xxy_yyyy, g_y_0_xxy_yyyz, g_y_0_xxy_yyzz, g_y_0_xxy_yzzz, g_y_0_xxy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxy_xxxx[k] = -g_y_0_xxy_xxxx[k] * cd_x[k] + g_y_0_xxy_xxxxx[k];

            g_y_0_xxxy_xxxy[k] = -g_y_0_xxy_xxxy[k] * cd_x[k] + g_y_0_xxy_xxxxy[k];

            g_y_0_xxxy_xxxz[k] = -g_y_0_xxy_xxxz[k] * cd_x[k] + g_y_0_xxy_xxxxz[k];

            g_y_0_xxxy_xxyy[k] = -g_y_0_xxy_xxyy[k] * cd_x[k] + g_y_0_xxy_xxxyy[k];

            g_y_0_xxxy_xxyz[k] = -g_y_0_xxy_xxyz[k] * cd_x[k] + g_y_0_xxy_xxxyz[k];

            g_y_0_xxxy_xxzz[k] = -g_y_0_xxy_xxzz[k] * cd_x[k] + g_y_0_xxy_xxxzz[k];

            g_y_0_xxxy_xyyy[k] = -g_y_0_xxy_xyyy[k] * cd_x[k] + g_y_0_xxy_xxyyy[k];

            g_y_0_xxxy_xyyz[k] = -g_y_0_xxy_xyyz[k] * cd_x[k] + g_y_0_xxy_xxyyz[k];

            g_y_0_xxxy_xyzz[k] = -g_y_0_xxy_xyzz[k] * cd_x[k] + g_y_0_xxy_xxyzz[k];

            g_y_0_xxxy_xzzz[k] = -g_y_0_xxy_xzzz[k] * cd_x[k] + g_y_0_xxy_xxzzz[k];

            g_y_0_xxxy_yyyy[k] = -g_y_0_xxy_yyyy[k] * cd_x[k] + g_y_0_xxy_xyyyy[k];

            g_y_0_xxxy_yyyz[k] = -g_y_0_xxy_yyyz[k] * cd_x[k] + g_y_0_xxy_xyyyz[k];

            g_y_0_xxxy_yyzz[k] = -g_y_0_xxy_yyzz[k] * cd_x[k] + g_y_0_xxy_xyyzz[k];

            g_y_0_xxxy_yzzz[k] = -g_y_0_xxy_yzzz[k] * cd_x[k] + g_y_0_xxy_xyzzz[k];

            g_y_0_xxxy_zzzz[k] = -g_y_0_xxy_zzzz[k] * cd_x[k] + g_y_0_xxy_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 30);

        auto g_y_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 31);

        auto g_y_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 32);

        auto g_y_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 33);

        auto g_y_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 34);

        auto g_y_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 35);

        auto g_y_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 36);

        auto g_y_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 37);

        auto g_y_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 38);

        auto g_y_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 39);

        auto g_y_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 40);

        auto g_y_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 41);

        auto g_y_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 42);

        auto g_y_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 43);

        auto g_y_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_y_0_xxxz_xxxx, g_y_0_xxxz_xxxy, g_y_0_xxxz_xxxz, g_y_0_xxxz_xxyy, g_y_0_xxxz_xxyz, g_y_0_xxxz_xxzz, g_y_0_xxxz_xyyy, g_y_0_xxxz_xyyz, g_y_0_xxxz_xyzz, g_y_0_xxxz_xzzz, g_y_0_xxxz_yyyy, g_y_0_xxxz_yyyz, g_y_0_xxxz_yyzz, g_y_0_xxxz_yzzz, g_y_0_xxxz_zzzz, g_y_0_xxz_xxxx, g_y_0_xxz_xxxxx, g_y_0_xxz_xxxxy, g_y_0_xxz_xxxxz, g_y_0_xxz_xxxy, g_y_0_xxz_xxxyy, g_y_0_xxz_xxxyz, g_y_0_xxz_xxxz, g_y_0_xxz_xxxzz, g_y_0_xxz_xxyy, g_y_0_xxz_xxyyy, g_y_0_xxz_xxyyz, g_y_0_xxz_xxyz, g_y_0_xxz_xxyzz, g_y_0_xxz_xxzz, g_y_0_xxz_xxzzz, g_y_0_xxz_xyyy, g_y_0_xxz_xyyyy, g_y_0_xxz_xyyyz, g_y_0_xxz_xyyz, g_y_0_xxz_xyyzz, g_y_0_xxz_xyzz, g_y_0_xxz_xyzzz, g_y_0_xxz_xzzz, g_y_0_xxz_xzzzz, g_y_0_xxz_yyyy, g_y_0_xxz_yyyz, g_y_0_xxz_yyzz, g_y_0_xxz_yzzz, g_y_0_xxz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxxz_xxxx[k] = -g_y_0_xxz_xxxx[k] * cd_x[k] + g_y_0_xxz_xxxxx[k];

            g_y_0_xxxz_xxxy[k] = -g_y_0_xxz_xxxy[k] * cd_x[k] + g_y_0_xxz_xxxxy[k];

            g_y_0_xxxz_xxxz[k] = -g_y_0_xxz_xxxz[k] * cd_x[k] + g_y_0_xxz_xxxxz[k];

            g_y_0_xxxz_xxyy[k] = -g_y_0_xxz_xxyy[k] * cd_x[k] + g_y_0_xxz_xxxyy[k];

            g_y_0_xxxz_xxyz[k] = -g_y_0_xxz_xxyz[k] * cd_x[k] + g_y_0_xxz_xxxyz[k];

            g_y_0_xxxz_xxzz[k] = -g_y_0_xxz_xxzz[k] * cd_x[k] + g_y_0_xxz_xxxzz[k];

            g_y_0_xxxz_xyyy[k] = -g_y_0_xxz_xyyy[k] * cd_x[k] + g_y_0_xxz_xxyyy[k];

            g_y_0_xxxz_xyyz[k] = -g_y_0_xxz_xyyz[k] * cd_x[k] + g_y_0_xxz_xxyyz[k];

            g_y_0_xxxz_xyzz[k] = -g_y_0_xxz_xyzz[k] * cd_x[k] + g_y_0_xxz_xxyzz[k];

            g_y_0_xxxz_xzzz[k] = -g_y_0_xxz_xzzz[k] * cd_x[k] + g_y_0_xxz_xxzzz[k];

            g_y_0_xxxz_yyyy[k] = -g_y_0_xxz_yyyy[k] * cd_x[k] + g_y_0_xxz_xyyyy[k];

            g_y_0_xxxz_yyyz[k] = -g_y_0_xxz_yyyz[k] * cd_x[k] + g_y_0_xxz_xyyyz[k];

            g_y_0_xxxz_yyzz[k] = -g_y_0_xxz_yyzz[k] * cd_x[k] + g_y_0_xxz_xyyzz[k];

            g_y_0_xxxz_yzzz[k] = -g_y_0_xxz_yzzz[k] * cd_x[k] + g_y_0_xxz_xyzzz[k];

            g_y_0_xxxz_zzzz[k] = -g_y_0_xxz_zzzz[k] * cd_x[k] + g_y_0_xxz_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 45);

        auto g_y_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 46);

        auto g_y_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 47);

        auto g_y_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 48);

        auto g_y_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 49);

        auto g_y_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 50);

        auto g_y_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 51);

        auto g_y_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 52);

        auto g_y_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 53);

        auto g_y_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 54);

        auto g_y_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 55);

        auto g_y_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 56);

        auto g_y_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 57);

        auto g_y_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 58);

        auto g_y_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_y_0_xxyy_xxxx, g_y_0_xxyy_xxxy, g_y_0_xxyy_xxxz, g_y_0_xxyy_xxyy, g_y_0_xxyy_xxyz, g_y_0_xxyy_xxzz, g_y_0_xxyy_xyyy, g_y_0_xxyy_xyyz, g_y_0_xxyy_xyzz, g_y_0_xxyy_xzzz, g_y_0_xxyy_yyyy, g_y_0_xxyy_yyyz, g_y_0_xxyy_yyzz, g_y_0_xxyy_yzzz, g_y_0_xxyy_zzzz, g_y_0_xyy_xxxx, g_y_0_xyy_xxxxx, g_y_0_xyy_xxxxy, g_y_0_xyy_xxxxz, g_y_0_xyy_xxxy, g_y_0_xyy_xxxyy, g_y_0_xyy_xxxyz, g_y_0_xyy_xxxz, g_y_0_xyy_xxxzz, g_y_0_xyy_xxyy, g_y_0_xyy_xxyyy, g_y_0_xyy_xxyyz, g_y_0_xyy_xxyz, g_y_0_xyy_xxyzz, g_y_0_xyy_xxzz, g_y_0_xyy_xxzzz, g_y_0_xyy_xyyy, g_y_0_xyy_xyyyy, g_y_0_xyy_xyyyz, g_y_0_xyy_xyyz, g_y_0_xyy_xyyzz, g_y_0_xyy_xyzz, g_y_0_xyy_xyzzz, g_y_0_xyy_xzzz, g_y_0_xyy_xzzzz, g_y_0_xyy_yyyy, g_y_0_xyy_yyyz, g_y_0_xyy_yyzz, g_y_0_xyy_yzzz, g_y_0_xyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyy_xxxx[k] = -g_y_0_xyy_xxxx[k] * cd_x[k] + g_y_0_xyy_xxxxx[k];

            g_y_0_xxyy_xxxy[k] = -g_y_0_xyy_xxxy[k] * cd_x[k] + g_y_0_xyy_xxxxy[k];

            g_y_0_xxyy_xxxz[k] = -g_y_0_xyy_xxxz[k] * cd_x[k] + g_y_0_xyy_xxxxz[k];

            g_y_0_xxyy_xxyy[k] = -g_y_0_xyy_xxyy[k] * cd_x[k] + g_y_0_xyy_xxxyy[k];

            g_y_0_xxyy_xxyz[k] = -g_y_0_xyy_xxyz[k] * cd_x[k] + g_y_0_xyy_xxxyz[k];

            g_y_0_xxyy_xxzz[k] = -g_y_0_xyy_xxzz[k] * cd_x[k] + g_y_0_xyy_xxxzz[k];

            g_y_0_xxyy_xyyy[k] = -g_y_0_xyy_xyyy[k] * cd_x[k] + g_y_0_xyy_xxyyy[k];

            g_y_0_xxyy_xyyz[k] = -g_y_0_xyy_xyyz[k] * cd_x[k] + g_y_0_xyy_xxyyz[k];

            g_y_0_xxyy_xyzz[k] = -g_y_0_xyy_xyzz[k] * cd_x[k] + g_y_0_xyy_xxyzz[k];

            g_y_0_xxyy_xzzz[k] = -g_y_0_xyy_xzzz[k] * cd_x[k] + g_y_0_xyy_xxzzz[k];

            g_y_0_xxyy_yyyy[k] = -g_y_0_xyy_yyyy[k] * cd_x[k] + g_y_0_xyy_xyyyy[k];

            g_y_0_xxyy_yyyz[k] = -g_y_0_xyy_yyyz[k] * cd_x[k] + g_y_0_xyy_xyyyz[k];

            g_y_0_xxyy_yyzz[k] = -g_y_0_xyy_yyzz[k] * cd_x[k] + g_y_0_xyy_xyyzz[k];

            g_y_0_xxyy_yzzz[k] = -g_y_0_xyy_yzzz[k] * cd_x[k] + g_y_0_xyy_xyzzz[k];

            g_y_0_xxyy_zzzz[k] = -g_y_0_xyy_zzzz[k] * cd_x[k] + g_y_0_xyy_xzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 60);

        auto g_y_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 61);

        auto g_y_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 62);

        auto g_y_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 63);

        auto g_y_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 64);

        auto g_y_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 65);

        auto g_y_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 66);

        auto g_y_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 67);

        auto g_y_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 68);

        auto g_y_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 69);

        auto g_y_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 70);

        auto g_y_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 71);

        auto g_y_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 72);

        auto g_y_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 73);

        auto g_y_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 74);

        #pragma omp simd aligned(cd_x, g_y_0_xxyz_xxxx, g_y_0_xxyz_xxxy, g_y_0_xxyz_xxxz, g_y_0_xxyz_xxyy, g_y_0_xxyz_xxyz, g_y_0_xxyz_xxzz, g_y_0_xxyz_xyyy, g_y_0_xxyz_xyyz, g_y_0_xxyz_xyzz, g_y_0_xxyz_xzzz, g_y_0_xxyz_yyyy, g_y_0_xxyz_yyyz, g_y_0_xxyz_yyzz, g_y_0_xxyz_yzzz, g_y_0_xxyz_zzzz, g_y_0_xyz_xxxx, g_y_0_xyz_xxxxx, g_y_0_xyz_xxxxy, g_y_0_xyz_xxxxz, g_y_0_xyz_xxxy, g_y_0_xyz_xxxyy, g_y_0_xyz_xxxyz, g_y_0_xyz_xxxz, g_y_0_xyz_xxxzz, g_y_0_xyz_xxyy, g_y_0_xyz_xxyyy, g_y_0_xyz_xxyyz, g_y_0_xyz_xxyz, g_y_0_xyz_xxyzz, g_y_0_xyz_xxzz, g_y_0_xyz_xxzzz, g_y_0_xyz_xyyy, g_y_0_xyz_xyyyy, g_y_0_xyz_xyyyz, g_y_0_xyz_xyyz, g_y_0_xyz_xyyzz, g_y_0_xyz_xyzz, g_y_0_xyz_xyzzz, g_y_0_xyz_xzzz, g_y_0_xyz_xzzzz, g_y_0_xyz_yyyy, g_y_0_xyz_yyyz, g_y_0_xyz_yyzz, g_y_0_xyz_yzzz, g_y_0_xyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxyz_xxxx[k] = -g_y_0_xyz_xxxx[k] * cd_x[k] + g_y_0_xyz_xxxxx[k];

            g_y_0_xxyz_xxxy[k] = -g_y_0_xyz_xxxy[k] * cd_x[k] + g_y_0_xyz_xxxxy[k];

            g_y_0_xxyz_xxxz[k] = -g_y_0_xyz_xxxz[k] * cd_x[k] + g_y_0_xyz_xxxxz[k];

            g_y_0_xxyz_xxyy[k] = -g_y_0_xyz_xxyy[k] * cd_x[k] + g_y_0_xyz_xxxyy[k];

            g_y_0_xxyz_xxyz[k] = -g_y_0_xyz_xxyz[k] * cd_x[k] + g_y_0_xyz_xxxyz[k];

            g_y_0_xxyz_xxzz[k] = -g_y_0_xyz_xxzz[k] * cd_x[k] + g_y_0_xyz_xxxzz[k];

            g_y_0_xxyz_xyyy[k] = -g_y_0_xyz_xyyy[k] * cd_x[k] + g_y_0_xyz_xxyyy[k];

            g_y_0_xxyz_xyyz[k] = -g_y_0_xyz_xyyz[k] * cd_x[k] + g_y_0_xyz_xxyyz[k];

            g_y_0_xxyz_xyzz[k] = -g_y_0_xyz_xyzz[k] * cd_x[k] + g_y_0_xyz_xxyzz[k];

            g_y_0_xxyz_xzzz[k] = -g_y_0_xyz_xzzz[k] * cd_x[k] + g_y_0_xyz_xxzzz[k];

            g_y_0_xxyz_yyyy[k] = -g_y_0_xyz_yyyy[k] * cd_x[k] + g_y_0_xyz_xyyyy[k];

            g_y_0_xxyz_yyyz[k] = -g_y_0_xyz_yyyz[k] * cd_x[k] + g_y_0_xyz_xyyyz[k];

            g_y_0_xxyz_yyzz[k] = -g_y_0_xyz_yyzz[k] * cd_x[k] + g_y_0_xyz_xyyzz[k];

            g_y_0_xxyz_yzzz[k] = -g_y_0_xyz_yzzz[k] * cd_x[k] + g_y_0_xyz_xyzzz[k];

            g_y_0_xxyz_zzzz[k] = -g_y_0_xyz_zzzz[k] * cd_x[k] + g_y_0_xyz_xzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 75);

        auto g_y_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 76);

        auto g_y_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 77);

        auto g_y_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 78);

        auto g_y_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 79);

        auto g_y_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 80);

        auto g_y_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 81);

        auto g_y_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 82);

        auto g_y_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 83);

        auto g_y_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 84);

        auto g_y_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 85);

        auto g_y_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 86);

        auto g_y_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 87);

        auto g_y_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 88);

        auto g_y_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 89);

        #pragma omp simd aligned(cd_x, g_y_0_xxzz_xxxx, g_y_0_xxzz_xxxy, g_y_0_xxzz_xxxz, g_y_0_xxzz_xxyy, g_y_0_xxzz_xxyz, g_y_0_xxzz_xxzz, g_y_0_xxzz_xyyy, g_y_0_xxzz_xyyz, g_y_0_xxzz_xyzz, g_y_0_xxzz_xzzz, g_y_0_xxzz_yyyy, g_y_0_xxzz_yyyz, g_y_0_xxzz_yyzz, g_y_0_xxzz_yzzz, g_y_0_xxzz_zzzz, g_y_0_xzz_xxxx, g_y_0_xzz_xxxxx, g_y_0_xzz_xxxxy, g_y_0_xzz_xxxxz, g_y_0_xzz_xxxy, g_y_0_xzz_xxxyy, g_y_0_xzz_xxxyz, g_y_0_xzz_xxxz, g_y_0_xzz_xxxzz, g_y_0_xzz_xxyy, g_y_0_xzz_xxyyy, g_y_0_xzz_xxyyz, g_y_0_xzz_xxyz, g_y_0_xzz_xxyzz, g_y_0_xzz_xxzz, g_y_0_xzz_xxzzz, g_y_0_xzz_xyyy, g_y_0_xzz_xyyyy, g_y_0_xzz_xyyyz, g_y_0_xzz_xyyz, g_y_0_xzz_xyyzz, g_y_0_xzz_xyzz, g_y_0_xzz_xyzzz, g_y_0_xzz_xzzz, g_y_0_xzz_xzzzz, g_y_0_xzz_yyyy, g_y_0_xzz_yyyz, g_y_0_xzz_yyzz, g_y_0_xzz_yzzz, g_y_0_xzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxzz_xxxx[k] = -g_y_0_xzz_xxxx[k] * cd_x[k] + g_y_0_xzz_xxxxx[k];

            g_y_0_xxzz_xxxy[k] = -g_y_0_xzz_xxxy[k] * cd_x[k] + g_y_0_xzz_xxxxy[k];

            g_y_0_xxzz_xxxz[k] = -g_y_0_xzz_xxxz[k] * cd_x[k] + g_y_0_xzz_xxxxz[k];

            g_y_0_xxzz_xxyy[k] = -g_y_0_xzz_xxyy[k] * cd_x[k] + g_y_0_xzz_xxxyy[k];

            g_y_0_xxzz_xxyz[k] = -g_y_0_xzz_xxyz[k] * cd_x[k] + g_y_0_xzz_xxxyz[k];

            g_y_0_xxzz_xxzz[k] = -g_y_0_xzz_xxzz[k] * cd_x[k] + g_y_0_xzz_xxxzz[k];

            g_y_0_xxzz_xyyy[k] = -g_y_0_xzz_xyyy[k] * cd_x[k] + g_y_0_xzz_xxyyy[k];

            g_y_0_xxzz_xyyz[k] = -g_y_0_xzz_xyyz[k] * cd_x[k] + g_y_0_xzz_xxyyz[k];

            g_y_0_xxzz_xyzz[k] = -g_y_0_xzz_xyzz[k] * cd_x[k] + g_y_0_xzz_xxyzz[k];

            g_y_0_xxzz_xzzz[k] = -g_y_0_xzz_xzzz[k] * cd_x[k] + g_y_0_xzz_xxzzz[k];

            g_y_0_xxzz_yyyy[k] = -g_y_0_xzz_yyyy[k] * cd_x[k] + g_y_0_xzz_xyyyy[k];

            g_y_0_xxzz_yyyz[k] = -g_y_0_xzz_yyyz[k] * cd_x[k] + g_y_0_xzz_xyyyz[k];

            g_y_0_xxzz_yyzz[k] = -g_y_0_xzz_yyzz[k] * cd_x[k] + g_y_0_xzz_xyyzz[k];

            g_y_0_xxzz_yzzz[k] = -g_y_0_xzz_yzzz[k] * cd_x[k] + g_y_0_xzz_xyzzz[k];

            g_y_0_xxzz_zzzz[k] = -g_y_0_xzz_zzzz[k] * cd_x[k] + g_y_0_xzz_xzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 90);

        auto g_y_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 91);

        auto g_y_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 92);

        auto g_y_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 93);

        auto g_y_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 94);

        auto g_y_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 95);

        auto g_y_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 96);

        auto g_y_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 97);

        auto g_y_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 98);

        auto g_y_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 99);

        auto g_y_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 100);

        auto g_y_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 101);

        auto g_y_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 102);

        auto g_y_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 103);

        auto g_y_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 104);

        #pragma omp simd aligned(cd_x, g_y_0_xyyy_xxxx, g_y_0_xyyy_xxxy, g_y_0_xyyy_xxxz, g_y_0_xyyy_xxyy, g_y_0_xyyy_xxyz, g_y_0_xyyy_xxzz, g_y_0_xyyy_xyyy, g_y_0_xyyy_xyyz, g_y_0_xyyy_xyzz, g_y_0_xyyy_xzzz, g_y_0_xyyy_yyyy, g_y_0_xyyy_yyyz, g_y_0_xyyy_yyzz, g_y_0_xyyy_yzzz, g_y_0_xyyy_zzzz, g_y_0_yyy_xxxx, g_y_0_yyy_xxxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyz, g_y_0_yyy_yyzz, g_y_0_yyy_yzzz, g_y_0_yyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyy_xxxx[k] = -g_y_0_yyy_xxxx[k] * cd_x[k] + g_y_0_yyy_xxxxx[k];

            g_y_0_xyyy_xxxy[k] = -g_y_0_yyy_xxxy[k] * cd_x[k] + g_y_0_yyy_xxxxy[k];

            g_y_0_xyyy_xxxz[k] = -g_y_0_yyy_xxxz[k] * cd_x[k] + g_y_0_yyy_xxxxz[k];

            g_y_0_xyyy_xxyy[k] = -g_y_0_yyy_xxyy[k] * cd_x[k] + g_y_0_yyy_xxxyy[k];

            g_y_0_xyyy_xxyz[k] = -g_y_0_yyy_xxyz[k] * cd_x[k] + g_y_0_yyy_xxxyz[k];

            g_y_0_xyyy_xxzz[k] = -g_y_0_yyy_xxzz[k] * cd_x[k] + g_y_0_yyy_xxxzz[k];

            g_y_0_xyyy_xyyy[k] = -g_y_0_yyy_xyyy[k] * cd_x[k] + g_y_0_yyy_xxyyy[k];

            g_y_0_xyyy_xyyz[k] = -g_y_0_yyy_xyyz[k] * cd_x[k] + g_y_0_yyy_xxyyz[k];

            g_y_0_xyyy_xyzz[k] = -g_y_0_yyy_xyzz[k] * cd_x[k] + g_y_0_yyy_xxyzz[k];

            g_y_0_xyyy_xzzz[k] = -g_y_0_yyy_xzzz[k] * cd_x[k] + g_y_0_yyy_xxzzz[k];

            g_y_0_xyyy_yyyy[k] = -g_y_0_yyy_yyyy[k] * cd_x[k] + g_y_0_yyy_xyyyy[k];

            g_y_0_xyyy_yyyz[k] = -g_y_0_yyy_yyyz[k] * cd_x[k] + g_y_0_yyy_xyyyz[k];

            g_y_0_xyyy_yyzz[k] = -g_y_0_yyy_yyzz[k] * cd_x[k] + g_y_0_yyy_xyyzz[k];

            g_y_0_xyyy_yzzz[k] = -g_y_0_yyy_yzzz[k] * cd_x[k] + g_y_0_yyy_xyzzz[k];

            g_y_0_xyyy_zzzz[k] = -g_y_0_yyy_zzzz[k] * cd_x[k] + g_y_0_yyy_xzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 105);

        auto g_y_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 106);

        auto g_y_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 107);

        auto g_y_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 108);

        auto g_y_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 109);

        auto g_y_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 110);

        auto g_y_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 111);

        auto g_y_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 112);

        auto g_y_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 113);

        auto g_y_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 114);

        auto g_y_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 115);

        auto g_y_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 116);

        auto g_y_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 117);

        auto g_y_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 118);

        auto g_y_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 119);

        #pragma omp simd aligned(cd_x, g_y_0_xyyz_xxxx, g_y_0_xyyz_xxxy, g_y_0_xyyz_xxxz, g_y_0_xyyz_xxyy, g_y_0_xyyz_xxyz, g_y_0_xyyz_xxzz, g_y_0_xyyz_xyyy, g_y_0_xyyz_xyyz, g_y_0_xyyz_xyzz, g_y_0_xyyz_xzzz, g_y_0_xyyz_yyyy, g_y_0_xyyz_yyyz, g_y_0_xyyz_yyzz, g_y_0_xyyz_yzzz, g_y_0_xyyz_zzzz, g_y_0_yyz_xxxx, g_y_0_yyz_xxxxx, g_y_0_yyz_xxxxy, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxy, g_y_0_yyz_xxxyy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyy, g_y_0_yyz_xxyyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyy, g_y_0_yyz_xyyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyy, g_y_0_yyz_yyyz, g_y_0_yyz_yyzz, g_y_0_yyz_yzzz, g_y_0_yyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyyz_xxxx[k] = -g_y_0_yyz_xxxx[k] * cd_x[k] + g_y_0_yyz_xxxxx[k];

            g_y_0_xyyz_xxxy[k] = -g_y_0_yyz_xxxy[k] * cd_x[k] + g_y_0_yyz_xxxxy[k];

            g_y_0_xyyz_xxxz[k] = -g_y_0_yyz_xxxz[k] * cd_x[k] + g_y_0_yyz_xxxxz[k];

            g_y_0_xyyz_xxyy[k] = -g_y_0_yyz_xxyy[k] * cd_x[k] + g_y_0_yyz_xxxyy[k];

            g_y_0_xyyz_xxyz[k] = -g_y_0_yyz_xxyz[k] * cd_x[k] + g_y_0_yyz_xxxyz[k];

            g_y_0_xyyz_xxzz[k] = -g_y_0_yyz_xxzz[k] * cd_x[k] + g_y_0_yyz_xxxzz[k];

            g_y_0_xyyz_xyyy[k] = -g_y_0_yyz_xyyy[k] * cd_x[k] + g_y_0_yyz_xxyyy[k];

            g_y_0_xyyz_xyyz[k] = -g_y_0_yyz_xyyz[k] * cd_x[k] + g_y_0_yyz_xxyyz[k];

            g_y_0_xyyz_xyzz[k] = -g_y_0_yyz_xyzz[k] * cd_x[k] + g_y_0_yyz_xxyzz[k];

            g_y_0_xyyz_xzzz[k] = -g_y_0_yyz_xzzz[k] * cd_x[k] + g_y_0_yyz_xxzzz[k];

            g_y_0_xyyz_yyyy[k] = -g_y_0_yyz_yyyy[k] * cd_x[k] + g_y_0_yyz_xyyyy[k];

            g_y_0_xyyz_yyyz[k] = -g_y_0_yyz_yyyz[k] * cd_x[k] + g_y_0_yyz_xyyyz[k];

            g_y_0_xyyz_yyzz[k] = -g_y_0_yyz_yyzz[k] * cd_x[k] + g_y_0_yyz_xyyzz[k];

            g_y_0_xyyz_yzzz[k] = -g_y_0_yyz_yzzz[k] * cd_x[k] + g_y_0_yyz_xyzzz[k];

            g_y_0_xyyz_zzzz[k] = -g_y_0_yyz_zzzz[k] * cd_x[k] + g_y_0_yyz_xzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 120);

        auto g_y_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 121);

        auto g_y_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 122);

        auto g_y_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 123);

        auto g_y_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 124);

        auto g_y_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 125);

        auto g_y_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 126);

        auto g_y_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 127);

        auto g_y_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 128);

        auto g_y_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 129);

        auto g_y_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 130);

        auto g_y_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 131);

        auto g_y_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 132);

        auto g_y_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 133);

        auto g_y_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 134);

        #pragma omp simd aligned(cd_x, g_y_0_xyzz_xxxx, g_y_0_xyzz_xxxy, g_y_0_xyzz_xxxz, g_y_0_xyzz_xxyy, g_y_0_xyzz_xxyz, g_y_0_xyzz_xxzz, g_y_0_xyzz_xyyy, g_y_0_xyzz_xyyz, g_y_0_xyzz_xyzz, g_y_0_xyzz_xzzz, g_y_0_xyzz_yyyy, g_y_0_xyzz_yyyz, g_y_0_xyzz_yyzz, g_y_0_xyzz_yzzz, g_y_0_xyzz_zzzz, g_y_0_yzz_xxxx, g_y_0_yzz_xxxxx, g_y_0_yzz_xxxxy, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxy, g_y_0_yzz_xxxyy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyy, g_y_0_yzz_xxyyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyy, g_y_0_yzz_xyyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyy, g_y_0_yzz_yyyz, g_y_0_yzz_yyzz, g_y_0_yzz_yzzz, g_y_0_yzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyzz_xxxx[k] = -g_y_0_yzz_xxxx[k] * cd_x[k] + g_y_0_yzz_xxxxx[k];

            g_y_0_xyzz_xxxy[k] = -g_y_0_yzz_xxxy[k] * cd_x[k] + g_y_0_yzz_xxxxy[k];

            g_y_0_xyzz_xxxz[k] = -g_y_0_yzz_xxxz[k] * cd_x[k] + g_y_0_yzz_xxxxz[k];

            g_y_0_xyzz_xxyy[k] = -g_y_0_yzz_xxyy[k] * cd_x[k] + g_y_0_yzz_xxxyy[k];

            g_y_0_xyzz_xxyz[k] = -g_y_0_yzz_xxyz[k] * cd_x[k] + g_y_0_yzz_xxxyz[k];

            g_y_0_xyzz_xxzz[k] = -g_y_0_yzz_xxzz[k] * cd_x[k] + g_y_0_yzz_xxxzz[k];

            g_y_0_xyzz_xyyy[k] = -g_y_0_yzz_xyyy[k] * cd_x[k] + g_y_0_yzz_xxyyy[k];

            g_y_0_xyzz_xyyz[k] = -g_y_0_yzz_xyyz[k] * cd_x[k] + g_y_0_yzz_xxyyz[k];

            g_y_0_xyzz_xyzz[k] = -g_y_0_yzz_xyzz[k] * cd_x[k] + g_y_0_yzz_xxyzz[k];

            g_y_0_xyzz_xzzz[k] = -g_y_0_yzz_xzzz[k] * cd_x[k] + g_y_0_yzz_xxzzz[k];

            g_y_0_xyzz_yyyy[k] = -g_y_0_yzz_yyyy[k] * cd_x[k] + g_y_0_yzz_xyyyy[k];

            g_y_0_xyzz_yyyz[k] = -g_y_0_yzz_yyyz[k] * cd_x[k] + g_y_0_yzz_xyyyz[k];

            g_y_0_xyzz_yyzz[k] = -g_y_0_yzz_yyzz[k] * cd_x[k] + g_y_0_yzz_xyyzz[k];

            g_y_0_xyzz_yzzz[k] = -g_y_0_yzz_yzzz[k] * cd_x[k] + g_y_0_yzz_xyzzz[k];

            g_y_0_xyzz_zzzz[k] = -g_y_0_yzz_zzzz[k] * cd_x[k] + g_y_0_yzz_xzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 135);

        auto g_y_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 136);

        auto g_y_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 137);

        auto g_y_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 138);

        auto g_y_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 139);

        auto g_y_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 140);

        auto g_y_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 141);

        auto g_y_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 142);

        auto g_y_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 143);

        auto g_y_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 144);

        auto g_y_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 145);

        auto g_y_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 146);

        auto g_y_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 147);

        auto g_y_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 148);

        auto g_y_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 149);

        #pragma omp simd aligned(cd_x, g_y_0_xzzz_xxxx, g_y_0_xzzz_xxxy, g_y_0_xzzz_xxxz, g_y_0_xzzz_xxyy, g_y_0_xzzz_xxyz, g_y_0_xzzz_xxzz, g_y_0_xzzz_xyyy, g_y_0_xzzz_xyyz, g_y_0_xzzz_xyzz, g_y_0_xzzz_xzzz, g_y_0_xzzz_yyyy, g_y_0_xzzz_yyyz, g_y_0_xzzz_yyzz, g_y_0_xzzz_yzzz, g_y_0_xzzz_zzzz, g_y_0_zzz_xxxx, g_y_0_zzz_xxxxx, g_y_0_zzz_xxxxy, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxy, g_y_0_zzz_xxxyy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyy, g_y_0_zzz_xxyyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyy, g_y_0_zzz_xyyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyy, g_y_0_zzz_yyyz, g_y_0_zzz_yyzz, g_y_0_zzz_yzzz, g_y_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzzz_xxxx[k] = -g_y_0_zzz_xxxx[k] * cd_x[k] + g_y_0_zzz_xxxxx[k];

            g_y_0_xzzz_xxxy[k] = -g_y_0_zzz_xxxy[k] * cd_x[k] + g_y_0_zzz_xxxxy[k];

            g_y_0_xzzz_xxxz[k] = -g_y_0_zzz_xxxz[k] * cd_x[k] + g_y_0_zzz_xxxxz[k];

            g_y_0_xzzz_xxyy[k] = -g_y_0_zzz_xxyy[k] * cd_x[k] + g_y_0_zzz_xxxyy[k];

            g_y_0_xzzz_xxyz[k] = -g_y_0_zzz_xxyz[k] * cd_x[k] + g_y_0_zzz_xxxyz[k];

            g_y_0_xzzz_xxzz[k] = -g_y_0_zzz_xxzz[k] * cd_x[k] + g_y_0_zzz_xxxzz[k];

            g_y_0_xzzz_xyyy[k] = -g_y_0_zzz_xyyy[k] * cd_x[k] + g_y_0_zzz_xxyyy[k];

            g_y_0_xzzz_xyyz[k] = -g_y_0_zzz_xyyz[k] * cd_x[k] + g_y_0_zzz_xxyyz[k];

            g_y_0_xzzz_xyzz[k] = -g_y_0_zzz_xyzz[k] * cd_x[k] + g_y_0_zzz_xxyzz[k];

            g_y_0_xzzz_xzzz[k] = -g_y_0_zzz_xzzz[k] * cd_x[k] + g_y_0_zzz_xxzzz[k];

            g_y_0_xzzz_yyyy[k] = -g_y_0_zzz_yyyy[k] * cd_x[k] + g_y_0_zzz_xyyyy[k];

            g_y_0_xzzz_yyyz[k] = -g_y_0_zzz_yyyz[k] * cd_x[k] + g_y_0_zzz_xyyyz[k];

            g_y_0_xzzz_yyzz[k] = -g_y_0_zzz_yyzz[k] * cd_x[k] + g_y_0_zzz_xyyzz[k];

            g_y_0_xzzz_yzzz[k] = -g_y_0_zzz_yzzz[k] * cd_x[k] + g_y_0_zzz_xyzzz[k];

            g_y_0_xzzz_zzzz[k] = -g_y_0_zzz_zzzz[k] * cd_x[k] + g_y_0_zzz_xzzzz[k];
        }

        /// Set up 150-165 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 150);

        auto g_y_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 151);

        auto g_y_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 152);

        auto g_y_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 153);

        auto g_y_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 154);

        auto g_y_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 155);

        auto g_y_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 156);

        auto g_y_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 157);

        auto g_y_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 158);

        auto g_y_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 159);

        auto g_y_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 160);

        auto g_y_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 161);

        auto g_y_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 162);

        auto g_y_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 163);

        auto g_y_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 164);

        #pragma omp simd aligned(cd_y, g_y_0_yyy_xxxx, g_y_0_yyy_xxxxy, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzz, g_y_0_yyyy_xxxx, g_y_0_yyyy_xxxy, g_y_0_yyyy_xxxz, g_y_0_yyyy_xxyy, g_y_0_yyyy_xxyz, g_y_0_yyyy_xxzz, g_y_0_yyyy_xyyy, g_y_0_yyyy_xyyz, g_y_0_yyyy_xyzz, g_y_0_yyyy_xzzz, g_y_0_yyyy_yyyy, g_y_0_yyyy_yyyz, g_y_0_yyyy_yyzz, g_y_0_yyyy_yzzz, g_y_0_yyyy_zzzz, g_yyy_xxxx, g_yyy_xxxy, g_yyy_xxxz, g_yyy_xxyy, g_yyy_xxyz, g_yyy_xxzz, g_yyy_xyyy, g_yyy_xyyz, g_yyy_xyzz, g_yyy_xzzz, g_yyy_yyyy, g_yyy_yyyz, g_yyy_yyzz, g_yyy_yzzz, g_yyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyy_xxxx[k] = -g_yyy_xxxx[k] - g_y_0_yyy_xxxx[k] * cd_y[k] + g_y_0_yyy_xxxxy[k];

            g_y_0_yyyy_xxxy[k] = -g_yyy_xxxy[k] - g_y_0_yyy_xxxy[k] * cd_y[k] + g_y_0_yyy_xxxyy[k];

            g_y_0_yyyy_xxxz[k] = -g_yyy_xxxz[k] - g_y_0_yyy_xxxz[k] * cd_y[k] + g_y_0_yyy_xxxyz[k];

            g_y_0_yyyy_xxyy[k] = -g_yyy_xxyy[k] - g_y_0_yyy_xxyy[k] * cd_y[k] + g_y_0_yyy_xxyyy[k];

            g_y_0_yyyy_xxyz[k] = -g_yyy_xxyz[k] - g_y_0_yyy_xxyz[k] * cd_y[k] + g_y_0_yyy_xxyyz[k];

            g_y_0_yyyy_xxzz[k] = -g_yyy_xxzz[k] - g_y_0_yyy_xxzz[k] * cd_y[k] + g_y_0_yyy_xxyzz[k];

            g_y_0_yyyy_xyyy[k] = -g_yyy_xyyy[k] - g_y_0_yyy_xyyy[k] * cd_y[k] + g_y_0_yyy_xyyyy[k];

            g_y_0_yyyy_xyyz[k] = -g_yyy_xyyz[k] - g_y_0_yyy_xyyz[k] * cd_y[k] + g_y_0_yyy_xyyyz[k];

            g_y_0_yyyy_xyzz[k] = -g_yyy_xyzz[k] - g_y_0_yyy_xyzz[k] * cd_y[k] + g_y_0_yyy_xyyzz[k];

            g_y_0_yyyy_xzzz[k] = -g_yyy_xzzz[k] - g_y_0_yyy_xzzz[k] * cd_y[k] + g_y_0_yyy_xyzzz[k];

            g_y_0_yyyy_yyyy[k] = -g_yyy_yyyy[k] - g_y_0_yyy_yyyy[k] * cd_y[k] + g_y_0_yyy_yyyyy[k];

            g_y_0_yyyy_yyyz[k] = -g_yyy_yyyz[k] - g_y_0_yyy_yyyz[k] * cd_y[k] + g_y_0_yyy_yyyyz[k];

            g_y_0_yyyy_yyzz[k] = -g_yyy_yyzz[k] - g_y_0_yyy_yyzz[k] * cd_y[k] + g_y_0_yyy_yyyzz[k];

            g_y_0_yyyy_yzzz[k] = -g_yyy_yzzz[k] - g_y_0_yyy_yzzz[k] * cd_y[k] + g_y_0_yyy_yyzzz[k];

            g_y_0_yyyy_zzzz[k] = -g_yyy_zzzz[k] - g_y_0_yyy_zzzz[k] * cd_y[k] + g_y_0_yyy_yzzzz[k];
        }

        /// Set up 165-180 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 165);

        auto g_y_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 166);

        auto g_y_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 167);

        auto g_y_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 168);

        auto g_y_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 169);

        auto g_y_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 170);

        auto g_y_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 171);

        auto g_y_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 172);

        auto g_y_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 173);

        auto g_y_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 174);

        auto g_y_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 175);

        auto g_y_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 176);

        auto g_y_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 177);

        auto g_y_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 178);

        auto g_y_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 179);

        #pragma omp simd aligned(cd_z, g_y_0_yyy_xxxx, g_y_0_yyy_xxxxz, g_y_0_yyy_xxxy, g_y_0_yyy_xxxyz, g_y_0_yyy_xxxz, g_y_0_yyy_xxxzz, g_y_0_yyy_xxyy, g_y_0_yyy_xxyyz, g_y_0_yyy_xxyz, g_y_0_yyy_xxyzz, g_y_0_yyy_xxzz, g_y_0_yyy_xxzzz, g_y_0_yyy_xyyy, g_y_0_yyy_xyyyz, g_y_0_yyy_xyyz, g_y_0_yyy_xyyzz, g_y_0_yyy_xyzz, g_y_0_yyy_xyzzz, g_y_0_yyy_xzzz, g_y_0_yyy_xzzzz, g_y_0_yyy_yyyy, g_y_0_yyy_yyyyz, g_y_0_yyy_yyyz, g_y_0_yyy_yyyzz, g_y_0_yyy_yyzz, g_y_0_yyy_yyzzz, g_y_0_yyy_yzzz, g_y_0_yyy_yzzzz, g_y_0_yyy_zzzz, g_y_0_yyy_zzzzz, g_y_0_yyyz_xxxx, g_y_0_yyyz_xxxy, g_y_0_yyyz_xxxz, g_y_0_yyyz_xxyy, g_y_0_yyyz_xxyz, g_y_0_yyyz_xxzz, g_y_0_yyyz_xyyy, g_y_0_yyyz_xyyz, g_y_0_yyyz_xyzz, g_y_0_yyyz_xzzz, g_y_0_yyyz_yyyy, g_y_0_yyyz_yyyz, g_y_0_yyyz_yyzz, g_y_0_yyyz_yzzz, g_y_0_yyyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyyz_xxxx[k] = -g_y_0_yyy_xxxx[k] * cd_z[k] + g_y_0_yyy_xxxxz[k];

            g_y_0_yyyz_xxxy[k] = -g_y_0_yyy_xxxy[k] * cd_z[k] + g_y_0_yyy_xxxyz[k];

            g_y_0_yyyz_xxxz[k] = -g_y_0_yyy_xxxz[k] * cd_z[k] + g_y_0_yyy_xxxzz[k];

            g_y_0_yyyz_xxyy[k] = -g_y_0_yyy_xxyy[k] * cd_z[k] + g_y_0_yyy_xxyyz[k];

            g_y_0_yyyz_xxyz[k] = -g_y_0_yyy_xxyz[k] * cd_z[k] + g_y_0_yyy_xxyzz[k];

            g_y_0_yyyz_xxzz[k] = -g_y_0_yyy_xxzz[k] * cd_z[k] + g_y_0_yyy_xxzzz[k];

            g_y_0_yyyz_xyyy[k] = -g_y_0_yyy_xyyy[k] * cd_z[k] + g_y_0_yyy_xyyyz[k];

            g_y_0_yyyz_xyyz[k] = -g_y_0_yyy_xyyz[k] * cd_z[k] + g_y_0_yyy_xyyzz[k];

            g_y_0_yyyz_xyzz[k] = -g_y_0_yyy_xyzz[k] * cd_z[k] + g_y_0_yyy_xyzzz[k];

            g_y_0_yyyz_xzzz[k] = -g_y_0_yyy_xzzz[k] * cd_z[k] + g_y_0_yyy_xzzzz[k];

            g_y_0_yyyz_yyyy[k] = -g_y_0_yyy_yyyy[k] * cd_z[k] + g_y_0_yyy_yyyyz[k];

            g_y_0_yyyz_yyyz[k] = -g_y_0_yyy_yyyz[k] * cd_z[k] + g_y_0_yyy_yyyzz[k];

            g_y_0_yyyz_yyzz[k] = -g_y_0_yyy_yyzz[k] * cd_z[k] + g_y_0_yyy_yyzzz[k];

            g_y_0_yyyz_yzzz[k] = -g_y_0_yyy_yzzz[k] * cd_z[k] + g_y_0_yyy_yzzzz[k];

            g_y_0_yyyz_zzzz[k] = -g_y_0_yyy_zzzz[k] * cd_z[k] + g_y_0_yyy_zzzzz[k];
        }

        /// Set up 180-195 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 180);

        auto g_y_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 181);

        auto g_y_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 182);

        auto g_y_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 183);

        auto g_y_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 184);

        auto g_y_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 185);

        auto g_y_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 186);

        auto g_y_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 187);

        auto g_y_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 188);

        auto g_y_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 189);

        auto g_y_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 190);

        auto g_y_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 191);

        auto g_y_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 192);

        auto g_y_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 193);

        auto g_y_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 194);

        #pragma omp simd aligned(cd_z, g_y_0_yyz_xxxx, g_y_0_yyz_xxxxz, g_y_0_yyz_xxxy, g_y_0_yyz_xxxyz, g_y_0_yyz_xxxz, g_y_0_yyz_xxxzz, g_y_0_yyz_xxyy, g_y_0_yyz_xxyyz, g_y_0_yyz_xxyz, g_y_0_yyz_xxyzz, g_y_0_yyz_xxzz, g_y_0_yyz_xxzzz, g_y_0_yyz_xyyy, g_y_0_yyz_xyyyz, g_y_0_yyz_xyyz, g_y_0_yyz_xyyzz, g_y_0_yyz_xyzz, g_y_0_yyz_xyzzz, g_y_0_yyz_xzzz, g_y_0_yyz_xzzzz, g_y_0_yyz_yyyy, g_y_0_yyz_yyyyz, g_y_0_yyz_yyyz, g_y_0_yyz_yyyzz, g_y_0_yyz_yyzz, g_y_0_yyz_yyzzz, g_y_0_yyz_yzzz, g_y_0_yyz_yzzzz, g_y_0_yyz_zzzz, g_y_0_yyz_zzzzz, g_y_0_yyzz_xxxx, g_y_0_yyzz_xxxy, g_y_0_yyzz_xxxz, g_y_0_yyzz_xxyy, g_y_0_yyzz_xxyz, g_y_0_yyzz_xxzz, g_y_0_yyzz_xyyy, g_y_0_yyzz_xyyz, g_y_0_yyzz_xyzz, g_y_0_yyzz_xzzz, g_y_0_yyzz_yyyy, g_y_0_yyzz_yyyz, g_y_0_yyzz_yyzz, g_y_0_yyzz_yzzz, g_y_0_yyzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyzz_xxxx[k] = -g_y_0_yyz_xxxx[k] * cd_z[k] + g_y_0_yyz_xxxxz[k];

            g_y_0_yyzz_xxxy[k] = -g_y_0_yyz_xxxy[k] * cd_z[k] + g_y_0_yyz_xxxyz[k];

            g_y_0_yyzz_xxxz[k] = -g_y_0_yyz_xxxz[k] * cd_z[k] + g_y_0_yyz_xxxzz[k];

            g_y_0_yyzz_xxyy[k] = -g_y_0_yyz_xxyy[k] * cd_z[k] + g_y_0_yyz_xxyyz[k];

            g_y_0_yyzz_xxyz[k] = -g_y_0_yyz_xxyz[k] * cd_z[k] + g_y_0_yyz_xxyzz[k];

            g_y_0_yyzz_xxzz[k] = -g_y_0_yyz_xxzz[k] * cd_z[k] + g_y_0_yyz_xxzzz[k];

            g_y_0_yyzz_xyyy[k] = -g_y_0_yyz_xyyy[k] * cd_z[k] + g_y_0_yyz_xyyyz[k];

            g_y_0_yyzz_xyyz[k] = -g_y_0_yyz_xyyz[k] * cd_z[k] + g_y_0_yyz_xyyzz[k];

            g_y_0_yyzz_xyzz[k] = -g_y_0_yyz_xyzz[k] * cd_z[k] + g_y_0_yyz_xyzzz[k];

            g_y_0_yyzz_xzzz[k] = -g_y_0_yyz_xzzz[k] * cd_z[k] + g_y_0_yyz_xzzzz[k];

            g_y_0_yyzz_yyyy[k] = -g_y_0_yyz_yyyy[k] * cd_z[k] + g_y_0_yyz_yyyyz[k];

            g_y_0_yyzz_yyyz[k] = -g_y_0_yyz_yyyz[k] * cd_z[k] + g_y_0_yyz_yyyzz[k];

            g_y_0_yyzz_yyzz[k] = -g_y_0_yyz_yyzz[k] * cd_z[k] + g_y_0_yyz_yyzzz[k];

            g_y_0_yyzz_yzzz[k] = -g_y_0_yyz_yzzz[k] * cd_z[k] + g_y_0_yyz_yzzzz[k];

            g_y_0_yyzz_zzzz[k] = -g_y_0_yyz_zzzz[k] * cd_z[k] + g_y_0_yyz_zzzzz[k];
        }

        /// Set up 195-210 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 195);

        auto g_y_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 196);

        auto g_y_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 197);

        auto g_y_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 198);

        auto g_y_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 199);

        auto g_y_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 200);

        auto g_y_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 201);

        auto g_y_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 202);

        auto g_y_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 203);

        auto g_y_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 204);

        auto g_y_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 205);

        auto g_y_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 206);

        auto g_y_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 207);

        auto g_y_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 208);

        auto g_y_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 209);

        #pragma omp simd aligned(cd_z, g_y_0_yzz_xxxx, g_y_0_yzz_xxxxz, g_y_0_yzz_xxxy, g_y_0_yzz_xxxyz, g_y_0_yzz_xxxz, g_y_0_yzz_xxxzz, g_y_0_yzz_xxyy, g_y_0_yzz_xxyyz, g_y_0_yzz_xxyz, g_y_0_yzz_xxyzz, g_y_0_yzz_xxzz, g_y_0_yzz_xxzzz, g_y_0_yzz_xyyy, g_y_0_yzz_xyyyz, g_y_0_yzz_xyyz, g_y_0_yzz_xyyzz, g_y_0_yzz_xyzz, g_y_0_yzz_xyzzz, g_y_0_yzz_xzzz, g_y_0_yzz_xzzzz, g_y_0_yzz_yyyy, g_y_0_yzz_yyyyz, g_y_0_yzz_yyyz, g_y_0_yzz_yyyzz, g_y_0_yzz_yyzz, g_y_0_yzz_yyzzz, g_y_0_yzz_yzzz, g_y_0_yzz_yzzzz, g_y_0_yzz_zzzz, g_y_0_yzz_zzzzz, g_y_0_yzzz_xxxx, g_y_0_yzzz_xxxy, g_y_0_yzzz_xxxz, g_y_0_yzzz_xxyy, g_y_0_yzzz_xxyz, g_y_0_yzzz_xxzz, g_y_0_yzzz_xyyy, g_y_0_yzzz_xyyz, g_y_0_yzzz_xyzz, g_y_0_yzzz_xzzz, g_y_0_yzzz_yyyy, g_y_0_yzzz_yyyz, g_y_0_yzzz_yyzz, g_y_0_yzzz_yzzz, g_y_0_yzzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzzz_xxxx[k] = -g_y_0_yzz_xxxx[k] * cd_z[k] + g_y_0_yzz_xxxxz[k];

            g_y_0_yzzz_xxxy[k] = -g_y_0_yzz_xxxy[k] * cd_z[k] + g_y_0_yzz_xxxyz[k];

            g_y_0_yzzz_xxxz[k] = -g_y_0_yzz_xxxz[k] * cd_z[k] + g_y_0_yzz_xxxzz[k];

            g_y_0_yzzz_xxyy[k] = -g_y_0_yzz_xxyy[k] * cd_z[k] + g_y_0_yzz_xxyyz[k];

            g_y_0_yzzz_xxyz[k] = -g_y_0_yzz_xxyz[k] * cd_z[k] + g_y_0_yzz_xxyzz[k];

            g_y_0_yzzz_xxzz[k] = -g_y_0_yzz_xxzz[k] * cd_z[k] + g_y_0_yzz_xxzzz[k];

            g_y_0_yzzz_xyyy[k] = -g_y_0_yzz_xyyy[k] * cd_z[k] + g_y_0_yzz_xyyyz[k];

            g_y_0_yzzz_xyyz[k] = -g_y_0_yzz_xyyz[k] * cd_z[k] + g_y_0_yzz_xyyzz[k];

            g_y_0_yzzz_xyzz[k] = -g_y_0_yzz_xyzz[k] * cd_z[k] + g_y_0_yzz_xyzzz[k];

            g_y_0_yzzz_xzzz[k] = -g_y_0_yzz_xzzz[k] * cd_z[k] + g_y_0_yzz_xzzzz[k];

            g_y_0_yzzz_yyyy[k] = -g_y_0_yzz_yyyy[k] * cd_z[k] + g_y_0_yzz_yyyyz[k];

            g_y_0_yzzz_yyyz[k] = -g_y_0_yzz_yyyz[k] * cd_z[k] + g_y_0_yzz_yyyzz[k];

            g_y_0_yzzz_yyzz[k] = -g_y_0_yzz_yyzz[k] * cd_z[k] + g_y_0_yzz_yyzzz[k];

            g_y_0_yzzz_yzzz[k] = -g_y_0_yzz_yzzz[k] * cd_z[k] + g_y_0_yzz_yzzzz[k];

            g_y_0_yzzz_zzzz[k] = -g_y_0_yzz_zzzz[k] * cd_z[k] + g_y_0_yzz_zzzzz[k];
        }

        /// Set up 210-225 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 225 * acomps  + 210);

        auto g_y_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 211);

        auto g_y_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 212);

        auto g_y_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 213);

        auto g_y_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 214);

        auto g_y_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 215);

        auto g_y_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 216);

        auto g_y_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 217);

        auto g_y_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 218);

        auto g_y_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 219);

        auto g_y_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 225 * acomps  + 220);

        auto g_y_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 221);

        auto g_y_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 222);

        auto g_y_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 223);

        auto g_y_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 225 * acomps  + 224);

        #pragma omp simd aligned(cd_z, g_y_0_zzz_xxxx, g_y_0_zzz_xxxxz, g_y_0_zzz_xxxy, g_y_0_zzz_xxxyz, g_y_0_zzz_xxxz, g_y_0_zzz_xxxzz, g_y_0_zzz_xxyy, g_y_0_zzz_xxyyz, g_y_0_zzz_xxyz, g_y_0_zzz_xxyzz, g_y_0_zzz_xxzz, g_y_0_zzz_xxzzz, g_y_0_zzz_xyyy, g_y_0_zzz_xyyyz, g_y_0_zzz_xyyz, g_y_0_zzz_xyyzz, g_y_0_zzz_xyzz, g_y_0_zzz_xyzzz, g_y_0_zzz_xzzz, g_y_0_zzz_xzzzz, g_y_0_zzz_yyyy, g_y_0_zzz_yyyyz, g_y_0_zzz_yyyz, g_y_0_zzz_yyyzz, g_y_0_zzz_yyzz, g_y_0_zzz_yyzzz, g_y_0_zzz_yzzz, g_y_0_zzz_yzzzz, g_y_0_zzz_zzzz, g_y_0_zzz_zzzzz, g_y_0_zzzz_xxxx, g_y_0_zzzz_xxxy, g_y_0_zzzz_xxxz, g_y_0_zzzz_xxyy, g_y_0_zzzz_xxyz, g_y_0_zzzz_xxzz, g_y_0_zzzz_xyyy, g_y_0_zzzz_xyyz, g_y_0_zzzz_xyzz, g_y_0_zzzz_xzzz, g_y_0_zzzz_yyyy, g_y_0_zzzz_yyyz, g_y_0_zzzz_yyzz, g_y_0_zzzz_yzzz, g_y_0_zzzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzzz_xxxx[k] = -g_y_0_zzz_xxxx[k] * cd_z[k] + g_y_0_zzz_xxxxz[k];

            g_y_0_zzzz_xxxy[k] = -g_y_0_zzz_xxxy[k] * cd_z[k] + g_y_0_zzz_xxxyz[k];

            g_y_0_zzzz_xxxz[k] = -g_y_0_zzz_xxxz[k] * cd_z[k] + g_y_0_zzz_xxxzz[k];

            g_y_0_zzzz_xxyy[k] = -g_y_0_zzz_xxyy[k] * cd_z[k] + g_y_0_zzz_xxyyz[k];

            g_y_0_zzzz_xxyz[k] = -g_y_0_zzz_xxyz[k] * cd_z[k] + g_y_0_zzz_xxyzz[k];

            g_y_0_zzzz_xxzz[k] = -g_y_0_zzz_xxzz[k] * cd_z[k] + g_y_0_zzz_xxzzz[k];

            g_y_0_zzzz_xyyy[k] = -g_y_0_zzz_xyyy[k] * cd_z[k] + g_y_0_zzz_xyyyz[k];

            g_y_0_zzzz_xyyz[k] = -g_y_0_zzz_xyyz[k] * cd_z[k] + g_y_0_zzz_xyyzz[k];

            g_y_0_zzzz_xyzz[k] = -g_y_0_zzz_xyzz[k] * cd_z[k] + g_y_0_zzz_xyzzz[k];

            g_y_0_zzzz_xzzz[k] = -g_y_0_zzz_xzzz[k] * cd_z[k] + g_y_0_zzz_xzzzz[k];

            g_y_0_zzzz_yyyy[k] = -g_y_0_zzz_yyyy[k] * cd_z[k] + g_y_0_zzz_yyyyz[k];

            g_y_0_zzzz_yyyz[k] = -g_y_0_zzz_yyyz[k] * cd_z[k] + g_y_0_zzz_yyyzz[k];

            g_y_0_zzzz_yyzz[k] = -g_y_0_zzz_yyzz[k] * cd_z[k] + g_y_0_zzz_yyzzz[k];

            g_y_0_zzzz_yzzz[k] = -g_y_0_zzz_yzzz[k] * cd_z[k] + g_y_0_zzz_yzzzz[k];

            g_y_0_zzzz_zzzz[k] = -g_y_0_zzz_zzzz[k] * cd_z[k] + g_y_0_zzz_zzzzz[k];
        }
        /// Set up 0-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxx_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 0);

        auto g_z_0_xxxx_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 1);

        auto g_z_0_xxxx_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 2);

        auto g_z_0_xxxx_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 3);

        auto g_z_0_xxxx_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 4);

        auto g_z_0_xxxx_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 5);

        auto g_z_0_xxxx_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 6);

        auto g_z_0_xxxx_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 7);

        auto g_z_0_xxxx_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 8);

        auto g_z_0_xxxx_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 9);

        auto g_z_0_xxxx_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 10);

        auto g_z_0_xxxx_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 11);

        auto g_z_0_xxxx_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 12);

        auto g_z_0_xxxx_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 13);

        auto g_z_0_xxxx_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_z_0_xxx_xxxx, g_z_0_xxx_xxxxx, g_z_0_xxx_xxxxy, g_z_0_xxx_xxxxz, g_z_0_xxx_xxxy, g_z_0_xxx_xxxyy, g_z_0_xxx_xxxyz, g_z_0_xxx_xxxz, g_z_0_xxx_xxxzz, g_z_0_xxx_xxyy, g_z_0_xxx_xxyyy, g_z_0_xxx_xxyyz, g_z_0_xxx_xxyz, g_z_0_xxx_xxyzz, g_z_0_xxx_xxzz, g_z_0_xxx_xxzzz, g_z_0_xxx_xyyy, g_z_0_xxx_xyyyy, g_z_0_xxx_xyyyz, g_z_0_xxx_xyyz, g_z_0_xxx_xyyzz, g_z_0_xxx_xyzz, g_z_0_xxx_xyzzz, g_z_0_xxx_xzzz, g_z_0_xxx_xzzzz, g_z_0_xxx_yyyy, g_z_0_xxx_yyyz, g_z_0_xxx_yyzz, g_z_0_xxx_yzzz, g_z_0_xxx_zzzz, g_z_0_xxxx_xxxx, g_z_0_xxxx_xxxy, g_z_0_xxxx_xxxz, g_z_0_xxxx_xxyy, g_z_0_xxxx_xxyz, g_z_0_xxxx_xxzz, g_z_0_xxxx_xyyy, g_z_0_xxxx_xyyz, g_z_0_xxxx_xyzz, g_z_0_xxxx_xzzz, g_z_0_xxxx_yyyy, g_z_0_xxxx_yyyz, g_z_0_xxxx_yyzz, g_z_0_xxxx_yzzz, g_z_0_xxxx_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxx_xxxx[k] = -g_z_0_xxx_xxxx[k] * cd_x[k] + g_z_0_xxx_xxxxx[k];

            g_z_0_xxxx_xxxy[k] = -g_z_0_xxx_xxxy[k] * cd_x[k] + g_z_0_xxx_xxxxy[k];

            g_z_0_xxxx_xxxz[k] = -g_z_0_xxx_xxxz[k] * cd_x[k] + g_z_0_xxx_xxxxz[k];

            g_z_0_xxxx_xxyy[k] = -g_z_0_xxx_xxyy[k] * cd_x[k] + g_z_0_xxx_xxxyy[k];

            g_z_0_xxxx_xxyz[k] = -g_z_0_xxx_xxyz[k] * cd_x[k] + g_z_0_xxx_xxxyz[k];

            g_z_0_xxxx_xxzz[k] = -g_z_0_xxx_xxzz[k] * cd_x[k] + g_z_0_xxx_xxxzz[k];

            g_z_0_xxxx_xyyy[k] = -g_z_0_xxx_xyyy[k] * cd_x[k] + g_z_0_xxx_xxyyy[k];

            g_z_0_xxxx_xyyz[k] = -g_z_0_xxx_xyyz[k] * cd_x[k] + g_z_0_xxx_xxyyz[k];

            g_z_0_xxxx_xyzz[k] = -g_z_0_xxx_xyzz[k] * cd_x[k] + g_z_0_xxx_xxyzz[k];

            g_z_0_xxxx_xzzz[k] = -g_z_0_xxx_xzzz[k] * cd_x[k] + g_z_0_xxx_xxzzz[k];

            g_z_0_xxxx_yyyy[k] = -g_z_0_xxx_yyyy[k] * cd_x[k] + g_z_0_xxx_xyyyy[k];

            g_z_0_xxxx_yyyz[k] = -g_z_0_xxx_yyyz[k] * cd_x[k] + g_z_0_xxx_xyyyz[k];

            g_z_0_xxxx_yyzz[k] = -g_z_0_xxx_yyzz[k] * cd_x[k] + g_z_0_xxx_xyyzz[k];

            g_z_0_xxxx_yzzz[k] = -g_z_0_xxx_yzzz[k] * cd_x[k] + g_z_0_xxx_xyzzz[k];

            g_z_0_xxxx_zzzz[k] = -g_z_0_xxx_zzzz[k] * cd_x[k] + g_z_0_xxx_xzzzz[k];
        }

        /// Set up 15-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 15);

        auto g_z_0_xxxy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 16);

        auto g_z_0_xxxy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 17);

        auto g_z_0_xxxy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 18);

        auto g_z_0_xxxy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 19);

        auto g_z_0_xxxy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 20);

        auto g_z_0_xxxy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 21);

        auto g_z_0_xxxy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 22);

        auto g_z_0_xxxy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 23);

        auto g_z_0_xxxy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 24);

        auto g_z_0_xxxy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 25);

        auto g_z_0_xxxy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 26);

        auto g_z_0_xxxy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 27);

        auto g_z_0_xxxy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 28);

        auto g_z_0_xxxy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xxxy_xxxx, g_z_0_xxxy_xxxy, g_z_0_xxxy_xxxz, g_z_0_xxxy_xxyy, g_z_0_xxxy_xxyz, g_z_0_xxxy_xxzz, g_z_0_xxxy_xyyy, g_z_0_xxxy_xyyz, g_z_0_xxxy_xyzz, g_z_0_xxxy_xzzz, g_z_0_xxxy_yyyy, g_z_0_xxxy_yyyz, g_z_0_xxxy_yyzz, g_z_0_xxxy_yzzz, g_z_0_xxxy_zzzz, g_z_0_xxy_xxxx, g_z_0_xxy_xxxxx, g_z_0_xxy_xxxxy, g_z_0_xxy_xxxxz, g_z_0_xxy_xxxy, g_z_0_xxy_xxxyy, g_z_0_xxy_xxxyz, g_z_0_xxy_xxxz, g_z_0_xxy_xxxzz, g_z_0_xxy_xxyy, g_z_0_xxy_xxyyy, g_z_0_xxy_xxyyz, g_z_0_xxy_xxyz, g_z_0_xxy_xxyzz, g_z_0_xxy_xxzz, g_z_0_xxy_xxzzz, g_z_0_xxy_xyyy, g_z_0_xxy_xyyyy, g_z_0_xxy_xyyyz, g_z_0_xxy_xyyz, g_z_0_xxy_xyyzz, g_z_0_xxy_xyzz, g_z_0_xxy_xyzzz, g_z_0_xxy_xzzz, g_z_0_xxy_xzzzz, g_z_0_xxy_yyyy, g_z_0_xxy_yyyz, g_z_0_xxy_yyzz, g_z_0_xxy_yzzz, g_z_0_xxy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxy_xxxx[k] = -g_z_0_xxy_xxxx[k] * cd_x[k] + g_z_0_xxy_xxxxx[k];

            g_z_0_xxxy_xxxy[k] = -g_z_0_xxy_xxxy[k] * cd_x[k] + g_z_0_xxy_xxxxy[k];

            g_z_0_xxxy_xxxz[k] = -g_z_0_xxy_xxxz[k] * cd_x[k] + g_z_0_xxy_xxxxz[k];

            g_z_0_xxxy_xxyy[k] = -g_z_0_xxy_xxyy[k] * cd_x[k] + g_z_0_xxy_xxxyy[k];

            g_z_0_xxxy_xxyz[k] = -g_z_0_xxy_xxyz[k] * cd_x[k] + g_z_0_xxy_xxxyz[k];

            g_z_0_xxxy_xxzz[k] = -g_z_0_xxy_xxzz[k] * cd_x[k] + g_z_0_xxy_xxxzz[k];

            g_z_0_xxxy_xyyy[k] = -g_z_0_xxy_xyyy[k] * cd_x[k] + g_z_0_xxy_xxyyy[k];

            g_z_0_xxxy_xyyz[k] = -g_z_0_xxy_xyyz[k] * cd_x[k] + g_z_0_xxy_xxyyz[k];

            g_z_0_xxxy_xyzz[k] = -g_z_0_xxy_xyzz[k] * cd_x[k] + g_z_0_xxy_xxyzz[k];

            g_z_0_xxxy_xzzz[k] = -g_z_0_xxy_xzzz[k] * cd_x[k] + g_z_0_xxy_xxzzz[k];

            g_z_0_xxxy_yyyy[k] = -g_z_0_xxy_yyyy[k] * cd_x[k] + g_z_0_xxy_xyyyy[k];

            g_z_0_xxxy_yyyz[k] = -g_z_0_xxy_yyyz[k] * cd_x[k] + g_z_0_xxy_xyyyz[k];

            g_z_0_xxxy_yyzz[k] = -g_z_0_xxy_yyzz[k] * cd_x[k] + g_z_0_xxy_xyyzz[k];

            g_z_0_xxxy_yzzz[k] = -g_z_0_xxy_yzzz[k] * cd_x[k] + g_z_0_xxy_xyzzz[k];

            g_z_0_xxxy_zzzz[k] = -g_z_0_xxy_zzzz[k] * cd_x[k] + g_z_0_xxy_xzzzz[k];
        }

        /// Set up 30-45 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxxz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 30);

        auto g_z_0_xxxz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 31);

        auto g_z_0_xxxz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 32);

        auto g_z_0_xxxz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 33);

        auto g_z_0_xxxz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 34);

        auto g_z_0_xxxz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 35);

        auto g_z_0_xxxz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 36);

        auto g_z_0_xxxz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 37);

        auto g_z_0_xxxz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 38);

        auto g_z_0_xxxz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 39);

        auto g_z_0_xxxz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 40);

        auto g_z_0_xxxz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 41);

        auto g_z_0_xxxz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 42);

        auto g_z_0_xxxz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 43);

        auto g_z_0_xxxz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 44);

        #pragma omp simd aligned(cd_x, g_z_0_xxxz_xxxx, g_z_0_xxxz_xxxy, g_z_0_xxxz_xxxz, g_z_0_xxxz_xxyy, g_z_0_xxxz_xxyz, g_z_0_xxxz_xxzz, g_z_0_xxxz_xyyy, g_z_0_xxxz_xyyz, g_z_0_xxxz_xyzz, g_z_0_xxxz_xzzz, g_z_0_xxxz_yyyy, g_z_0_xxxz_yyyz, g_z_0_xxxz_yyzz, g_z_0_xxxz_yzzz, g_z_0_xxxz_zzzz, g_z_0_xxz_xxxx, g_z_0_xxz_xxxxx, g_z_0_xxz_xxxxy, g_z_0_xxz_xxxxz, g_z_0_xxz_xxxy, g_z_0_xxz_xxxyy, g_z_0_xxz_xxxyz, g_z_0_xxz_xxxz, g_z_0_xxz_xxxzz, g_z_0_xxz_xxyy, g_z_0_xxz_xxyyy, g_z_0_xxz_xxyyz, g_z_0_xxz_xxyz, g_z_0_xxz_xxyzz, g_z_0_xxz_xxzz, g_z_0_xxz_xxzzz, g_z_0_xxz_xyyy, g_z_0_xxz_xyyyy, g_z_0_xxz_xyyyz, g_z_0_xxz_xyyz, g_z_0_xxz_xyyzz, g_z_0_xxz_xyzz, g_z_0_xxz_xyzzz, g_z_0_xxz_xzzz, g_z_0_xxz_xzzzz, g_z_0_xxz_yyyy, g_z_0_xxz_yyyz, g_z_0_xxz_yyzz, g_z_0_xxz_yzzz, g_z_0_xxz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxxz_xxxx[k] = -g_z_0_xxz_xxxx[k] * cd_x[k] + g_z_0_xxz_xxxxx[k];

            g_z_0_xxxz_xxxy[k] = -g_z_0_xxz_xxxy[k] * cd_x[k] + g_z_0_xxz_xxxxy[k];

            g_z_0_xxxz_xxxz[k] = -g_z_0_xxz_xxxz[k] * cd_x[k] + g_z_0_xxz_xxxxz[k];

            g_z_0_xxxz_xxyy[k] = -g_z_0_xxz_xxyy[k] * cd_x[k] + g_z_0_xxz_xxxyy[k];

            g_z_0_xxxz_xxyz[k] = -g_z_0_xxz_xxyz[k] * cd_x[k] + g_z_0_xxz_xxxyz[k];

            g_z_0_xxxz_xxzz[k] = -g_z_0_xxz_xxzz[k] * cd_x[k] + g_z_0_xxz_xxxzz[k];

            g_z_0_xxxz_xyyy[k] = -g_z_0_xxz_xyyy[k] * cd_x[k] + g_z_0_xxz_xxyyy[k];

            g_z_0_xxxz_xyyz[k] = -g_z_0_xxz_xyyz[k] * cd_x[k] + g_z_0_xxz_xxyyz[k];

            g_z_0_xxxz_xyzz[k] = -g_z_0_xxz_xyzz[k] * cd_x[k] + g_z_0_xxz_xxyzz[k];

            g_z_0_xxxz_xzzz[k] = -g_z_0_xxz_xzzz[k] * cd_x[k] + g_z_0_xxz_xxzzz[k];

            g_z_0_xxxz_yyyy[k] = -g_z_0_xxz_yyyy[k] * cd_x[k] + g_z_0_xxz_xyyyy[k];

            g_z_0_xxxz_yyyz[k] = -g_z_0_xxz_yyyz[k] * cd_x[k] + g_z_0_xxz_xyyyz[k];

            g_z_0_xxxz_yyzz[k] = -g_z_0_xxz_yyzz[k] * cd_x[k] + g_z_0_xxz_xyyzz[k];

            g_z_0_xxxz_yzzz[k] = -g_z_0_xxz_yzzz[k] * cd_x[k] + g_z_0_xxz_xyzzz[k];

            g_z_0_xxxz_zzzz[k] = -g_z_0_xxz_zzzz[k] * cd_x[k] + g_z_0_xxz_xzzzz[k];
        }

        /// Set up 45-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 45);

        auto g_z_0_xxyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 46);

        auto g_z_0_xxyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 47);

        auto g_z_0_xxyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 48);

        auto g_z_0_xxyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 49);

        auto g_z_0_xxyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 50);

        auto g_z_0_xxyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 51);

        auto g_z_0_xxyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 52);

        auto g_z_0_xxyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 53);

        auto g_z_0_xxyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 54);

        auto g_z_0_xxyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 55);

        auto g_z_0_xxyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 56);

        auto g_z_0_xxyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 57);

        auto g_z_0_xxyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 58);

        auto g_z_0_xxyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 59);

        #pragma omp simd aligned(cd_x, g_z_0_xxyy_xxxx, g_z_0_xxyy_xxxy, g_z_0_xxyy_xxxz, g_z_0_xxyy_xxyy, g_z_0_xxyy_xxyz, g_z_0_xxyy_xxzz, g_z_0_xxyy_xyyy, g_z_0_xxyy_xyyz, g_z_0_xxyy_xyzz, g_z_0_xxyy_xzzz, g_z_0_xxyy_yyyy, g_z_0_xxyy_yyyz, g_z_0_xxyy_yyzz, g_z_0_xxyy_yzzz, g_z_0_xxyy_zzzz, g_z_0_xyy_xxxx, g_z_0_xyy_xxxxx, g_z_0_xyy_xxxxy, g_z_0_xyy_xxxxz, g_z_0_xyy_xxxy, g_z_0_xyy_xxxyy, g_z_0_xyy_xxxyz, g_z_0_xyy_xxxz, g_z_0_xyy_xxxzz, g_z_0_xyy_xxyy, g_z_0_xyy_xxyyy, g_z_0_xyy_xxyyz, g_z_0_xyy_xxyz, g_z_0_xyy_xxyzz, g_z_0_xyy_xxzz, g_z_0_xyy_xxzzz, g_z_0_xyy_xyyy, g_z_0_xyy_xyyyy, g_z_0_xyy_xyyyz, g_z_0_xyy_xyyz, g_z_0_xyy_xyyzz, g_z_0_xyy_xyzz, g_z_0_xyy_xyzzz, g_z_0_xyy_xzzz, g_z_0_xyy_xzzzz, g_z_0_xyy_yyyy, g_z_0_xyy_yyyz, g_z_0_xyy_yyzz, g_z_0_xyy_yzzz, g_z_0_xyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyy_xxxx[k] = -g_z_0_xyy_xxxx[k] * cd_x[k] + g_z_0_xyy_xxxxx[k];

            g_z_0_xxyy_xxxy[k] = -g_z_0_xyy_xxxy[k] * cd_x[k] + g_z_0_xyy_xxxxy[k];

            g_z_0_xxyy_xxxz[k] = -g_z_0_xyy_xxxz[k] * cd_x[k] + g_z_0_xyy_xxxxz[k];

            g_z_0_xxyy_xxyy[k] = -g_z_0_xyy_xxyy[k] * cd_x[k] + g_z_0_xyy_xxxyy[k];

            g_z_0_xxyy_xxyz[k] = -g_z_0_xyy_xxyz[k] * cd_x[k] + g_z_0_xyy_xxxyz[k];

            g_z_0_xxyy_xxzz[k] = -g_z_0_xyy_xxzz[k] * cd_x[k] + g_z_0_xyy_xxxzz[k];

            g_z_0_xxyy_xyyy[k] = -g_z_0_xyy_xyyy[k] * cd_x[k] + g_z_0_xyy_xxyyy[k];

            g_z_0_xxyy_xyyz[k] = -g_z_0_xyy_xyyz[k] * cd_x[k] + g_z_0_xyy_xxyyz[k];

            g_z_0_xxyy_xyzz[k] = -g_z_0_xyy_xyzz[k] * cd_x[k] + g_z_0_xyy_xxyzz[k];

            g_z_0_xxyy_xzzz[k] = -g_z_0_xyy_xzzz[k] * cd_x[k] + g_z_0_xyy_xxzzz[k];

            g_z_0_xxyy_yyyy[k] = -g_z_0_xyy_yyyy[k] * cd_x[k] + g_z_0_xyy_xyyyy[k];

            g_z_0_xxyy_yyyz[k] = -g_z_0_xyy_yyyz[k] * cd_x[k] + g_z_0_xyy_xyyyz[k];

            g_z_0_xxyy_yyzz[k] = -g_z_0_xyy_yyzz[k] * cd_x[k] + g_z_0_xyy_xyyzz[k];

            g_z_0_xxyy_yzzz[k] = -g_z_0_xyy_yzzz[k] * cd_x[k] + g_z_0_xyy_xyzzz[k];

            g_z_0_xxyy_zzzz[k] = -g_z_0_xyy_zzzz[k] * cd_x[k] + g_z_0_xyy_xzzzz[k];
        }

        /// Set up 60-75 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 60);

        auto g_z_0_xxyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 61);

        auto g_z_0_xxyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 62);

        auto g_z_0_xxyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 63);

        auto g_z_0_xxyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 64);

        auto g_z_0_xxyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 65);

        auto g_z_0_xxyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 66);

        auto g_z_0_xxyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 67);

        auto g_z_0_xxyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 68);

        auto g_z_0_xxyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 69);

        auto g_z_0_xxyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 70);

        auto g_z_0_xxyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 71);

        auto g_z_0_xxyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 72);

        auto g_z_0_xxyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 73);

        auto g_z_0_xxyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 74);

        #pragma omp simd aligned(cd_x, g_z_0_xxyz_xxxx, g_z_0_xxyz_xxxy, g_z_0_xxyz_xxxz, g_z_0_xxyz_xxyy, g_z_0_xxyz_xxyz, g_z_0_xxyz_xxzz, g_z_0_xxyz_xyyy, g_z_0_xxyz_xyyz, g_z_0_xxyz_xyzz, g_z_0_xxyz_xzzz, g_z_0_xxyz_yyyy, g_z_0_xxyz_yyyz, g_z_0_xxyz_yyzz, g_z_0_xxyz_yzzz, g_z_0_xxyz_zzzz, g_z_0_xyz_xxxx, g_z_0_xyz_xxxxx, g_z_0_xyz_xxxxy, g_z_0_xyz_xxxxz, g_z_0_xyz_xxxy, g_z_0_xyz_xxxyy, g_z_0_xyz_xxxyz, g_z_0_xyz_xxxz, g_z_0_xyz_xxxzz, g_z_0_xyz_xxyy, g_z_0_xyz_xxyyy, g_z_0_xyz_xxyyz, g_z_0_xyz_xxyz, g_z_0_xyz_xxyzz, g_z_0_xyz_xxzz, g_z_0_xyz_xxzzz, g_z_0_xyz_xyyy, g_z_0_xyz_xyyyy, g_z_0_xyz_xyyyz, g_z_0_xyz_xyyz, g_z_0_xyz_xyyzz, g_z_0_xyz_xyzz, g_z_0_xyz_xyzzz, g_z_0_xyz_xzzz, g_z_0_xyz_xzzzz, g_z_0_xyz_yyyy, g_z_0_xyz_yyyz, g_z_0_xyz_yyzz, g_z_0_xyz_yzzz, g_z_0_xyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxyz_xxxx[k] = -g_z_0_xyz_xxxx[k] * cd_x[k] + g_z_0_xyz_xxxxx[k];

            g_z_0_xxyz_xxxy[k] = -g_z_0_xyz_xxxy[k] * cd_x[k] + g_z_0_xyz_xxxxy[k];

            g_z_0_xxyz_xxxz[k] = -g_z_0_xyz_xxxz[k] * cd_x[k] + g_z_0_xyz_xxxxz[k];

            g_z_0_xxyz_xxyy[k] = -g_z_0_xyz_xxyy[k] * cd_x[k] + g_z_0_xyz_xxxyy[k];

            g_z_0_xxyz_xxyz[k] = -g_z_0_xyz_xxyz[k] * cd_x[k] + g_z_0_xyz_xxxyz[k];

            g_z_0_xxyz_xxzz[k] = -g_z_0_xyz_xxzz[k] * cd_x[k] + g_z_0_xyz_xxxzz[k];

            g_z_0_xxyz_xyyy[k] = -g_z_0_xyz_xyyy[k] * cd_x[k] + g_z_0_xyz_xxyyy[k];

            g_z_0_xxyz_xyyz[k] = -g_z_0_xyz_xyyz[k] * cd_x[k] + g_z_0_xyz_xxyyz[k];

            g_z_0_xxyz_xyzz[k] = -g_z_0_xyz_xyzz[k] * cd_x[k] + g_z_0_xyz_xxyzz[k];

            g_z_0_xxyz_xzzz[k] = -g_z_0_xyz_xzzz[k] * cd_x[k] + g_z_0_xyz_xxzzz[k];

            g_z_0_xxyz_yyyy[k] = -g_z_0_xyz_yyyy[k] * cd_x[k] + g_z_0_xyz_xyyyy[k];

            g_z_0_xxyz_yyyz[k] = -g_z_0_xyz_yyyz[k] * cd_x[k] + g_z_0_xyz_xyyyz[k];

            g_z_0_xxyz_yyzz[k] = -g_z_0_xyz_yyzz[k] * cd_x[k] + g_z_0_xyz_xyyzz[k];

            g_z_0_xxyz_yzzz[k] = -g_z_0_xyz_yzzz[k] * cd_x[k] + g_z_0_xyz_xyzzz[k];

            g_z_0_xxyz_zzzz[k] = -g_z_0_xyz_zzzz[k] * cd_x[k] + g_z_0_xyz_xzzzz[k];
        }

        /// Set up 75-90 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 75);

        auto g_z_0_xxzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 76);

        auto g_z_0_xxzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 77);

        auto g_z_0_xxzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 78);

        auto g_z_0_xxzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 79);

        auto g_z_0_xxzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 80);

        auto g_z_0_xxzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 81);

        auto g_z_0_xxzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 82);

        auto g_z_0_xxzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 83);

        auto g_z_0_xxzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 84);

        auto g_z_0_xxzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 85);

        auto g_z_0_xxzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 86);

        auto g_z_0_xxzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 87);

        auto g_z_0_xxzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 88);

        auto g_z_0_xxzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 89);

        #pragma omp simd aligned(cd_x, g_z_0_xxzz_xxxx, g_z_0_xxzz_xxxy, g_z_0_xxzz_xxxz, g_z_0_xxzz_xxyy, g_z_0_xxzz_xxyz, g_z_0_xxzz_xxzz, g_z_0_xxzz_xyyy, g_z_0_xxzz_xyyz, g_z_0_xxzz_xyzz, g_z_0_xxzz_xzzz, g_z_0_xxzz_yyyy, g_z_0_xxzz_yyyz, g_z_0_xxzz_yyzz, g_z_0_xxzz_yzzz, g_z_0_xxzz_zzzz, g_z_0_xzz_xxxx, g_z_0_xzz_xxxxx, g_z_0_xzz_xxxxy, g_z_0_xzz_xxxxz, g_z_0_xzz_xxxy, g_z_0_xzz_xxxyy, g_z_0_xzz_xxxyz, g_z_0_xzz_xxxz, g_z_0_xzz_xxxzz, g_z_0_xzz_xxyy, g_z_0_xzz_xxyyy, g_z_0_xzz_xxyyz, g_z_0_xzz_xxyz, g_z_0_xzz_xxyzz, g_z_0_xzz_xxzz, g_z_0_xzz_xxzzz, g_z_0_xzz_xyyy, g_z_0_xzz_xyyyy, g_z_0_xzz_xyyyz, g_z_0_xzz_xyyz, g_z_0_xzz_xyyzz, g_z_0_xzz_xyzz, g_z_0_xzz_xyzzz, g_z_0_xzz_xzzz, g_z_0_xzz_xzzzz, g_z_0_xzz_yyyy, g_z_0_xzz_yyyz, g_z_0_xzz_yyzz, g_z_0_xzz_yzzz, g_z_0_xzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxzz_xxxx[k] = -g_z_0_xzz_xxxx[k] * cd_x[k] + g_z_0_xzz_xxxxx[k];

            g_z_0_xxzz_xxxy[k] = -g_z_0_xzz_xxxy[k] * cd_x[k] + g_z_0_xzz_xxxxy[k];

            g_z_0_xxzz_xxxz[k] = -g_z_0_xzz_xxxz[k] * cd_x[k] + g_z_0_xzz_xxxxz[k];

            g_z_0_xxzz_xxyy[k] = -g_z_0_xzz_xxyy[k] * cd_x[k] + g_z_0_xzz_xxxyy[k];

            g_z_0_xxzz_xxyz[k] = -g_z_0_xzz_xxyz[k] * cd_x[k] + g_z_0_xzz_xxxyz[k];

            g_z_0_xxzz_xxzz[k] = -g_z_0_xzz_xxzz[k] * cd_x[k] + g_z_0_xzz_xxxzz[k];

            g_z_0_xxzz_xyyy[k] = -g_z_0_xzz_xyyy[k] * cd_x[k] + g_z_0_xzz_xxyyy[k];

            g_z_0_xxzz_xyyz[k] = -g_z_0_xzz_xyyz[k] * cd_x[k] + g_z_0_xzz_xxyyz[k];

            g_z_0_xxzz_xyzz[k] = -g_z_0_xzz_xyzz[k] * cd_x[k] + g_z_0_xzz_xxyzz[k];

            g_z_0_xxzz_xzzz[k] = -g_z_0_xzz_xzzz[k] * cd_x[k] + g_z_0_xzz_xxzzz[k];

            g_z_0_xxzz_yyyy[k] = -g_z_0_xzz_yyyy[k] * cd_x[k] + g_z_0_xzz_xyyyy[k];

            g_z_0_xxzz_yyyz[k] = -g_z_0_xzz_yyyz[k] * cd_x[k] + g_z_0_xzz_xyyyz[k];

            g_z_0_xxzz_yyzz[k] = -g_z_0_xzz_yyzz[k] * cd_x[k] + g_z_0_xzz_xyyzz[k];

            g_z_0_xxzz_yzzz[k] = -g_z_0_xzz_yzzz[k] * cd_x[k] + g_z_0_xzz_xyzzz[k];

            g_z_0_xxzz_zzzz[k] = -g_z_0_xzz_zzzz[k] * cd_x[k] + g_z_0_xzz_xzzzz[k];
        }

        /// Set up 90-105 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 90);

        auto g_z_0_xyyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 91);

        auto g_z_0_xyyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 92);

        auto g_z_0_xyyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 93);

        auto g_z_0_xyyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 94);

        auto g_z_0_xyyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 95);

        auto g_z_0_xyyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 96);

        auto g_z_0_xyyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 97);

        auto g_z_0_xyyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 98);

        auto g_z_0_xyyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 99);

        auto g_z_0_xyyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 100);

        auto g_z_0_xyyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 101);

        auto g_z_0_xyyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 102);

        auto g_z_0_xyyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 103);

        auto g_z_0_xyyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 104);

        #pragma omp simd aligned(cd_x, g_z_0_xyyy_xxxx, g_z_0_xyyy_xxxy, g_z_0_xyyy_xxxz, g_z_0_xyyy_xxyy, g_z_0_xyyy_xxyz, g_z_0_xyyy_xxzz, g_z_0_xyyy_xyyy, g_z_0_xyyy_xyyz, g_z_0_xyyy_xyzz, g_z_0_xyyy_xzzz, g_z_0_xyyy_yyyy, g_z_0_xyyy_yyyz, g_z_0_xyyy_yyzz, g_z_0_xyyy_yzzz, g_z_0_xyyy_zzzz, g_z_0_yyy_xxxx, g_z_0_yyy_xxxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxxz, g_z_0_yyy_xxxy, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxz, g_z_0_yyy_xxxzz, g_z_0_yyy_xxyy, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzz, g_z_0_yyy_xxzzz, g_z_0_yyy_xyyy, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzz, g_z_0_yyy_xzzzz, g_z_0_yyy_yyyy, g_z_0_yyy_yyyz, g_z_0_yyy_yyzz, g_z_0_yyy_yzzz, g_z_0_yyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyy_xxxx[k] = -g_z_0_yyy_xxxx[k] * cd_x[k] + g_z_0_yyy_xxxxx[k];

            g_z_0_xyyy_xxxy[k] = -g_z_0_yyy_xxxy[k] * cd_x[k] + g_z_0_yyy_xxxxy[k];

            g_z_0_xyyy_xxxz[k] = -g_z_0_yyy_xxxz[k] * cd_x[k] + g_z_0_yyy_xxxxz[k];

            g_z_0_xyyy_xxyy[k] = -g_z_0_yyy_xxyy[k] * cd_x[k] + g_z_0_yyy_xxxyy[k];

            g_z_0_xyyy_xxyz[k] = -g_z_0_yyy_xxyz[k] * cd_x[k] + g_z_0_yyy_xxxyz[k];

            g_z_0_xyyy_xxzz[k] = -g_z_0_yyy_xxzz[k] * cd_x[k] + g_z_0_yyy_xxxzz[k];

            g_z_0_xyyy_xyyy[k] = -g_z_0_yyy_xyyy[k] * cd_x[k] + g_z_0_yyy_xxyyy[k];

            g_z_0_xyyy_xyyz[k] = -g_z_0_yyy_xyyz[k] * cd_x[k] + g_z_0_yyy_xxyyz[k];

            g_z_0_xyyy_xyzz[k] = -g_z_0_yyy_xyzz[k] * cd_x[k] + g_z_0_yyy_xxyzz[k];

            g_z_0_xyyy_xzzz[k] = -g_z_0_yyy_xzzz[k] * cd_x[k] + g_z_0_yyy_xxzzz[k];

            g_z_0_xyyy_yyyy[k] = -g_z_0_yyy_yyyy[k] * cd_x[k] + g_z_0_yyy_xyyyy[k];

            g_z_0_xyyy_yyyz[k] = -g_z_0_yyy_yyyz[k] * cd_x[k] + g_z_0_yyy_xyyyz[k];

            g_z_0_xyyy_yyzz[k] = -g_z_0_yyy_yyzz[k] * cd_x[k] + g_z_0_yyy_xyyzz[k];

            g_z_0_xyyy_yzzz[k] = -g_z_0_yyy_yzzz[k] * cd_x[k] + g_z_0_yyy_xyzzz[k];

            g_z_0_xyyy_zzzz[k] = -g_z_0_yyy_zzzz[k] * cd_x[k] + g_z_0_yyy_xzzzz[k];
        }

        /// Set up 105-120 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 105);

        auto g_z_0_xyyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 106);

        auto g_z_0_xyyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 107);

        auto g_z_0_xyyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 108);

        auto g_z_0_xyyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 109);

        auto g_z_0_xyyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 110);

        auto g_z_0_xyyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 111);

        auto g_z_0_xyyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 112);

        auto g_z_0_xyyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 113);

        auto g_z_0_xyyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 114);

        auto g_z_0_xyyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 115);

        auto g_z_0_xyyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 116);

        auto g_z_0_xyyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 117);

        auto g_z_0_xyyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 118);

        auto g_z_0_xyyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 119);

        #pragma omp simd aligned(cd_x, g_z_0_xyyz_xxxx, g_z_0_xyyz_xxxy, g_z_0_xyyz_xxxz, g_z_0_xyyz_xxyy, g_z_0_xyyz_xxyz, g_z_0_xyyz_xxzz, g_z_0_xyyz_xyyy, g_z_0_xyyz_xyyz, g_z_0_xyyz_xyzz, g_z_0_xyyz_xzzz, g_z_0_xyyz_yyyy, g_z_0_xyyz_yyyz, g_z_0_xyyz_yyzz, g_z_0_xyyz_yzzz, g_z_0_xyyz_zzzz, g_z_0_yyz_xxxx, g_z_0_yyz_xxxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxxz, g_z_0_yyz_xxxy, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxz, g_z_0_yyz_xxxzz, g_z_0_yyz_xxyy, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzz, g_z_0_yyz_xxzzz, g_z_0_yyz_xyyy, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzz, g_z_0_yyz_xzzzz, g_z_0_yyz_yyyy, g_z_0_yyz_yyyz, g_z_0_yyz_yyzz, g_z_0_yyz_yzzz, g_z_0_yyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyyz_xxxx[k] = -g_z_0_yyz_xxxx[k] * cd_x[k] + g_z_0_yyz_xxxxx[k];

            g_z_0_xyyz_xxxy[k] = -g_z_0_yyz_xxxy[k] * cd_x[k] + g_z_0_yyz_xxxxy[k];

            g_z_0_xyyz_xxxz[k] = -g_z_0_yyz_xxxz[k] * cd_x[k] + g_z_0_yyz_xxxxz[k];

            g_z_0_xyyz_xxyy[k] = -g_z_0_yyz_xxyy[k] * cd_x[k] + g_z_0_yyz_xxxyy[k];

            g_z_0_xyyz_xxyz[k] = -g_z_0_yyz_xxyz[k] * cd_x[k] + g_z_0_yyz_xxxyz[k];

            g_z_0_xyyz_xxzz[k] = -g_z_0_yyz_xxzz[k] * cd_x[k] + g_z_0_yyz_xxxzz[k];

            g_z_0_xyyz_xyyy[k] = -g_z_0_yyz_xyyy[k] * cd_x[k] + g_z_0_yyz_xxyyy[k];

            g_z_0_xyyz_xyyz[k] = -g_z_0_yyz_xyyz[k] * cd_x[k] + g_z_0_yyz_xxyyz[k];

            g_z_0_xyyz_xyzz[k] = -g_z_0_yyz_xyzz[k] * cd_x[k] + g_z_0_yyz_xxyzz[k];

            g_z_0_xyyz_xzzz[k] = -g_z_0_yyz_xzzz[k] * cd_x[k] + g_z_0_yyz_xxzzz[k];

            g_z_0_xyyz_yyyy[k] = -g_z_0_yyz_yyyy[k] * cd_x[k] + g_z_0_yyz_xyyyy[k];

            g_z_0_xyyz_yyyz[k] = -g_z_0_yyz_yyyz[k] * cd_x[k] + g_z_0_yyz_xyyyz[k];

            g_z_0_xyyz_yyzz[k] = -g_z_0_yyz_yyzz[k] * cd_x[k] + g_z_0_yyz_xyyzz[k];

            g_z_0_xyyz_yzzz[k] = -g_z_0_yyz_yzzz[k] * cd_x[k] + g_z_0_yyz_xyzzz[k];

            g_z_0_xyyz_zzzz[k] = -g_z_0_yyz_zzzz[k] * cd_x[k] + g_z_0_yyz_xzzzz[k];
        }

        /// Set up 120-135 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 120);

        auto g_z_0_xyzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 121);

        auto g_z_0_xyzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 122);

        auto g_z_0_xyzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 123);

        auto g_z_0_xyzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 124);

        auto g_z_0_xyzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 125);

        auto g_z_0_xyzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 126);

        auto g_z_0_xyzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 127);

        auto g_z_0_xyzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 128);

        auto g_z_0_xyzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 129);

        auto g_z_0_xyzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 130);

        auto g_z_0_xyzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 131);

        auto g_z_0_xyzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 132);

        auto g_z_0_xyzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 133);

        auto g_z_0_xyzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 134);

        #pragma omp simd aligned(cd_x, g_z_0_xyzz_xxxx, g_z_0_xyzz_xxxy, g_z_0_xyzz_xxxz, g_z_0_xyzz_xxyy, g_z_0_xyzz_xxyz, g_z_0_xyzz_xxzz, g_z_0_xyzz_xyyy, g_z_0_xyzz_xyyz, g_z_0_xyzz_xyzz, g_z_0_xyzz_xzzz, g_z_0_xyzz_yyyy, g_z_0_xyzz_yyyz, g_z_0_xyzz_yyzz, g_z_0_xyzz_yzzz, g_z_0_xyzz_zzzz, g_z_0_yzz_xxxx, g_z_0_yzz_xxxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxxz, g_z_0_yzz_xxxy, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxz, g_z_0_yzz_xxxzz, g_z_0_yzz_xxyy, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzz, g_z_0_yzz_xxzzz, g_z_0_yzz_xyyy, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzz, g_z_0_yzz_xzzzz, g_z_0_yzz_yyyy, g_z_0_yzz_yyyz, g_z_0_yzz_yyzz, g_z_0_yzz_yzzz, g_z_0_yzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyzz_xxxx[k] = -g_z_0_yzz_xxxx[k] * cd_x[k] + g_z_0_yzz_xxxxx[k];

            g_z_0_xyzz_xxxy[k] = -g_z_0_yzz_xxxy[k] * cd_x[k] + g_z_0_yzz_xxxxy[k];

            g_z_0_xyzz_xxxz[k] = -g_z_0_yzz_xxxz[k] * cd_x[k] + g_z_0_yzz_xxxxz[k];

            g_z_0_xyzz_xxyy[k] = -g_z_0_yzz_xxyy[k] * cd_x[k] + g_z_0_yzz_xxxyy[k];

            g_z_0_xyzz_xxyz[k] = -g_z_0_yzz_xxyz[k] * cd_x[k] + g_z_0_yzz_xxxyz[k];

            g_z_0_xyzz_xxzz[k] = -g_z_0_yzz_xxzz[k] * cd_x[k] + g_z_0_yzz_xxxzz[k];

            g_z_0_xyzz_xyyy[k] = -g_z_0_yzz_xyyy[k] * cd_x[k] + g_z_0_yzz_xxyyy[k];

            g_z_0_xyzz_xyyz[k] = -g_z_0_yzz_xyyz[k] * cd_x[k] + g_z_0_yzz_xxyyz[k];

            g_z_0_xyzz_xyzz[k] = -g_z_0_yzz_xyzz[k] * cd_x[k] + g_z_0_yzz_xxyzz[k];

            g_z_0_xyzz_xzzz[k] = -g_z_0_yzz_xzzz[k] * cd_x[k] + g_z_0_yzz_xxzzz[k];

            g_z_0_xyzz_yyyy[k] = -g_z_0_yzz_yyyy[k] * cd_x[k] + g_z_0_yzz_xyyyy[k];

            g_z_0_xyzz_yyyz[k] = -g_z_0_yzz_yyyz[k] * cd_x[k] + g_z_0_yzz_xyyyz[k];

            g_z_0_xyzz_yyzz[k] = -g_z_0_yzz_yyzz[k] * cd_x[k] + g_z_0_yzz_xyyzz[k];

            g_z_0_xyzz_yzzz[k] = -g_z_0_yzz_yzzz[k] * cd_x[k] + g_z_0_yzz_xyzzz[k];

            g_z_0_xyzz_zzzz[k] = -g_z_0_yzz_zzzz[k] * cd_x[k] + g_z_0_yzz_xzzzz[k];
        }

        /// Set up 135-150 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 135);

        auto g_z_0_xzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 136);

        auto g_z_0_xzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 137);

        auto g_z_0_xzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 138);

        auto g_z_0_xzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 139);

        auto g_z_0_xzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 140);

        auto g_z_0_xzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 141);

        auto g_z_0_xzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 142);

        auto g_z_0_xzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 143);

        auto g_z_0_xzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 144);

        auto g_z_0_xzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 145);

        auto g_z_0_xzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 146);

        auto g_z_0_xzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 147);

        auto g_z_0_xzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 148);

        auto g_z_0_xzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 149);

        #pragma omp simd aligned(cd_x, g_z_0_xzzz_xxxx, g_z_0_xzzz_xxxy, g_z_0_xzzz_xxxz, g_z_0_xzzz_xxyy, g_z_0_xzzz_xxyz, g_z_0_xzzz_xxzz, g_z_0_xzzz_xyyy, g_z_0_xzzz_xyyz, g_z_0_xzzz_xyzz, g_z_0_xzzz_xzzz, g_z_0_xzzz_yyyy, g_z_0_xzzz_yyyz, g_z_0_xzzz_yyzz, g_z_0_xzzz_yzzz, g_z_0_xzzz_zzzz, g_z_0_zzz_xxxx, g_z_0_zzz_xxxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyz, g_z_0_zzz_yyzz, g_z_0_zzz_yzzz, g_z_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzzz_xxxx[k] = -g_z_0_zzz_xxxx[k] * cd_x[k] + g_z_0_zzz_xxxxx[k];

            g_z_0_xzzz_xxxy[k] = -g_z_0_zzz_xxxy[k] * cd_x[k] + g_z_0_zzz_xxxxy[k];

            g_z_0_xzzz_xxxz[k] = -g_z_0_zzz_xxxz[k] * cd_x[k] + g_z_0_zzz_xxxxz[k];

            g_z_0_xzzz_xxyy[k] = -g_z_0_zzz_xxyy[k] * cd_x[k] + g_z_0_zzz_xxxyy[k];

            g_z_0_xzzz_xxyz[k] = -g_z_0_zzz_xxyz[k] * cd_x[k] + g_z_0_zzz_xxxyz[k];

            g_z_0_xzzz_xxzz[k] = -g_z_0_zzz_xxzz[k] * cd_x[k] + g_z_0_zzz_xxxzz[k];

            g_z_0_xzzz_xyyy[k] = -g_z_0_zzz_xyyy[k] * cd_x[k] + g_z_0_zzz_xxyyy[k];

            g_z_0_xzzz_xyyz[k] = -g_z_0_zzz_xyyz[k] * cd_x[k] + g_z_0_zzz_xxyyz[k];

            g_z_0_xzzz_xyzz[k] = -g_z_0_zzz_xyzz[k] * cd_x[k] + g_z_0_zzz_xxyzz[k];

            g_z_0_xzzz_xzzz[k] = -g_z_0_zzz_xzzz[k] * cd_x[k] + g_z_0_zzz_xxzzz[k];

            g_z_0_xzzz_yyyy[k] = -g_z_0_zzz_yyyy[k] * cd_x[k] + g_z_0_zzz_xyyyy[k];

            g_z_0_xzzz_yyyz[k] = -g_z_0_zzz_yyyz[k] * cd_x[k] + g_z_0_zzz_xyyyz[k];

            g_z_0_xzzz_yyzz[k] = -g_z_0_zzz_yyzz[k] * cd_x[k] + g_z_0_zzz_xyyzz[k];

            g_z_0_xzzz_yzzz[k] = -g_z_0_zzz_yzzz[k] * cd_x[k] + g_z_0_zzz_xyzzz[k];

            g_z_0_xzzz_zzzz[k] = -g_z_0_zzz_zzzz[k] * cd_x[k] + g_z_0_zzz_xzzzz[k];
        }

        /// Set up 150-165 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyy_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 150);

        auto g_z_0_yyyy_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 151);

        auto g_z_0_yyyy_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 152);

        auto g_z_0_yyyy_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 153);

        auto g_z_0_yyyy_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 154);

        auto g_z_0_yyyy_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 155);

        auto g_z_0_yyyy_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 156);

        auto g_z_0_yyyy_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 157);

        auto g_z_0_yyyy_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 158);

        auto g_z_0_yyyy_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 159);

        auto g_z_0_yyyy_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 160);

        auto g_z_0_yyyy_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 161);

        auto g_z_0_yyyy_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 162);

        auto g_z_0_yyyy_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 163);

        auto g_z_0_yyyy_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 164);

        #pragma omp simd aligned(cd_y, g_z_0_yyy_xxxx, g_z_0_yyy_xxxxy, g_z_0_yyy_xxxy, g_z_0_yyy_xxxyy, g_z_0_yyy_xxxyz, g_z_0_yyy_xxxz, g_z_0_yyy_xxyy, g_z_0_yyy_xxyyy, g_z_0_yyy_xxyyz, g_z_0_yyy_xxyz, g_z_0_yyy_xxyzz, g_z_0_yyy_xxzz, g_z_0_yyy_xyyy, g_z_0_yyy_xyyyy, g_z_0_yyy_xyyyz, g_z_0_yyy_xyyz, g_z_0_yyy_xyyzz, g_z_0_yyy_xyzz, g_z_0_yyy_xyzzz, g_z_0_yyy_xzzz, g_z_0_yyy_yyyy, g_z_0_yyy_yyyyy, g_z_0_yyy_yyyyz, g_z_0_yyy_yyyz, g_z_0_yyy_yyyzz, g_z_0_yyy_yyzz, g_z_0_yyy_yyzzz, g_z_0_yyy_yzzz, g_z_0_yyy_yzzzz, g_z_0_yyy_zzzz, g_z_0_yyyy_xxxx, g_z_0_yyyy_xxxy, g_z_0_yyyy_xxxz, g_z_0_yyyy_xxyy, g_z_0_yyyy_xxyz, g_z_0_yyyy_xxzz, g_z_0_yyyy_xyyy, g_z_0_yyyy_xyyz, g_z_0_yyyy_xyzz, g_z_0_yyyy_xzzz, g_z_0_yyyy_yyyy, g_z_0_yyyy_yyyz, g_z_0_yyyy_yyzz, g_z_0_yyyy_yzzz, g_z_0_yyyy_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyy_xxxx[k] = -g_z_0_yyy_xxxx[k] * cd_y[k] + g_z_0_yyy_xxxxy[k];

            g_z_0_yyyy_xxxy[k] = -g_z_0_yyy_xxxy[k] * cd_y[k] + g_z_0_yyy_xxxyy[k];

            g_z_0_yyyy_xxxz[k] = -g_z_0_yyy_xxxz[k] * cd_y[k] + g_z_0_yyy_xxxyz[k];

            g_z_0_yyyy_xxyy[k] = -g_z_0_yyy_xxyy[k] * cd_y[k] + g_z_0_yyy_xxyyy[k];

            g_z_0_yyyy_xxyz[k] = -g_z_0_yyy_xxyz[k] * cd_y[k] + g_z_0_yyy_xxyyz[k];

            g_z_0_yyyy_xxzz[k] = -g_z_0_yyy_xxzz[k] * cd_y[k] + g_z_0_yyy_xxyzz[k];

            g_z_0_yyyy_xyyy[k] = -g_z_0_yyy_xyyy[k] * cd_y[k] + g_z_0_yyy_xyyyy[k];

            g_z_0_yyyy_xyyz[k] = -g_z_0_yyy_xyyz[k] * cd_y[k] + g_z_0_yyy_xyyyz[k];

            g_z_0_yyyy_xyzz[k] = -g_z_0_yyy_xyzz[k] * cd_y[k] + g_z_0_yyy_xyyzz[k];

            g_z_0_yyyy_xzzz[k] = -g_z_0_yyy_xzzz[k] * cd_y[k] + g_z_0_yyy_xyzzz[k];

            g_z_0_yyyy_yyyy[k] = -g_z_0_yyy_yyyy[k] * cd_y[k] + g_z_0_yyy_yyyyy[k];

            g_z_0_yyyy_yyyz[k] = -g_z_0_yyy_yyyz[k] * cd_y[k] + g_z_0_yyy_yyyyz[k];

            g_z_0_yyyy_yyzz[k] = -g_z_0_yyy_yyzz[k] * cd_y[k] + g_z_0_yyy_yyyzz[k];

            g_z_0_yyyy_yzzz[k] = -g_z_0_yyy_yzzz[k] * cd_y[k] + g_z_0_yyy_yyzzz[k];

            g_z_0_yyyy_zzzz[k] = -g_z_0_yyy_zzzz[k] * cd_y[k] + g_z_0_yyy_yzzzz[k];
        }

        /// Set up 165-180 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyyz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 165);

        auto g_z_0_yyyz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 166);

        auto g_z_0_yyyz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 167);

        auto g_z_0_yyyz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 168);

        auto g_z_0_yyyz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 169);

        auto g_z_0_yyyz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 170);

        auto g_z_0_yyyz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 171);

        auto g_z_0_yyyz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 172);

        auto g_z_0_yyyz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 173);

        auto g_z_0_yyyz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 174);

        auto g_z_0_yyyz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 175);

        auto g_z_0_yyyz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 176);

        auto g_z_0_yyyz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 177);

        auto g_z_0_yyyz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 178);

        auto g_z_0_yyyz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 179);

        #pragma omp simd aligned(cd_y, g_z_0_yyyz_xxxx, g_z_0_yyyz_xxxy, g_z_0_yyyz_xxxz, g_z_0_yyyz_xxyy, g_z_0_yyyz_xxyz, g_z_0_yyyz_xxzz, g_z_0_yyyz_xyyy, g_z_0_yyyz_xyyz, g_z_0_yyyz_xyzz, g_z_0_yyyz_xzzz, g_z_0_yyyz_yyyy, g_z_0_yyyz_yyyz, g_z_0_yyyz_yyzz, g_z_0_yyyz_yzzz, g_z_0_yyyz_zzzz, g_z_0_yyz_xxxx, g_z_0_yyz_xxxxy, g_z_0_yyz_xxxy, g_z_0_yyz_xxxyy, g_z_0_yyz_xxxyz, g_z_0_yyz_xxxz, g_z_0_yyz_xxyy, g_z_0_yyz_xxyyy, g_z_0_yyz_xxyyz, g_z_0_yyz_xxyz, g_z_0_yyz_xxyzz, g_z_0_yyz_xxzz, g_z_0_yyz_xyyy, g_z_0_yyz_xyyyy, g_z_0_yyz_xyyyz, g_z_0_yyz_xyyz, g_z_0_yyz_xyyzz, g_z_0_yyz_xyzz, g_z_0_yyz_xyzzz, g_z_0_yyz_xzzz, g_z_0_yyz_yyyy, g_z_0_yyz_yyyyy, g_z_0_yyz_yyyyz, g_z_0_yyz_yyyz, g_z_0_yyz_yyyzz, g_z_0_yyz_yyzz, g_z_0_yyz_yyzzz, g_z_0_yyz_yzzz, g_z_0_yyz_yzzzz, g_z_0_yyz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyyz_xxxx[k] = -g_z_0_yyz_xxxx[k] * cd_y[k] + g_z_0_yyz_xxxxy[k];

            g_z_0_yyyz_xxxy[k] = -g_z_0_yyz_xxxy[k] * cd_y[k] + g_z_0_yyz_xxxyy[k];

            g_z_0_yyyz_xxxz[k] = -g_z_0_yyz_xxxz[k] * cd_y[k] + g_z_0_yyz_xxxyz[k];

            g_z_0_yyyz_xxyy[k] = -g_z_0_yyz_xxyy[k] * cd_y[k] + g_z_0_yyz_xxyyy[k];

            g_z_0_yyyz_xxyz[k] = -g_z_0_yyz_xxyz[k] * cd_y[k] + g_z_0_yyz_xxyyz[k];

            g_z_0_yyyz_xxzz[k] = -g_z_0_yyz_xxzz[k] * cd_y[k] + g_z_0_yyz_xxyzz[k];

            g_z_0_yyyz_xyyy[k] = -g_z_0_yyz_xyyy[k] * cd_y[k] + g_z_0_yyz_xyyyy[k];

            g_z_0_yyyz_xyyz[k] = -g_z_0_yyz_xyyz[k] * cd_y[k] + g_z_0_yyz_xyyyz[k];

            g_z_0_yyyz_xyzz[k] = -g_z_0_yyz_xyzz[k] * cd_y[k] + g_z_0_yyz_xyyzz[k];

            g_z_0_yyyz_xzzz[k] = -g_z_0_yyz_xzzz[k] * cd_y[k] + g_z_0_yyz_xyzzz[k];

            g_z_0_yyyz_yyyy[k] = -g_z_0_yyz_yyyy[k] * cd_y[k] + g_z_0_yyz_yyyyy[k];

            g_z_0_yyyz_yyyz[k] = -g_z_0_yyz_yyyz[k] * cd_y[k] + g_z_0_yyz_yyyyz[k];

            g_z_0_yyyz_yyzz[k] = -g_z_0_yyz_yyzz[k] * cd_y[k] + g_z_0_yyz_yyyzz[k];

            g_z_0_yyyz_yzzz[k] = -g_z_0_yyz_yzzz[k] * cd_y[k] + g_z_0_yyz_yyzzz[k];

            g_z_0_yyyz_zzzz[k] = -g_z_0_yyz_zzzz[k] * cd_y[k] + g_z_0_yyz_yzzzz[k];
        }

        /// Set up 180-195 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 180);

        auto g_z_0_yyzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 181);

        auto g_z_0_yyzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 182);

        auto g_z_0_yyzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 183);

        auto g_z_0_yyzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 184);

        auto g_z_0_yyzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 185);

        auto g_z_0_yyzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 186);

        auto g_z_0_yyzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 187);

        auto g_z_0_yyzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 188);

        auto g_z_0_yyzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 189);

        auto g_z_0_yyzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 190);

        auto g_z_0_yyzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 191);

        auto g_z_0_yyzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 192);

        auto g_z_0_yyzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 193);

        auto g_z_0_yyzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 194);

        #pragma omp simd aligned(cd_y, g_z_0_yyzz_xxxx, g_z_0_yyzz_xxxy, g_z_0_yyzz_xxxz, g_z_0_yyzz_xxyy, g_z_0_yyzz_xxyz, g_z_0_yyzz_xxzz, g_z_0_yyzz_xyyy, g_z_0_yyzz_xyyz, g_z_0_yyzz_xyzz, g_z_0_yyzz_xzzz, g_z_0_yyzz_yyyy, g_z_0_yyzz_yyyz, g_z_0_yyzz_yyzz, g_z_0_yyzz_yzzz, g_z_0_yyzz_zzzz, g_z_0_yzz_xxxx, g_z_0_yzz_xxxxy, g_z_0_yzz_xxxy, g_z_0_yzz_xxxyy, g_z_0_yzz_xxxyz, g_z_0_yzz_xxxz, g_z_0_yzz_xxyy, g_z_0_yzz_xxyyy, g_z_0_yzz_xxyyz, g_z_0_yzz_xxyz, g_z_0_yzz_xxyzz, g_z_0_yzz_xxzz, g_z_0_yzz_xyyy, g_z_0_yzz_xyyyy, g_z_0_yzz_xyyyz, g_z_0_yzz_xyyz, g_z_0_yzz_xyyzz, g_z_0_yzz_xyzz, g_z_0_yzz_xyzzz, g_z_0_yzz_xzzz, g_z_0_yzz_yyyy, g_z_0_yzz_yyyyy, g_z_0_yzz_yyyyz, g_z_0_yzz_yyyz, g_z_0_yzz_yyyzz, g_z_0_yzz_yyzz, g_z_0_yzz_yyzzz, g_z_0_yzz_yzzz, g_z_0_yzz_yzzzz, g_z_0_yzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyzz_xxxx[k] = -g_z_0_yzz_xxxx[k] * cd_y[k] + g_z_0_yzz_xxxxy[k];

            g_z_0_yyzz_xxxy[k] = -g_z_0_yzz_xxxy[k] * cd_y[k] + g_z_0_yzz_xxxyy[k];

            g_z_0_yyzz_xxxz[k] = -g_z_0_yzz_xxxz[k] * cd_y[k] + g_z_0_yzz_xxxyz[k];

            g_z_0_yyzz_xxyy[k] = -g_z_0_yzz_xxyy[k] * cd_y[k] + g_z_0_yzz_xxyyy[k];

            g_z_0_yyzz_xxyz[k] = -g_z_0_yzz_xxyz[k] * cd_y[k] + g_z_0_yzz_xxyyz[k];

            g_z_0_yyzz_xxzz[k] = -g_z_0_yzz_xxzz[k] * cd_y[k] + g_z_0_yzz_xxyzz[k];

            g_z_0_yyzz_xyyy[k] = -g_z_0_yzz_xyyy[k] * cd_y[k] + g_z_0_yzz_xyyyy[k];

            g_z_0_yyzz_xyyz[k] = -g_z_0_yzz_xyyz[k] * cd_y[k] + g_z_0_yzz_xyyyz[k];

            g_z_0_yyzz_xyzz[k] = -g_z_0_yzz_xyzz[k] * cd_y[k] + g_z_0_yzz_xyyzz[k];

            g_z_0_yyzz_xzzz[k] = -g_z_0_yzz_xzzz[k] * cd_y[k] + g_z_0_yzz_xyzzz[k];

            g_z_0_yyzz_yyyy[k] = -g_z_0_yzz_yyyy[k] * cd_y[k] + g_z_0_yzz_yyyyy[k];

            g_z_0_yyzz_yyyz[k] = -g_z_0_yzz_yyyz[k] * cd_y[k] + g_z_0_yzz_yyyyz[k];

            g_z_0_yyzz_yyzz[k] = -g_z_0_yzz_yyzz[k] * cd_y[k] + g_z_0_yzz_yyyzz[k];

            g_z_0_yyzz_yzzz[k] = -g_z_0_yzz_yzzz[k] * cd_y[k] + g_z_0_yzz_yyzzz[k];

            g_z_0_yyzz_zzzz[k] = -g_z_0_yzz_zzzz[k] * cd_y[k] + g_z_0_yzz_yzzzz[k];
        }

        /// Set up 195-210 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 195);

        auto g_z_0_yzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 196);

        auto g_z_0_yzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 197);

        auto g_z_0_yzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 198);

        auto g_z_0_yzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 199);

        auto g_z_0_yzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 200);

        auto g_z_0_yzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 201);

        auto g_z_0_yzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 202);

        auto g_z_0_yzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 203);

        auto g_z_0_yzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 204);

        auto g_z_0_yzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 205);

        auto g_z_0_yzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 206);

        auto g_z_0_yzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 207);

        auto g_z_0_yzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 208);

        auto g_z_0_yzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 209);

        #pragma omp simd aligned(cd_y, g_z_0_yzzz_xxxx, g_z_0_yzzz_xxxy, g_z_0_yzzz_xxxz, g_z_0_yzzz_xxyy, g_z_0_yzzz_xxyz, g_z_0_yzzz_xxzz, g_z_0_yzzz_xyyy, g_z_0_yzzz_xyyz, g_z_0_yzzz_xyzz, g_z_0_yzzz_xzzz, g_z_0_yzzz_yyyy, g_z_0_yzzz_yyyz, g_z_0_yzzz_yyzz, g_z_0_yzzz_yzzz, g_z_0_yzzz_zzzz, g_z_0_zzz_xxxx, g_z_0_zzz_xxxxy, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzzz_xxxx[k] = -g_z_0_zzz_xxxx[k] * cd_y[k] + g_z_0_zzz_xxxxy[k];

            g_z_0_yzzz_xxxy[k] = -g_z_0_zzz_xxxy[k] * cd_y[k] + g_z_0_zzz_xxxyy[k];

            g_z_0_yzzz_xxxz[k] = -g_z_0_zzz_xxxz[k] * cd_y[k] + g_z_0_zzz_xxxyz[k];

            g_z_0_yzzz_xxyy[k] = -g_z_0_zzz_xxyy[k] * cd_y[k] + g_z_0_zzz_xxyyy[k];

            g_z_0_yzzz_xxyz[k] = -g_z_0_zzz_xxyz[k] * cd_y[k] + g_z_0_zzz_xxyyz[k];

            g_z_0_yzzz_xxzz[k] = -g_z_0_zzz_xxzz[k] * cd_y[k] + g_z_0_zzz_xxyzz[k];

            g_z_0_yzzz_xyyy[k] = -g_z_0_zzz_xyyy[k] * cd_y[k] + g_z_0_zzz_xyyyy[k];

            g_z_0_yzzz_xyyz[k] = -g_z_0_zzz_xyyz[k] * cd_y[k] + g_z_0_zzz_xyyyz[k];

            g_z_0_yzzz_xyzz[k] = -g_z_0_zzz_xyzz[k] * cd_y[k] + g_z_0_zzz_xyyzz[k];

            g_z_0_yzzz_xzzz[k] = -g_z_0_zzz_xzzz[k] * cd_y[k] + g_z_0_zzz_xyzzz[k];

            g_z_0_yzzz_yyyy[k] = -g_z_0_zzz_yyyy[k] * cd_y[k] + g_z_0_zzz_yyyyy[k];

            g_z_0_yzzz_yyyz[k] = -g_z_0_zzz_yyyz[k] * cd_y[k] + g_z_0_zzz_yyyyz[k];

            g_z_0_yzzz_yyzz[k] = -g_z_0_zzz_yyzz[k] * cd_y[k] + g_z_0_zzz_yyyzz[k];

            g_z_0_yzzz_yzzz[k] = -g_z_0_zzz_yzzz[k] * cd_y[k] + g_z_0_zzz_yyzzz[k];

            g_z_0_yzzz_zzzz[k] = -g_z_0_zzz_zzzz[k] * cd_y[k] + g_z_0_zzz_yzzzz[k];
        }

        /// Set up 210-225 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzzz_xxxx = cbuffer.data(gg_geom_10_off + 450 * acomps  + 210);

        auto g_z_0_zzzz_xxxy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 211);

        auto g_z_0_zzzz_xxxz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 212);

        auto g_z_0_zzzz_xxyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 213);

        auto g_z_0_zzzz_xxyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 214);

        auto g_z_0_zzzz_xxzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 215);

        auto g_z_0_zzzz_xyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 216);

        auto g_z_0_zzzz_xyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 217);

        auto g_z_0_zzzz_xyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 218);

        auto g_z_0_zzzz_xzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 219);

        auto g_z_0_zzzz_yyyy = cbuffer.data(gg_geom_10_off + 450 * acomps  + 220);

        auto g_z_0_zzzz_yyyz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 221);

        auto g_z_0_zzzz_yyzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 222);

        auto g_z_0_zzzz_yzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 223);

        auto g_z_0_zzzz_zzzz = cbuffer.data(gg_geom_10_off + 450 * acomps  + 224);

        #pragma omp simd aligned(cd_z, g_z_0_zzz_xxxx, g_z_0_zzz_xxxxz, g_z_0_zzz_xxxy, g_z_0_zzz_xxxyz, g_z_0_zzz_xxxz, g_z_0_zzz_xxxzz, g_z_0_zzz_xxyy, g_z_0_zzz_xxyyz, g_z_0_zzz_xxyz, g_z_0_zzz_xxyzz, g_z_0_zzz_xxzz, g_z_0_zzz_xxzzz, g_z_0_zzz_xyyy, g_z_0_zzz_xyyyz, g_z_0_zzz_xyyz, g_z_0_zzz_xyyzz, g_z_0_zzz_xyzz, g_z_0_zzz_xyzzz, g_z_0_zzz_xzzz, g_z_0_zzz_xzzzz, g_z_0_zzz_yyyy, g_z_0_zzz_yyyyz, g_z_0_zzz_yyyz, g_z_0_zzz_yyyzz, g_z_0_zzz_yyzz, g_z_0_zzz_yyzzz, g_z_0_zzz_yzzz, g_z_0_zzz_yzzzz, g_z_0_zzz_zzzz, g_z_0_zzz_zzzzz, g_z_0_zzzz_xxxx, g_z_0_zzzz_xxxy, g_z_0_zzzz_xxxz, g_z_0_zzzz_xxyy, g_z_0_zzzz_xxyz, g_z_0_zzzz_xxzz, g_z_0_zzzz_xyyy, g_z_0_zzzz_xyyz, g_z_0_zzzz_xyzz, g_z_0_zzzz_xzzz, g_z_0_zzzz_yyyy, g_z_0_zzzz_yyyz, g_z_0_zzzz_yyzz, g_z_0_zzzz_yzzz, g_z_0_zzzz_zzzz, g_zzz_xxxx, g_zzz_xxxy, g_zzz_xxxz, g_zzz_xxyy, g_zzz_xxyz, g_zzz_xxzz, g_zzz_xyyy, g_zzz_xyyz, g_zzz_xyzz, g_zzz_xzzz, g_zzz_yyyy, g_zzz_yyyz, g_zzz_yyzz, g_zzz_yzzz, g_zzz_zzzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzzz_xxxx[k] = -g_zzz_xxxx[k] - g_z_0_zzz_xxxx[k] * cd_z[k] + g_z_0_zzz_xxxxz[k];

            g_z_0_zzzz_xxxy[k] = -g_zzz_xxxy[k] - g_z_0_zzz_xxxy[k] * cd_z[k] + g_z_0_zzz_xxxyz[k];

            g_z_0_zzzz_xxxz[k] = -g_zzz_xxxz[k] - g_z_0_zzz_xxxz[k] * cd_z[k] + g_z_0_zzz_xxxzz[k];

            g_z_0_zzzz_xxyy[k] = -g_zzz_xxyy[k] - g_z_0_zzz_xxyy[k] * cd_z[k] + g_z_0_zzz_xxyyz[k];

            g_z_0_zzzz_xxyz[k] = -g_zzz_xxyz[k] - g_z_0_zzz_xxyz[k] * cd_z[k] + g_z_0_zzz_xxyzz[k];

            g_z_0_zzzz_xxzz[k] = -g_zzz_xxzz[k] - g_z_0_zzz_xxzz[k] * cd_z[k] + g_z_0_zzz_xxzzz[k];

            g_z_0_zzzz_xyyy[k] = -g_zzz_xyyy[k] - g_z_0_zzz_xyyy[k] * cd_z[k] + g_z_0_zzz_xyyyz[k];

            g_z_0_zzzz_xyyz[k] = -g_zzz_xyyz[k] - g_z_0_zzz_xyyz[k] * cd_z[k] + g_z_0_zzz_xyyzz[k];

            g_z_0_zzzz_xyzz[k] = -g_zzz_xyzz[k] - g_z_0_zzz_xyzz[k] * cd_z[k] + g_z_0_zzz_xyzzz[k];

            g_z_0_zzzz_xzzz[k] = -g_zzz_xzzz[k] - g_z_0_zzz_xzzz[k] * cd_z[k] + g_z_0_zzz_xzzzz[k];

            g_z_0_zzzz_yyyy[k] = -g_zzz_yyyy[k] - g_z_0_zzz_yyyy[k] * cd_z[k] + g_z_0_zzz_yyyyz[k];

            g_z_0_zzzz_yyyz[k] = -g_zzz_yyyz[k] - g_z_0_zzz_yyyz[k] * cd_z[k] + g_z_0_zzz_yyyzz[k];

            g_z_0_zzzz_yyzz[k] = -g_zzz_yyzz[k] - g_z_0_zzz_yyzz[k] * cd_z[k] + g_z_0_zzz_yyzzz[k];

            g_z_0_zzzz_yzzz[k] = -g_zzz_yzzz[k] - g_z_0_zzz_yzzz[k] * cd_z[k] + g_z_0_zzz_yzzzz[k];

            g_z_0_zzzz_zzzz[k] = -g_zzz_zzzz[k] - g_z_0_zzz_zzzz[k] * cd_z[k] + g_z_0_zzz_zzzzz[k];
        }
    }
}

} // t3ceri namespace

