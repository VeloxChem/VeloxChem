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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXFD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xfd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfd,
                                        const size_t idx_xdd,
                                        const size_t idx_geom_10_xdd,
                                        const size_t idx_geom_10_xdf,
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
        /// Set up components of auxilary buffer : SDD

        const auto dd_off = idx_xdd + i * 36;

        auto g_xx_xx = cbuffer.data(dd_off + 0);

        auto g_xx_xy = cbuffer.data(dd_off + 1);

        auto g_xx_xz = cbuffer.data(dd_off + 2);

        auto g_xx_yy = cbuffer.data(dd_off + 3);

        auto g_xx_yz = cbuffer.data(dd_off + 4);

        auto g_xx_zz = cbuffer.data(dd_off + 5);

        auto g_yy_xx = cbuffer.data(dd_off + 18);

        auto g_yy_xy = cbuffer.data(dd_off + 19);

        auto g_yy_xz = cbuffer.data(dd_off + 20);

        auto g_yy_yy = cbuffer.data(dd_off + 21);

        auto g_yy_yz = cbuffer.data(dd_off + 22);

        auto g_yy_zz = cbuffer.data(dd_off + 23);

        auto g_zz_xx = cbuffer.data(dd_off + 30);

        auto g_zz_xy = cbuffer.data(dd_off + 31);

        auto g_zz_xz = cbuffer.data(dd_off + 32);

        auto g_zz_yy = cbuffer.data(dd_off + 33);

        auto g_zz_yz = cbuffer.data(dd_off + 34);

        auto g_zz_zz = cbuffer.data(dd_off + 35);

        /// Set up components of auxilary buffer : SDD

        const auto dd_geom_10_off = idx_geom_10_xdd + i * 36;

        auto g_x_0_xx_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xy_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xy_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xy_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 12);

        auto g_x_0_xz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 17);

        auto g_x_0_yy_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 18);

        auto g_x_0_yy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 19);

        auto g_x_0_yy_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 20);

        auto g_x_0_yy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 21);

        auto g_x_0_yy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 22);

        auto g_x_0_yy_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 23);

        auto g_x_0_yz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 24);

        auto g_x_0_yz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 25);

        auto g_x_0_yz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 26);

        auto g_x_0_yz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 27);

        auto g_x_0_yz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 28);

        auto g_x_0_yz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 29);

        auto g_x_0_zz_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 30);

        auto g_x_0_zz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 31);

        auto g_x_0_zz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 32);

        auto g_x_0_zz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 33);

        auto g_x_0_zz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 34);

        auto g_x_0_zz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 35);

        auto g_y_0_xx_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 0);

        auto g_y_0_xx_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 1);

        auto g_y_0_xx_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 2);

        auto g_y_0_xx_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 3);

        auto g_y_0_xx_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 4);

        auto g_y_0_xx_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 5);

        auto g_y_0_xy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 6);

        auto g_y_0_xy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 7);

        auto g_y_0_xy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 8);

        auto g_y_0_xy_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 9);

        auto g_y_0_xy_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 10);

        auto g_y_0_xy_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 11);

        auto g_y_0_xz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 12);

        auto g_y_0_xz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 13);

        auto g_y_0_xz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 14);

        auto g_y_0_xz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 15);

        auto g_y_0_xz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 16);

        auto g_y_0_xz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 17);

        auto g_y_0_yy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 18);

        auto g_y_0_yy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 19);

        auto g_y_0_yy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 20);

        auto g_y_0_yy_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 21);

        auto g_y_0_yy_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 22);

        auto g_y_0_yy_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 23);

        auto g_y_0_yz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 24);

        auto g_y_0_yz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 25);

        auto g_y_0_yz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 26);

        auto g_y_0_yz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 27);

        auto g_y_0_yz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 28);

        auto g_y_0_yz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 29);

        auto g_y_0_zz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 30);

        auto g_y_0_zz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 31);

        auto g_y_0_zz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 32);

        auto g_y_0_zz_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 33);

        auto g_y_0_zz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 34);

        auto g_y_0_zz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 35);

        auto g_z_0_xx_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 0);

        auto g_z_0_xx_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 1);

        auto g_z_0_xx_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 2);

        auto g_z_0_xx_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 3);

        auto g_z_0_xx_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 4);

        auto g_z_0_xx_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 5);

        auto g_z_0_xy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 6);

        auto g_z_0_xy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 7);

        auto g_z_0_xy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 8);

        auto g_z_0_xy_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 9);

        auto g_z_0_xy_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 10);

        auto g_z_0_xy_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 11);

        auto g_z_0_xz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 12);

        auto g_z_0_xz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 13);

        auto g_z_0_xz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 14);

        auto g_z_0_xz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 15);

        auto g_z_0_xz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 16);

        auto g_z_0_xz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 17);

        auto g_z_0_yy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 18);

        auto g_z_0_yy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 19);

        auto g_z_0_yy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 20);

        auto g_z_0_yy_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 21);

        auto g_z_0_yy_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 22);

        auto g_z_0_yy_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 23);

        auto g_z_0_yz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 24);

        auto g_z_0_yz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 25);

        auto g_z_0_yz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 26);

        auto g_z_0_yz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 27);

        auto g_z_0_yz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 28);

        auto g_z_0_yz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 29);

        auto g_z_0_zz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 30);

        auto g_z_0_zz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 31);

        auto g_z_0_zz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 32);

        auto g_z_0_zz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 33);

        auto g_z_0_zz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 34);

        auto g_z_0_zz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 35);

        /// Set up components of auxilary buffer : SDF

        const auto df_geom_10_off = idx_geom_10_xdf + i * 60;

        auto g_x_0_xx_xxx = cbuffer.data(df_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xx_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xx_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xx_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 8);

        auto g_x_0_xx_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 11);

        auto g_x_0_xy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 17);

        auto g_x_0_xy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 18);

        auto g_x_0_xz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 21);

        auto g_x_0_xz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 22);

        auto g_x_0_xz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 23);

        auto g_x_0_xz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 24);

        auto g_x_0_xz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 25);

        auto g_x_0_xz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 26);

        auto g_x_0_xz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 27);

        auto g_x_0_xz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 28);

        auto g_x_0_xz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 29);

        auto g_x_0_yy_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 31);

        auto g_x_0_yy_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 33);

        auto g_x_0_yy_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 34);

        auto g_x_0_yy_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 36);

        auto g_x_0_yy_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 37);

        auto g_x_0_yy_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 38);

        auto g_x_0_yz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 41);

        auto g_x_0_yz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 43);

        auto g_x_0_yz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 44);

        auto g_x_0_yz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 46);

        auto g_x_0_yz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 47);

        auto g_x_0_yz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 48);

        auto g_x_0_zz_xxy = cbuffer.data(df_geom_10_off + 0 * acomps + 51);

        auto g_x_0_zz_xxz = cbuffer.data(df_geom_10_off + 0 * acomps + 52);

        auto g_x_0_zz_xyy = cbuffer.data(df_geom_10_off + 0 * acomps + 53);

        auto g_x_0_zz_xyz = cbuffer.data(df_geom_10_off + 0 * acomps + 54);

        auto g_x_0_zz_xzz = cbuffer.data(df_geom_10_off + 0 * acomps + 55);

        auto g_x_0_zz_yyy = cbuffer.data(df_geom_10_off + 0 * acomps + 56);

        auto g_x_0_zz_yyz = cbuffer.data(df_geom_10_off + 0 * acomps + 57);

        auto g_x_0_zz_yzz = cbuffer.data(df_geom_10_off + 0 * acomps + 58);

        auto g_x_0_zz_zzz = cbuffer.data(df_geom_10_off + 0 * acomps + 59);

        auto g_y_0_xx_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 0);

        auto g_y_0_xx_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 1);

        auto g_y_0_xx_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 2);

        auto g_y_0_xx_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 3);

        auto g_y_0_xx_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 4);

        auto g_y_0_xx_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 5);

        auto g_y_0_xy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 10);

        auto g_y_0_xy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 11);

        auto g_y_0_xy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 12);

        auto g_y_0_xy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 13);

        auto g_y_0_xy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 14);

        auto g_y_0_xy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 15);

        auto g_y_0_xz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 20);

        auto g_y_0_xz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 21);

        auto g_y_0_xz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 22);

        auto g_y_0_xz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 23);

        auto g_y_0_xz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 24);

        auto g_y_0_xz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 25);

        auto g_y_0_yy_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 30);

        auto g_y_0_yy_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 31);

        auto g_y_0_yy_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 32);

        auto g_y_0_yy_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 33);

        auto g_y_0_yy_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 34);

        auto g_y_0_yy_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 35);

        auto g_y_0_yy_yyy = cbuffer.data(df_geom_10_off + 60 * acomps + 36);

        auto g_y_0_yy_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 37);

        auto g_y_0_yy_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 38);

        auto g_y_0_yy_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 39);

        auto g_y_0_yz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 40);

        auto g_y_0_yz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 41);

        auto g_y_0_yz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 42);

        auto g_y_0_yz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 43);

        auto g_y_0_yz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 44);

        auto g_y_0_yz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 45);

        auto g_y_0_yz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 47);

        auto g_y_0_yz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 48);

        auto g_y_0_yz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 49);

        auto g_y_0_zz_xxx = cbuffer.data(df_geom_10_off + 60 * acomps + 50);

        auto g_y_0_zz_xxy = cbuffer.data(df_geom_10_off + 60 * acomps + 51);

        auto g_y_0_zz_xxz = cbuffer.data(df_geom_10_off + 60 * acomps + 52);

        auto g_y_0_zz_xyy = cbuffer.data(df_geom_10_off + 60 * acomps + 53);

        auto g_y_0_zz_xyz = cbuffer.data(df_geom_10_off + 60 * acomps + 54);

        auto g_y_0_zz_xzz = cbuffer.data(df_geom_10_off + 60 * acomps + 55);

        auto g_y_0_zz_yyz = cbuffer.data(df_geom_10_off + 60 * acomps + 57);

        auto g_y_0_zz_yzz = cbuffer.data(df_geom_10_off + 60 * acomps + 58);

        auto g_y_0_zz_zzz = cbuffer.data(df_geom_10_off + 60 * acomps + 59);

        auto g_z_0_xx_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 0);

        auto g_z_0_xx_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 1);

        auto g_z_0_xx_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 2);

        auto g_z_0_xx_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 3);

        auto g_z_0_xx_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 4);

        auto g_z_0_xx_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 5);

        auto g_z_0_xy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 10);

        auto g_z_0_xy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 11);

        auto g_z_0_xy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 12);

        auto g_z_0_xy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 13);

        auto g_z_0_xy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 14);

        auto g_z_0_xy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 15);

        auto g_z_0_xz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 20);

        auto g_z_0_xz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 21);

        auto g_z_0_xz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 22);

        auto g_z_0_xz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 23);

        auto g_z_0_xz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 24);

        auto g_z_0_xz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 25);

        auto g_z_0_yy_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 30);

        auto g_z_0_yy_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 31);

        auto g_z_0_yy_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 32);

        auto g_z_0_yy_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 33);

        auto g_z_0_yy_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 34);

        auto g_z_0_yy_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 35);

        auto g_z_0_yy_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 36);

        auto g_z_0_yy_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 37);

        auto g_z_0_yy_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 38);

        auto g_z_0_yz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 40);

        auto g_z_0_yz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 41);

        auto g_z_0_yz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 42);

        auto g_z_0_yz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 43);

        auto g_z_0_yz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 44);

        auto g_z_0_yz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 45);

        auto g_z_0_yz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 46);

        auto g_z_0_yz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 47);

        auto g_z_0_yz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 48);

        auto g_z_0_zz_xxx = cbuffer.data(df_geom_10_off + 120 * acomps + 50);

        auto g_z_0_zz_xxy = cbuffer.data(df_geom_10_off + 120 * acomps + 51);

        auto g_z_0_zz_xxz = cbuffer.data(df_geom_10_off + 120 * acomps + 52);

        auto g_z_0_zz_xyy = cbuffer.data(df_geom_10_off + 120 * acomps + 53);

        auto g_z_0_zz_xyz = cbuffer.data(df_geom_10_off + 120 * acomps + 54);

        auto g_z_0_zz_xzz = cbuffer.data(df_geom_10_off + 120 * acomps + 55);

        auto g_z_0_zz_yyy = cbuffer.data(df_geom_10_off + 120 * acomps + 56);

        auto g_z_0_zz_yyz = cbuffer.data(df_geom_10_off + 120 * acomps + 57);

        auto g_z_0_zz_yzz = cbuffer.data(df_geom_10_off + 120 * acomps + 58);

        auto g_z_0_zz_zzz = cbuffer.data(df_geom_10_off + 120 * acomps + 59);

        /// set up bra offset for contr_buffer_xxfd

        const auto fd_geom_10_off = idx_geom_10_xfd + i * 60;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxx_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxx_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_xxx_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxx_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxx_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_x_0_xx_xx, g_x_0_xx_xxx, g_x_0_xx_xxy, g_x_0_xx_xxz, g_x_0_xx_xy, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_xzz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_0_xxx_xx, g_x_0_xxx_xy, g_x_0_xxx_xz, g_x_0_xxx_yy, g_x_0_xxx_yz, g_x_0_xxx_zz, g_xx_xx, g_xx_xy, g_xx_xz, g_xx_yy, g_xx_yz, g_xx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_xx[k] = -g_xx_xx[k] - g_x_0_xx_xx[k] * cd_x[k] + g_x_0_xx_xxx[k];

            g_x_0_xxx_xy[k] = -g_xx_xy[k] - g_x_0_xx_xy[k] * cd_x[k] + g_x_0_xx_xxy[k];

            g_x_0_xxx_xz[k] = -g_xx_xz[k] - g_x_0_xx_xz[k] * cd_x[k] + g_x_0_xx_xxz[k];

            g_x_0_xxx_yy[k] = -g_xx_yy[k] - g_x_0_xx_yy[k] * cd_x[k] + g_x_0_xx_xyy[k];

            g_x_0_xxx_yz[k] = -g_xx_yz[k] - g_x_0_xx_yz[k] * cd_x[k] + g_x_0_xx_xyz[k];

            g_x_0_xxx_zz[k] = -g_xx_zz[k] - g_x_0_xx_zz[k] * cd_x[k] + g_x_0_xx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_xxy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xxy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xxy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_xx_xx, g_x_0_xx_xxy, g_x_0_xx_xy, g_x_0_xx_xyy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yyy, g_x_0_xx_yyz, g_x_0_xx_yz, g_x_0_xx_yzz, g_x_0_xx_zz, g_x_0_xxy_xx, g_x_0_xxy_xy, g_x_0_xxy_xz, g_x_0_xxy_yy, g_x_0_xxy_yz, g_x_0_xxy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_xx[k] = -g_x_0_xx_xx[k] * cd_y[k] + g_x_0_xx_xxy[k];

            g_x_0_xxy_xy[k] = -g_x_0_xx_xy[k] * cd_y[k] + g_x_0_xx_xyy[k];

            g_x_0_xxy_xz[k] = -g_x_0_xx_xz[k] * cd_y[k] + g_x_0_xx_xyz[k];

            g_x_0_xxy_yy[k] = -g_x_0_xx_yy[k] * cd_y[k] + g_x_0_xx_yyy[k];

            g_x_0_xxy_yz[k] = -g_x_0_xx_yz[k] * cd_y[k] + g_x_0_xx_yyz[k];

            g_x_0_xxy_zz[k] = -g_x_0_xx_zz[k] * cd_y[k] + g_x_0_xx_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xxz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xxz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_xxz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xxz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xxz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_xx_xx, g_x_0_xx_xxz, g_x_0_xx_xy, g_x_0_xx_xyz, g_x_0_xx_xz, g_x_0_xx_xzz, g_x_0_xx_yy, g_x_0_xx_yyz, g_x_0_xx_yz, g_x_0_xx_yzz, g_x_0_xx_zz, g_x_0_xx_zzz, g_x_0_xxz_xx, g_x_0_xxz_xy, g_x_0_xxz_xz, g_x_0_xxz_yy, g_x_0_xxz_yz, g_x_0_xxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_xx[k] = -g_x_0_xx_xx[k] * cd_z[k] + g_x_0_xx_xxz[k];

            g_x_0_xxz_xy[k] = -g_x_0_xx_xy[k] * cd_z[k] + g_x_0_xx_xyz[k];

            g_x_0_xxz_xz[k] = -g_x_0_xx_xz[k] * cd_z[k] + g_x_0_xx_xzz[k];

            g_x_0_xxz_yy[k] = -g_x_0_xx_yy[k] * cd_z[k] + g_x_0_xx_yyz[k];

            g_x_0_xxz_yz[k] = -g_x_0_xx_yz[k] * cd_z[k] + g_x_0_xx_yzz[k];

            g_x_0_xxz_zz[k] = -g_x_0_xx_zz[k] * cd_z[k] + g_x_0_xx_zzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_xyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_xyy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 20);

        auto g_x_0_xyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_xyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_xyy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_x_0_xy_xx, g_x_0_xy_xxy, g_x_0_xy_xy, g_x_0_xy_xyy, g_x_0_xy_xyz, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yyy, g_x_0_xy_yyz, g_x_0_xy_yz, g_x_0_xy_yzz, g_x_0_xy_zz, g_x_0_xyy_xx, g_x_0_xyy_xy, g_x_0_xyy_xz, g_x_0_xyy_yy, g_x_0_xyy_yz, g_x_0_xyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_xx[k] = -g_x_0_xy_xx[k] * cd_y[k] + g_x_0_xy_xxy[k];

            g_x_0_xyy_xy[k] = -g_x_0_xy_xy[k] * cd_y[k] + g_x_0_xy_xyy[k];

            g_x_0_xyy_xz[k] = -g_x_0_xy_xz[k] * cd_y[k] + g_x_0_xy_xyz[k];

            g_x_0_xyy_yy[k] = -g_x_0_xy_yy[k] * cd_y[k] + g_x_0_xy_yyy[k];

            g_x_0_xyy_yz[k] = -g_x_0_xy_yz[k] * cd_y[k] + g_x_0_xy_yyz[k];

            g_x_0_xyy_zz[k] = -g_x_0_xy_zz[k] * cd_y[k] + g_x_0_xy_yzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_xyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_xyz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 26);

        auto g_x_0_xyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_xyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_xyz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_xx, g_x_0_xyz_xy, g_x_0_xyz_xz, g_x_0_xyz_yy, g_x_0_xyz_yz, g_x_0_xyz_zz, g_x_0_xz_xx, g_x_0_xz_xxy, g_x_0_xz_xy, g_x_0_xz_xyy, g_x_0_xz_xyz, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yyy, g_x_0_xz_yyz, g_x_0_xz_yz, g_x_0_xz_yzz, g_x_0_xz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_xx[k] = -g_x_0_xz_xx[k] * cd_y[k] + g_x_0_xz_xxy[k];

            g_x_0_xyz_xy[k] = -g_x_0_xz_xy[k] * cd_y[k] + g_x_0_xz_xyy[k];

            g_x_0_xyz_xz[k] = -g_x_0_xz_xz[k] * cd_y[k] + g_x_0_xz_xyz[k];

            g_x_0_xyz_yy[k] = -g_x_0_xz_yy[k] * cd_y[k] + g_x_0_xz_yyy[k];

            g_x_0_xyz_yz[k] = -g_x_0_xz_yz[k] * cd_y[k] + g_x_0_xz_yyz[k];

            g_x_0_xyz_zz[k] = -g_x_0_xz_zz[k] * cd_y[k] + g_x_0_xz_yzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 30);

        auto g_x_0_xzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 31);

        auto g_x_0_xzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 32);

        auto g_x_0_xzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 33);

        auto g_x_0_xzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 34);

        auto g_x_0_xzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 35);

        #pragma omp simd aligned(cd_z, g_x_0_xz_xx, g_x_0_xz_xxz, g_x_0_xz_xy, g_x_0_xz_xyz, g_x_0_xz_xz, g_x_0_xz_xzz, g_x_0_xz_yy, g_x_0_xz_yyz, g_x_0_xz_yz, g_x_0_xz_yzz, g_x_0_xz_zz, g_x_0_xz_zzz, g_x_0_xzz_xx, g_x_0_xzz_xy, g_x_0_xzz_xz, g_x_0_xzz_yy, g_x_0_xzz_yz, g_x_0_xzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_xx[k] = -g_x_0_xz_xx[k] * cd_z[k] + g_x_0_xz_xxz[k];

            g_x_0_xzz_xy[k] = -g_x_0_xz_xy[k] * cd_z[k] + g_x_0_xz_xyz[k];

            g_x_0_xzz_xz[k] = -g_x_0_xz_xz[k] * cd_z[k] + g_x_0_xz_xzz[k];

            g_x_0_xzz_yy[k] = -g_x_0_xz_yy[k] * cd_z[k] + g_x_0_xz_yyz[k];

            g_x_0_xzz_yz[k] = -g_x_0_xz_yz[k] * cd_z[k] + g_x_0_xz_yzz[k];

            g_x_0_xzz_zz[k] = -g_x_0_xz_zz[k] * cd_z[k] + g_x_0_xz_zzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 36);

        auto g_x_0_yyy_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 37);

        auto g_x_0_yyy_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 38);

        auto g_x_0_yyy_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 39);

        auto g_x_0_yyy_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 40);

        auto g_x_0_yyy_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_x_0_yy_xx, g_x_0_yy_xxy, g_x_0_yy_xy, g_x_0_yy_xyy, g_x_0_yy_xyz, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yyy, g_x_0_yy_yyz, g_x_0_yy_yz, g_x_0_yy_yzz, g_x_0_yy_zz, g_x_0_yyy_xx, g_x_0_yyy_xy, g_x_0_yyy_xz, g_x_0_yyy_yy, g_x_0_yyy_yz, g_x_0_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_xx[k] = -g_x_0_yy_xx[k] * cd_y[k] + g_x_0_yy_xxy[k];

            g_x_0_yyy_xy[k] = -g_x_0_yy_xy[k] * cd_y[k] + g_x_0_yy_xyy[k];

            g_x_0_yyy_xz[k] = -g_x_0_yy_xz[k] * cd_y[k] + g_x_0_yy_xyz[k];

            g_x_0_yyy_yy[k] = -g_x_0_yy_yy[k] * cd_y[k] + g_x_0_yy_yyy[k];

            g_x_0_yyy_yz[k] = -g_x_0_yy_yz[k] * cd_y[k] + g_x_0_yy_yyz[k];

            g_x_0_yyy_zz[k] = -g_x_0_yy_zz[k] * cd_y[k] + g_x_0_yy_yzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 42);

        auto g_x_0_yyz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 43);

        auto g_x_0_yyz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 44);

        auto g_x_0_yyz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 45);

        auto g_x_0_yyz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 46);

        auto g_x_0_yyz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 47);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_xx, g_x_0_yyz_xy, g_x_0_yyz_xz, g_x_0_yyz_yy, g_x_0_yyz_yz, g_x_0_yyz_zz, g_x_0_yz_xx, g_x_0_yz_xxy, g_x_0_yz_xy, g_x_0_yz_xyy, g_x_0_yz_xyz, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yyy, g_x_0_yz_yyz, g_x_0_yz_yz, g_x_0_yz_yzz, g_x_0_yz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_xx[k] = -g_x_0_yz_xx[k] * cd_y[k] + g_x_0_yz_xxy[k];

            g_x_0_yyz_xy[k] = -g_x_0_yz_xy[k] * cd_y[k] + g_x_0_yz_xyy[k];

            g_x_0_yyz_xz[k] = -g_x_0_yz_xz[k] * cd_y[k] + g_x_0_yz_xyz[k];

            g_x_0_yyz_yy[k] = -g_x_0_yz_yy[k] * cd_y[k] + g_x_0_yz_yyy[k];

            g_x_0_yyz_yz[k] = -g_x_0_yz_yz[k] * cd_y[k] + g_x_0_yz_yyz[k];

            g_x_0_yyz_zz[k] = -g_x_0_yz_zz[k] * cd_y[k] + g_x_0_yz_yzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 48);

        auto g_x_0_yzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 49);

        auto g_x_0_yzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 50);

        auto g_x_0_yzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 51);

        auto g_x_0_yzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 52);

        auto g_x_0_yzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 53);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_xx, g_x_0_yzz_xy, g_x_0_yzz_xz, g_x_0_yzz_yy, g_x_0_yzz_yz, g_x_0_yzz_zz, g_x_0_zz_xx, g_x_0_zz_xxy, g_x_0_zz_xy, g_x_0_zz_xyy, g_x_0_zz_xyz, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yyy, g_x_0_zz_yyz, g_x_0_zz_yz, g_x_0_zz_yzz, g_x_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_xx[k] = -g_x_0_zz_xx[k] * cd_y[k] + g_x_0_zz_xxy[k];

            g_x_0_yzz_xy[k] = -g_x_0_zz_xy[k] * cd_y[k] + g_x_0_zz_xyy[k];

            g_x_0_yzz_xz[k] = -g_x_0_zz_xz[k] * cd_y[k] + g_x_0_zz_xyz[k];

            g_x_0_yzz_yy[k] = -g_x_0_zz_yy[k] * cd_y[k] + g_x_0_zz_yyy[k];

            g_x_0_yzz_yz[k] = -g_x_0_zz_yz[k] * cd_y[k] + g_x_0_zz_yyz[k];

            g_x_0_yzz_zz[k] = -g_x_0_zz_zz[k] * cd_y[k] + g_x_0_zz_yzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_xx = cbuffer.data(fd_geom_10_off + 0 * acomps  + 54);

        auto g_x_0_zzz_xy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 55);

        auto g_x_0_zzz_xz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 56);

        auto g_x_0_zzz_yy = cbuffer.data(fd_geom_10_off + 0 * acomps  + 57);

        auto g_x_0_zzz_yz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 58);

        auto g_x_0_zzz_zz = cbuffer.data(fd_geom_10_off + 0 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_x_0_zz_xx, g_x_0_zz_xxz, g_x_0_zz_xy, g_x_0_zz_xyz, g_x_0_zz_xz, g_x_0_zz_xzz, g_x_0_zz_yy, g_x_0_zz_yyz, g_x_0_zz_yz, g_x_0_zz_yzz, g_x_0_zz_zz, g_x_0_zz_zzz, g_x_0_zzz_xx, g_x_0_zzz_xy, g_x_0_zzz_xz, g_x_0_zzz_yy, g_x_0_zzz_yz, g_x_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_xx[k] = -g_x_0_zz_xx[k] * cd_z[k] + g_x_0_zz_xxz[k];

            g_x_0_zzz_xy[k] = -g_x_0_zz_xy[k] * cd_z[k] + g_x_0_zz_xyz[k];

            g_x_0_zzz_xz[k] = -g_x_0_zz_xz[k] * cd_z[k] + g_x_0_zz_xzz[k];

            g_x_0_zzz_yy[k] = -g_x_0_zz_yy[k] * cd_z[k] + g_x_0_zz_yyz[k];

            g_x_0_zzz_yz[k] = -g_x_0_zz_yz[k] * cd_z[k] + g_x_0_zz_yzz[k];

            g_x_0_zzz_zz[k] = -g_x_0_zz_zz[k] * cd_z[k] + g_x_0_zz_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 0);

        auto g_y_0_xxx_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 1);

        auto g_y_0_xxx_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 2);

        auto g_y_0_xxx_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 3);

        auto g_y_0_xxx_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 4);

        auto g_y_0_xxx_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xx_xx, g_y_0_xx_xxx, g_y_0_xx_xxy, g_y_0_xx_xxz, g_y_0_xx_xy, g_y_0_xx_xyy, g_y_0_xx_xyz, g_y_0_xx_xz, g_y_0_xx_xzz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz, g_y_0_xxx_xx, g_y_0_xxx_xy, g_y_0_xxx_xz, g_y_0_xxx_yy, g_y_0_xxx_yz, g_y_0_xxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_xx[k] = -g_y_0_xx_xx[k] * cd_x[k] + g_y_0_xx_xxx[k];

            g_y_0_xxx_xy[k] = -g_y_0_xx_xy[k] * cd_x[k] + g_y_0_xx_xxy[k];

            g_y_0_xxx_xz[k] = -g_y_0_xx_xz[k] * cd_x[k] + g_y_0_xx_xxz[k];

            g_y_0_xxx_yy[k] = -g_y_0_xx_yy[k] * cd_x[k] + g_y_0_xx_xyy[k];

            g_y_0_xxx_yz[k] = -g_y_0_xx_yz[k] * cd_x[k] + g_y_0_xx_xyz[k];

            g_y_0_xxx_zz[k] = -g_y_0_xx_zz[k] * cd_x[k] + g_y_0_xx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 6);

        auto g_y_0_xxy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 7);

        auto g_y_0_xxy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 8);

        auto g_y_0_xxy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 9);

        auto g_y_0_xxy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 10);

        auto g_y_0_xxy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_xx, g_y_0_xxy_xy, g_y_0_xxy_xz, g_y_0_xxy_yy, g_y_0_xxy_yz, g_y_0_xxy_zz, g_y_0_xy_xx, g_y_0_xy_xxx, g_y_0_xy_xxy, g_y_0_xy_xxz, g_y_0_xy_xy, g_y_0_xy_xyy, g_y_0_xy_xyz, g_y_0_xy_xz, g_y_0_xy_xzz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_xx[k] = -g_y_0_xy_xx[k] * cd_x[k] + g_y_0_xy_xxx[k];

            g_y_0_xxy_xy[k] = -g_y_0_xy_xy[k] * cd_x[k] + g_y_0_xy_xxy[k];

            g_y_0_xxy_xz[k] = -g_y_0_xy_xz[k] * cd_x[k] + g_y_0_xy_xxz[k];

            g_y_0_xxy_yy[k] = -g_y_0_xy_yy[k] * cd_x[k] + g_y_0_xy_xyy[k];

            g_y_0_xxy_yz[k] = -g_y_0_xy_yz[k] * cd_x[k] + g_y_0_xy_xyz[k];

            g_y_0_xxy_zz[k] = -g_y_0_xy_zz[k] * cd_x[k] + g_y_0_xy_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 12);

        auto g_y_0_xxz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 13);

        auto g_y_0_xxz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 14);

        auto g_y_0_xxz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 15);

        auto g_y_0_xxz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 16);

        auto g_y_0_xxz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_xx, g_y_0_xxz_xy, g_y_0_xxz_xz, g_y_0_xxz_yy, g_y_0_xxz_yz, g_y_0_xxz_zz, g_y_0_xz_xx, g_y_0_xz_xxx, g_y_0_xz_xxy, g_y_0_xz_xxz, g_y_0_xz_xy, g_y_0_xz_xyy, g_y_0_xz_xyz, g_y_0_xz_xz, g_y_0_xz_xzz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_xx[k] = -g_y_0_xz_xx[k] * cd_x[k] + g_y_0_xz_xxx[k];

            g_y_0_xxz_xy[k] = -g_y_0_xz_xy[k] * cd_x[k] + g_y_0_xz_xxy[k];

            g_y_0_xxz_xz[k] = -g_y_0_xz_xz[k] * cd_x[k] + g_y_0_xz_xxz[k];

            g_y_0_xxz_yy[k] = -g_y_0_xz_yy[k] * cd_x[k] + g_y_0_xz_xyy[k];

            g_y_0_xxz_yz[k] = -g_y_0_xz_yz[k] * cd_x[k] + g_y_0_xz_xyz[k];

            g_y_0_xxz_zz[k] = -g_y_0_xz_zz[k] * cd_x[k] + g_y_0_xz_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 18);

        auto g_y_0_xyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 19);

        auto g_y_0_xyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 20);

        auto g_y_0_xyy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 21);

        auto g_y_0_xyy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 22);

        auto g_y_0_xyy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 23);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_xx, g_y_0_xyy_xy, g_y_0_xyy_xz, g_y_0_xyy_yy, g_y_0_xyy_yz, g_y_0_xyy_zz, g_y_0_yy_xx, g_y_0_yy_xxx, g_y_0_yy_xxy, g_y_0_yy_xxz, g_y_0_yy_xy, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_xzz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_xx[k] = -g_y_0_yy_xx[k] * cd_x[k] + g_y_0_yy_xxx[k];

            g_y_0_xyy_xy[k] = -g_y_0_yy_xy[k] * cd_x[k] + g_y_0_yy_xxy[k];

            g_y_0_xyy_xz[k] = -g_y_0_yy_xz[k] * cd_x[k] + g_y_0_yy_xxz[k];

            g_y_0_xyy_yy[k] = -g_y_0_yy_yy[k] * cd_x[k] + g_y_0_yy_xyy[k];

            g_y_0_xyy_yz[k] = -g_y_0_yy_yz[k] * cd_x[k] + g_y_0_yy_xyz[k];

            g_y_0_xyy_zz[k] = -g_y_0_yy_zz[k] * cd_x[k] + g_y_0_yy_xzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 24);

        auto g_y_0_xyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 25);

        auto g_y_0_xyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 26);

        auto g_y_0_xyz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 27);

        auto g_y_0_xyz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 28);

        auto g_y_0_xyz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_xx, g_y_0_xyz_xy, g_y_0_xyz_xz, g_y_0_xyz_yy, g_y_0_xyz_yz, g_y_0_xyz_zz, g_y_0_yz_xx, g_y_0_yz_xxx, g_y_0_yz_xxy, g_y_0_yz_xxz, g_y_0_yz_xy, g_y_0_yz_xyy, g_y_0_yz_xyz, g_y_0_yz_xz, g_y_0_yz_xzz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_xx[k] = -g_y_0_yz_xx[k] * cd_x[k] + g_y_0_yz_xxx[k];

            g_y_0_xyz_xy[k] = -g_y_0_yz_xy[k] * cd_x[k] + g_y_0_yz_xxy[k];

            g_y_0_xyz_xz[k] = -g_y_0_yz_xz[k] * cd_x[k] + g_y_0_yz_xxz[k];

            g_y_0_xyz_yy[k] = -g_y_0_yz_yy[k] * cd_x[k] + g_y_0_yz_xyy[k];

            g_y_0_xyz_yz[k] = -g_y_0_yz_yz[k] * cd_x[k] + g_y_0_yz_xyz[k];

            g_y_0_xyz_zz[k] = -g_y_0_yz_zz[k] * cd_x[k] + g_y_0_yz_xzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 30);

        auto g_y_0_xzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 31);

        auto g_y_0_xzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 32);

        auto g_y_0_xzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 33);

        auto g_y_0_xzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 34);

        auto g_y_0_xzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_xx, g_y_0_xzz_xy, g_y_0_xzz_xz, g_y_0_xzz_yy, g_y_0_xzz_yz, g_y_0_xzz_zz, g_y_0_zz_xx, g_y_0_zz_xxx, g_y_0_zz_xxy, g_y_0_zz_xxz, g_y_0_zz_xy, g_y_0_zz_xyy, g_y_0_zz_xyz, g_y_0_zz_xz, g_y_0_zz_xzz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_xx[k] = -g_y_0_zz_xx[k] * cd_x[k] + g_y_0_zz_xxx[k];

            g_y_0_xzz_xy[k] = -g_y_0_zz_xy[k] * cd_x[k] + g_y_0_zz_xxy[k];

            g_y_0_xzz_xz[k] = -g_y_0_zz_xz[k] * cd_x[k] + g_y_0_zz_xxz[k];

            g_y_0_xzz_yy[k] = -g_y_0_zz_yy[k] * cd_x[k] + g_y_0_zz_xyy[k];

            g_y_0_xzz_yz[k] = -g_y_0_zz_yz[k] * cd_x[k] + g_y_0_zz_xyz[k];

            g_y_0_xzz_zz[k] = -g_y_0_zz_zz[k] * cd_x[k] + g_y_0_zz_xzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 36);

        auto g_y_0_yyy_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 37);

        auto g_y_0_yyy_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 38);

        auto g_y_0_yyy_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 39);

        auto g_y_0_yyy_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 40);

        auto g_y_0_yyy_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_y_0_yy_xx, g_y_0_yy_xxy, g_y_0_yy_xy, g_y_0_yy_xyy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yyy, g_y_0_yy_yyz, g_y_0_yy_yz, g_y_0_yy_yzz, g_y_0_yy_zz, g_y_0_yyy_xx, g_y_0_yyy_xy, g_y_0_yyy_xz, g_y_0_yyy_yy, g_y_0_yyy_yz, g_y_0_yyy_zz, g_yy_xx, g_yy_xy, g_yy_xz, g_yy_yy, g_yy_yz, g_yy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_xx[k] = -g_yy_xx[k] - g_y_0_yy_xx[k] * cd_y[k] + g_y_0_yy_xxy[k];

            g_y_0_yyy_xy[k] = -g_yy_xy[k] - g_y_0_yy_xy[k] * cd_y[k] + g_y_0_yy_xyy[k];

            g_y_0_yyy_xz[k] = -g_yy_xz[k] - g_y_0_yy_xz[k] * cd_y[k] + g_y_0_yy_xyz[k];

            g_y_0_yyy_yy[k] = -g_yy_yy[k] - g_y_0_yy_yy[k] * cd_y[k] + g_y_0_yy_yyy[k];

            g_y_0_yyy_yz[k] = -g_yy_yz[k] - g_y_0_yy_yz[k] * cd_y[k] + g_y_0_yy_yyz[k];

            g_y_0_yyy_zz[k] = -g_yy_zz[k] - g_y_0_yy_zz[k] * cd_y[k] + g_y_0_yy_yzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 42);

        auto g_y_0_yyz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 43);

        auto g_y_0_yyz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 44);

        auto g_y_0_yyz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 45);

        auto g_y_0_yyz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 46);

        auto g_y_0_yyz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 47);

        #pragma omp simd aligned(cd_z, g_y_0_yy_xx, g_y_0_yy_xxz, g_y_0_yy_xy, g_y_0_yy_xyz, g_y_0_yy_xz, g_y_0_yy_xzz, g_y_0_yy_yy, g_y_0_yy_yyz, g_y_0_yy_yz, g_y_0_yy_yzz, g_y_0_yy_zz, g_y_0_yy_zzz, g_y_0_yyz_xx, g_y_0_yyz_xy, g_y_0_yyz_xz, g_y_0_yyz_yy, g_y_0_yyz_yz, g_y_0_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_xx[k] = -g_y_0_yy_xx[k] * cd_z[k] + g_y_0_yy_xxz[k];

            g_y_0_yyz_xy[k] = -g_y_0_yy_xy[k] * cd_z[k] + g_y_0_yy_xyz[k];

            g_y_0_yyz_xz[k] = -g_y_0_yy_xz[k] * cd_z[k] + g_y_0_yy_xzz[k];

            g_y_0_yyz_yy[k] = -g_y_0_yy_yy[k] * cd_z[k] + g_y_0_yy_yyz[k];

            g_y_0_yyz_yz[k] = -g_y_0_yy_yz[k] * cd_z[k] + g_y_0_yy_yzz[k];

            g_y_0_yyz_zz[k] = -g_y_0_yy_zz[k] * cd_z[k] + g_y_0_yy_zzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 48);

        auto g_y_0_yzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 49);

        auto g_y_0_yzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 50);

        auto g_y_0_yzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 51);

        auto g_y_0_yzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 52);

        auto g_y_0_yzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 53);

        #pragma omp simd aligned(cd_z, g_y_0_yz_xx, g_y_0_yz_xxz, g_y_0_yz_xy, g_y_0_yz_xyz, g_y_0_yz_xz, g_y_0_yz_xzz, g_y_0_yz_yy, g_y_0_yz_yyz, g_y_0_yz_yz, g_y_0_yz_yzz, g_y_0_yz_zz, g_y_0_yz_zzz, g_y_0_yzz_xx, g_y_0_yzz_xy, g_y_0_yzz_xz, g_y_0_yzz_yy, g_y_0_yzz_yz, g_y_0_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_xx[k] = -g_y_0_yz_xx[k] * cd_z[k] + g_y_0_yz_xxz[k];

            g_y_0_yzz_xy[k] = -g_y_0_yz_xy[k] * cd_z[k] + g_y_0_yz_xyz[k];

            g_y_0_yzz_xz[k] = -g_y_0_yz_xz[k] * cd_z[k] + g_y_0_yz_xzz[k];

            g_y_0_yzz_yy[k] = -g_y_0_yz_yy[k] * cd_z[k] + g_y_0_yz_yyz[k];

            g_y_0_yzz_yz[k] = -g_y_0_yz_yz[k] * cd_z[k] + g_y_0_yz_yzz[k];

            g_y_0_yzz_zz[k] = -g_y_0_yz_zz[k] * cd_z[k] + g_y_0_yz_zzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_xx = cbuffer.data(fd_geom_10_off + 60 * acomps  + 54);

        auto g_y_0_zzz_xy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 55);

        auto g_y_0_zzz_xz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 56);

        auto g_y_0_zzz_yy = cbuffer.data(fd_geom_10_off + 60 * acomps  + 57);

        auto g_y_0_zzz_yz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 58);

        auto g_y_0_zzz_zz = cbuffer.data(fd_geom_10_off + 60 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_y_0_zz_xx, g_y_0_zz_xxz, g_y_0_zz_xy, g_y_0_zz_xyz, g_y_0_zz_xz, g_y_0_zz_xzz, g_y_0_zz_yy, g_y_0_zz_yyz, g_y_0_zz_yz, g_y_0_zz_yzz, g_y_0_zz_zz, g_y_0_zz_zzz, g_y_0_zzz_xx, g_y_0_zzz_xy, g_y_0_zzz_xz, g_y_0_zzz_yy, g_y_0_zzz_yz, g_y_0_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_xx[k] = -g_y_0_zz_xx[k] * cd_z[k] + g_y_0_zz_xxz[k];

            g_y_0_zzz_xy[k] = -g_y_0_zz_xy[k] * cd_z[k] + g_y_0_zz_xyz[k];

            g_y_0_zzz_xz[k] = -g_y_0_zz_xz[k] * cd_z[k] + g_y_0_zz_xzz[k];

            g_y_0_zzz_yy[k] = -g_y_0_zz_yy[k] * cd_z[k] + g_y_0_zz_yyz[k];

            g_y_0_zzz_yz[k] = -g_y_0_zz_yz[k] * cd_z[k] + g_y_0_zz_yzz[k];

            g_y_0_zzz_zz[k] = -g_y_0_zz_zz[k] * cd_z[k] + g_y_0_zz_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 0);

        auto g_z_0_xxx_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 1);

        auto g_z_0_xxx_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 2);

        auto g_z_0_xxx_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 3);

        auto g_z_0_xxx_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 4);

        auto g_z_0_xxx_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xx_xx, g_z_0_xx_xxx, g_z_0_xx_xxy, g_z_0_xx_xxz, g_z_0_xx_xy, g_z_0_xx_xyy, g_z_0_xx_xyz, g_z_0_xx_xz, g_z_0_xx_xzz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz, g_z_0_xxx_xx, g_z_0_xxx_xy, g_z_0_xxx_xz, g_z_0_xxx_yy, g_z_0_xxx_yz, g_z_0_xxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_xx[k] = -g_z_0_xx_xx[k] * cd_x[k] + g_z_0_xx_xxx[k];

            g_z_0_xxx_xy[k] = -g_z_0_xx_xy[k] * cd_x[k] + g_z_0_xx_xxy[k];

            g_z_0_xxx_xz[k] = -g_z_0_xx_xz[k] * cd_x[k] + g_z_0_xx_xxz[k];

            g_z_0_xxx_yy[k] = -g_z_0_xx_yy[k] * cd_x[k] + g_z_0_xx_xyy[k];

            g_z_0_xxx_yz[k] = -g_z_0_xx_yz[k] * cd_x[k] + g_z_0_xx_xyz[k];

            g_z_0_xxx_zz[k] = -g_z_0_xx_zz[k] * cd_x[k] + g_z_0_xx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 6);

        auto g_z_0_xxy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 7);

        auto g_z_0_xxy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 8);

        auto g_z_0_xxy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 9);

        auto g_z_0_xxy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 10);

        auto g_z_0_xxy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_xx, g_z_0_xxy_xy, g_z_0_xxy_xz, g_z_0_xxy_yy, g_z_0_xxy_yz, g_z_0_xxy_zz, g_z_0_xy_xx, g_z_0_xy_xxx, g_z_0_xy_xxy, g_z_0_xy_xxz, g_z_0_xy_xy, g_z_0_xy_xyy, g_z_0_xy_xyz, g_z_0_xy_xz, g_z_0_xy_xzz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_xx[k] = -g_z_0_xy_xx[k] * cd_x[k] + g_z_0_xy_xxx[k];

            g_z_0_xxy_xy[k] = -g_z_0_xy_xy[k] * cd_x[k] + g_z_0_xy_xxy[k];

            g_z_0_xxy_xz[k] = -g_z_0_xy_xz[k] * cd_x[k] + g_z_0_xy_xxz[k];

            g_z_0_xxy_yy[k] = -g_z_0_xy_yy[k] * cd_x[k] + g_z_0_xy_xyy[k];

            g_z_0_xxy_yz[k] = -g_z_0_xy_yz[k] * cd_x[k] + g_z_0_xy_xyz[k];

            g_z_0_xxy_zz[k] = -g_z_0_xy_zz[k] * cd_x[k] + g_z_0_xy_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 12);

        auto g_z_0_xxz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 13);

        auto g_z_0_xxz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 14);

        auto g_z_0_xxz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 15);

        auto g_z_0_xxz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 16);

        auto g_z_0_xxz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_xx, g_z_0_xxz_xy, g_z_0_xxz_xz, g_z_0_xxz_yy, g_z_0_xxz_yz, g_z_0_xxz_zz, g_z_0_xz_xx, g_z_0_xz_xxx, g_z_0_xz_xxy, g_z_0_xz_xxz, g_z_0_xz_xy, g_z_0_xz_xyy, g_z_0_xz_xyz, g_z_0_xz_xz, g_z_0_xz_xzz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_xx[k] = -g_z_0_xz_xx[k] * cd_x[k] + g_z_0_xz_xxx[k];

            g_z_0_xxz_xy[k] = -g_z_0_xz_xy[k] * cd_x[k] + g_z_0_xz_xxy[k];

            g_z_0_xxz_xz[k] = -g_z_0_xz_xz[k] * cd_x[k] + g_z_0_xz_xxz[k];

            g_z_0_xxz_yy[k] = -g_z_0_xz_yy[k] * cd_x[k] + g_z_0_xz_xyy[k];

            g_z_0_xxz_yz[k] = -g_z_0_xz_yz[k] * cd_x[k] + g_z_0_xz_xyz[k];

            g_z_0_xxz_zz[k] = -g_z_0_xz_zz[k] * cd_x[k] + g_z_0_xz_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 18);

        auto g_z_0_xyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 19);

        auto g_z_0_xyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 20);

        auto g_z_0_xyy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 21);

        auto g_z_0_xyy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 22);

        auto g_z_0_xyy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 23);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_xx, g_z_0_xyy_xy, g_z_0_xyy_xz, g_z_0_xyy_yy, g_z_0_xyy_yz, g_z_0_xyy_zz, g_z_0_yy_xx, g_z_0_yy_xxx, g_z_0_yy_xxy, g_z_0_yy_xxz, g_z_0_yy_xy, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xz, g_z_0_yy_xzz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_xx[k] = -g_z_0_yy_xx[k] * cd_x[k] + g_z_0_yy_xxx[k];

            g_z_0_xyy_xy[k] = -g_z_0_yy_xy[k] * cd_x[k] + g_z_0_yy_xxy[k];

            g_z_0_xyy_xz[k] = -g_z_0_yy_xz[k] * cd_x[k] + g_z_0_yy_xxz[k];

            g_z_0_xyy_yy[k] = -g_z_0_yy_yy[k] * cd_x[k] + g_z_0_yy_xyy[k];

            g_z_0_xyy_yz[k] = -g_z_0_yy_yz[k] * cd_x[k] + g_z_0_yy_xyz[k];

            g_z_0_xyy_zz[k] = -g_z_0_yy_zz[k] * cd_x[k] + g_z_0_yy_xzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 24);

        auto g_z_0_xyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 25);

        auto g_z_0_xyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 26);

        auto g_z_0_xyz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 27);

        auto g_z_0_xyz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 28);

        auto g_z_0_xyz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 29);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_xx, g_z_0_xyz_xy, g_z_0_xyz_xz, g_z_0_xyz_yy, g_z_0_xyz_yz, g_z_0_xyz_zz, g_z_0_yz_xx, g_z_0_yz_xxx, g_z_0_yz_xxy, g_z_0_yz_xxz, g_z_0_yz_xy, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xz, g_z_0_yz_xzz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_xx[k] = -g_z_0_yz_xx[k] * cd_x[k] + g_z_0_yz_xxx[k];

            g_z_0_xyz_xy[k] = -g_z_0_yz_xy[k] * cd_x[k] + g_z_0_yz_xxy[k];

            g_z_0_xyz_xz[k] = -g_z_0_yz_xz[k] * cd_x[k] + g_z_0_yz_xxz[k];

            g_z_0_xyz_yy[k] = -g_z_0_yz_yy[k] * cd_x[k] + g_z_0_yz_xyy[k];

            g_z_0_xyz_yz[k] = -g_z_0_yz_yz[k] * cd_x[k] + g_z_0_yz_xyz[k];

            g_z_0_xyz_zz[k] = -g_z_0_yz_zz[k] * cd_x[k] + g_z_0_yz_xzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 30);

        auto g_z_0_xzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 31);

        auto g_z_0_xzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 32);

        auto g_z_0_xzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 33);

        auto g_z_0_xzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 34);

        auto g_z_0_xzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 35);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_xx, g_z_0_xzz_xy, g_z_0_xzz_xz, g_z_0_xzz_yy, g_z_0_xzz_yz, g_z_0_xzz_zz, g_z_0_zz_xx, g_z_0_zz_xxx, g_z_0_zz_xxy, g_z_0_zz_xxz, g_z_0_zz_xy, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_xzz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_xx[k] = -g_z_0_zz_xx[k] * cd_x[k] + g_z_0_zz_xxx[k];

            g_z_0_xzz_xy[k] = -g_z_0_zz_xy[k] * cd_x[k] + g_z_0_zz_xxy[k];

            g_z_0_xzz_xz[k] = -g_z_0_zz_xz[k] * cd_x[k] + g_z_0_zz_xxz[k];

            g_z_0_xzz_yy[k] = -g_z_0_zz_yy[k] * cd_x[k] + g_z_0_zz_xyy[k];

            g_z_0_xzz_yz[k] = -g_z_0_zz_yz[k] * cd_x[k] + g_z_0_zz_xyz[k];

            g_z_0_xzz_zz[k] = -g_z_0_zz_zz[k] * cd_x[k] + g_z_0_zz_xzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 36);

        auto g_z_0_yyy_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 37);

        auto g_z_0_yyy_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 38);

        auto g_z_0_yyy_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 39);

        auto g_z_0_yyy_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 40);

        auto g_z_0_yyy_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 41);

        #pragma omp simd aligned(cd_y, g_z_0_yy_xx, g_z_0_yy_xxy, g_z_0_yy_xy, g_z_0_yy_xyy, g_z_0_yy_xyz, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yyy, g_z_0_yy_yyz, g_z_0_yy_yz, g_z_0_yy_yzz, g_z_0_yy_zz, g_z_0_yyy_xx, g_z_0_yyy_xy, g_z_0_yyy_xz, g_z_0_yyy_yy, g_z_0_yyy_yz, g_z_0_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_xx[k] = -g_z_0_yy_xx[k] * cd_y[k] + g_z_0_yy_xxy[k];

            g_z_0_yyy_xy[k] = -g_z_0_yy_xy[k] * cd_y[k] + g_z_0_yy_xyy[k];

            g_z_0_yyy_xz[k] = -g_z_0_yy_xz[k] * cd_y[k] + g_z_0_yy_xyz[k];

            g_z_0_yyy_yy[k] = -g_z_0_yy_yy[k] * cd_y[k] + g_z_0_yy_yyy[k];

            g_z_0_yyy_yz[k] = -g_z_0_yy_yz[k] * cd_y[k] + g_z_0_yy_yyz[k];

            g_z_0_yyy_zz[k] = -g_z_0_yy_zz[k] * cd_y[k] + g_z_0_yy_yzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 42);

        auto g_z_0_yyz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 43);

        auto g_z_0_yyz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 44);

        auto g_z_0_yyz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 45);

        auto g_z_0_yyz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 46);

        auto g_z_0_yyz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 47);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_xx, g_z_0_yyz_xy, g_z_0_yyz_xz, g_z_0_yyz_yy, g_z_0_yyz_yz, g_z_0_yyz_zz, g_z_0_yz_xx, g_z_0_yz_xxy, g_z_0_yz_xy, g_z_0_yz_xyy, g_z_0_yz_xyz, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yyy, g_z_0_yz_yyz, g_z_0_yz_yz, g_z_0_yz_yzz, g_z_0_yz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_xx[k] = -g_z_0_yz_xx[k] * cd_y[k] + g_z_0_yz_xxy[k];

            g_z_0_yyz_xy[k] = -g_z_0_yz_xy[k] * cd_y[k] + g_z_0_yz_xyy[k];

            g_z_0_yyz_xz[k] = -g_z_0_yz_xz[k] * cd_y[k] + g_z_0_yz_xyz[k];

            g_z_0_yyz_yy[k] = -g_z_0_yz_yy[k] * cd_y[k] + g_z_0_yz_yyy[k];

            g_z_0_yyz_yz[k] = -g_z_0_yz_yz[k] * cd_y[k] + g_z_0_yz_yyz[k];

            g_z_0_yyz_zz[k] = -g_z_0_yz_zz[k] * cd_y[k] + g_z_0_yz_yzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 48);

        auto g_z_0_yzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 49);

        auto g_z_0_yzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 50);

        auto g_z_0_yzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 51);

        auto g_z_0_yzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 52);

        auto g_z_0_yzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 53);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_xx, g_z_0_yzz_xy, g_z_0_yzz_xz, g_z_0_yzz_yy, g_z_0_yzz_yz, g_z_0_yzz_zz, g_z_0_zz_xx, g_z_0_zz_xxy, g_z_0_zz_xy, g_z_0_zz_xyy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yyy, g_z_0_zz_yyz, g_z_0_zz_yz, g_z_0_zz_yzz, g_z_0_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_xx[k] = -g_z_0_zz_xx[k] * cd_y[k] + g_z_0_zz_xxy[k];

            g_z_0_yzz_xy[k] = -g_z_0_zz_xy[k] * cd_y[k] + g_z_0_zz_xyy[k];

            g_z_0_yzz_xz[k] = -g_z_0_zz_xz[k] * cd_y[k] + g_z_0_zz_xyz[k];

            g_z_0_yzz_yy[k] = -g_z_0_zz_yy[k] * cd_y[k] + g_z_0_zz_yyy[k];

            g_z_0_yzz_yz[k] = -g_z_0_zz_yz[k] * cd_y[k] + g_z_0_zz_yyz[k];

            g_z_0_yzz_zz[k] = -g_z_0_zz_zz[k] * cd_y[k] + g_z_0_zz_yzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_xx = cbuffer.data(fd_geom_10_off + 120 * acomps  + 54);

        auto g_z_0_zzz_xy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 55);

        auto g_z_0_zzz_xz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 56);

        auto g_z_0_zzz_yy = cbuffer.data(fd_geom_10_off + 120 * acomps  + 57);

        auto g_z_0_zzz_yz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 58);

        auto g_z_0_zzz_zz = cbuffer.data(fd_geom_10_off + 120 * acomps  + 59);

        #pragma omp simd aligned(cd_z, g_z_0_zz_xx, g_z_0_zz_xxz, g_z_0_zz_xy, g_z_0_zz_xyz, g_z_0_zz_xz, g_z_0_zz_xzz, g_z_0_zz_yy, g_z_0_zz_yyz, g_z_0_zz_yz, g_z_0_zz_yzz, g_z_0_zz_zz, g_z_0_zz_zzz, g_z_0_zzz_xx, g_z_0_zzz_xy, g_z_0_zzz_xz, g_z_0_zzz_yy, g_z_0_zzz_yz, g_z_0_zzz_zz, g_zz_xx, g_zz_xy, g_zz_xz, g_zz_yy, g_zz_yz, g_zz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_xx[k] = -g_zz_xx[k] - g_z_0_zz_xx[k] * cd_z[k] + g_z_0_zz_xxz[k];

            g_z_0_zzz_xy[k] = -g_zz_xy[k] - g_z_0_zz_xy[k] * cd_z[k] + g_z_0_zz_xyz[k];

            g_z_0_zzz_xz[k] = -g_zz_xz[k] - g_z_0_zz_xz[k] * cd_z[k] + g_z_0_zz_xzz[k];

            g_z_0_zzz_yy[k] = -g_zz_yy[k] - g_z_0_zz_yy[k] * cd_z[k] + g_z_0_zz_yyz[k];

            g_z_0_zzz_yz[k] = -g_zz_yz[k] - g_z_0_zz_yz[k] * cd_z[k] + g_z_0_zz_yzz[k];

            g_z_0_zzz_zz[k] = -g_zz_zz[k] - g_z_0_zz_zz[k] * cd_z[k] + g_z_0_zz_zzz[k];
        }
    }
}

} // t3ceri namespace

