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

#include "ThreeCenterElectronRepulsionContrRecXFF.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xff(CSimdArray<double>& cbuffer,
                                const size_t idx_xff,
                                const size_t idx_xdf,
                                const size_t idx_xdg,
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
        /// Set up components of auxilary buffer : SDF

        const auto df_off = idx_xdf + i * 60;

        auto g_xx_xxx = cbuffer.data(df_off + 0);

        auto g_xx_xxy = cbuffer.data(df_off + 1);

        auto g_xx_xxz = cbuffer.data(df_off + 2);

        auto g_xx_xyy = cbuffer.data(df_off + 3);

        auto g_xx_xyz = cbuffer.data(df_off + 4);

        auto g_xx_xzz = cbuffer.data(df_off + 5);

        auto g_xx_yyy = cbuffer.data(df_off + 6);

        auto g_xx_yyz = cbuffer.data(df_off + 7);

        auto g_xx_yzz = cbuffer.data(df_off + 8);

        auto g_xx_zzz = cbuffer.data(df_off + 9);

        auto g_xy_xxx = cbuffer.data(df_off + 10);

        auto g_xy_xxy = cbuffer.data(df_off + 11);

        auto g_xy_xxz = cbuffer.data(df_off + 12);

        auto g_xy_xyy = cbuffer.data(df_off + 13);

        auto g_xy_xyz = cbuffer.data(df_off + 14);

        auto g_xy_xzz = cbuffer.data(df_off + 15);

        auto g_xy_yyy = cbuffer.data(df_off + 16);

        auto g_xy_yyz = cbuffer.data(df_off + 17);

        auto g_xy_yzz = cbuffer.data(df_off + 18);

        auto g_xy_zzz = cbuffer.data(df_off + 19);

        auto g_xz_xxx = cbuffer.data(df_off + 20);

        auto g_xz_xxy = cbuffer.data(df_off + 21);

        auto g_xz_xxz = cbuffer.data(df_off + 22);

        auto g_xz_xyy = cbuffer.data(df_off + 23);

        auto g_xz_xyz = cbuffer.data(df_off + 24);

        auto g_xz_xzz = cbuffer.data(df_off + 25);

        auto g_xz_yyy = cbuffer.data(df_off + 26);

        auto g_xz_yyz = cbuffer.data(df_off + 27);

        auto g_xz_yzz = cbuffer.data(df_off + 28);

        auto g_xz_zzz = cbuffer.data(df_off + 29);

        auto g_yy_xxx = cbuffer.data(df_off + 30);

        auto g_yy_xxy = cbuffer.data(df_off + 31);

        auto g_yy_xxz = cbuffer.data(df_off + 32);

        auto g_yy_xyy = cbuffer.data(df_off + 33);

        auto g_yy_xyz = cbuffer.data(df_off + 34);

        auto g_yy_xzz = cbuffer.data(df_off + 35);

        auto g_yy_yyy = cbuffer.data(df_off + 36);

        auto g_yy_yyz = cbuffer.data(df_off + 37);

        auto g_yy_yzz = cbuffer.data(df_off + 38);

        auto g_yy_zzz = cbuffer.data(df_off + 39);

        auto g_yz_xxx = cbuffer.data(df_off + 40);

        auto g_yz_xxy = cbuffer.data(df_off + 41);

        auto g_yz_xxz = cbuffer.data(df_off + 42);

        auto g_yz_xyy = cbuffer.data(df_off + 43);

        auto g_yz_xyz = cbuffer.data(df_off + 44);

        auto g_yz_xzz = cbuffer.data(df_off + 45);

        auto g_yz_yyy = cbuffer.data(df_off + 46);

        auto g_yz_yyz = cbuffer.data(df_off + 47);

        auto g_yz_yzz = cbuffer.data(df_off + 48);

        auto g_yz_zzz = cbuffer.data(df_off + 49);

        auto g_zz_xxx = cbuffer.data(df_off + 50);

        auto g_zz_xxy = cbuffer.data(df_off + 51);

        auto g_zz_xxz = cbuffer.data(df_off + 52);

        auto g_zz_xyy = cbuffer.data(df_off + 53);

        auto g_zz_xyz = cbuffer.data(df_off + 54);

        auto g_zz_xzz = cbuffer.data(df_off + 55);

        auto g_zz_yyy = cbuffer.data(df_off + 56);

        auto g_zz_yyz = cbuffer.data(df_off + 57);

        auto g_zz_yzz = cbuffer.data(df_off + 58);

        auto g_zz_zzz = cbuffer.data(df_off + 59);

        /// Set up components of auxilary buffer : SDG

        const auto dg_off = idx_xdg + i * 90;

        auto g_xx_xxxx = cbuffer.data(dg_off + 0);

        auto g_xx_xxxy = cbuffer.data(dg_off + 1);

        auto g_xx_xxxz = cbuffer.data(dg_off + 2);

        auto g_xx_xxyy = cbuffer.data(dg_off + 3);

        auto g_xx_xxyz = cbuffer.data(dg_off + 4);

        auto g_xx_xxzz = cbuffer.data(dg_off + 5);

        auto g_xx_xyyy = cbuffer.data(dg_off + 6);

        auto g_xx_xyyz = cbuffer.data(dg_off + 7);

        auto g_xx_xyzz = cbuffer.data(dg_off + 8);

        auto g_xx_xzzz = cbuffer.data(dg_off + 9);

        auto g_xy_xxxx = cbuffer.data(dg_off + 15);

        auto g_xy_xxxy = cbuffer.data(dg_off + 16);

        auto g_xy_xxxz = cbuffer.data(dg_off + 17);

        auto g_xy_xxyy = cbuffer.data(dg_off + 18);

        auto g_xy_xxyz = cbuffer.data(dg_off + 19);

        auto g_xy_xxzz = cbuffer.data(dg_off + 20);

        auto g_xy_xyyy = cbuffer.data(dg_off + 21);

        auto g_xy_xyyz = cbuffer.data(dg_off + 22);

        auto g_xy_xyzz = cbuffer.data(dg_off + 23);

        auto g_xy_xzzz = cbuffer.data(dg_off + 24);

        auto g_xz_xxxx = cbuffer.data(dg_off + 30);

        auto g_xz_xxxy = cbuffer.data(dg_off + 31);

        auto g_xz_xxxz = cbuffer.data(dg_off + 32);

        auto g_xz_xxyy = cbuffer.data(dg_off + 33);

        auto g_xz_xxyz = cbuffer.data(dg_off + 34);

        auto g_xz_xxzz = cbuffer.data(dg_off + 35);

        auto g_xz_xyyy = cbuffer.data(dg_off + 36);

        auto g_xz_xyyz = cbuffer.data(dg_off + 37);

        auto g_xz_xyzz = cbuffer.data(dg_off + 38);

        auto g_xz_xzzz = cbuffer.data(dg_off + 39);

        auto g_yy_xxxx = cbuffer.data(dg_off + 45);

        auto g_yy_xxxy = cbuffer.data(dg_off + 46);

        auto g_yy_xxxz = cbuffer.data(dg_off + 47);

        auto g_yy_xxyy = cbuffer.data(dg_off + 48);

        auto g_yy_xxyz = cbuffer.data(dg_off + 49);

        auto g_yy_xxzz = cbuffer.data(dg_off + 50);

        auto g_yy_xyyy = cbuffer.data(dg_off + 51);

        auto g_yy_xyyz = cbuffer.data(dg_off + 52);

        auto g_yy_xyzz = cbuffer.data(dg_off + 53);

        auto g_yy_xzzz = cbuffer.data(dg_off + 54);

        auto g_yy_yyyy = cbuffer.data(dg_off + 55);

        auto g_yy_yyyz = cbuffer.data(dg_off + 56);

        auto g_yy_yyzz = cbuffer.data(dg_off + 57);

        auto g_yy_yzzz = cbuffer.data(dg_off + 58);

        auto g_yz_xxxx = cbuffer.data(dg_off + 60);

        auto g_yz_xxxy = cbuffer.data(dg_off + 61);

        auto g_yz_xxxz = cbuffer.data(dg_off + 62);

        auto g_yz_xxyy = cbuffer.data(dg_off + 63);

        auto g_yz_xxyz = cbuffer.data(dg_off + 64);

        auto g_yz_xxzz = cbuffer.data(dg_off + 65);

        auto g_yz_xyyy = cbuffer.data(dg_off + 66);

        auto g_yz_xyyz = cbuffer.data(dg_off + 67);

        auto g_yz_xyzz = cbuffer.data(dg_off + 68);

        auto g_yz_xzzz = cbuffer.data(dg_off + 69);

        auto g_yz_yyyy = cbuffer.data(dg_off + 70);

        auto g_yz_yyyz = cbuffer.data(dg_off + 71);

        auto g_yz_yyzz = cbuffer.data(dg_off + 72);

        auto g_yz_yzzz = cbuffer.data(dg_off + 73);

        auto g_zz_xxxx = cbuffer.data(dg_off + 75);

        auto g_zz_xxxy = cbuffer.data(dg_off + 76);

        auto g_zz_xxxz = cbuffer.data(dg_off + 77);

        auto g_zz_xxyy = cbuffer.data(dg_off + 78);

        auto g_zz_xxyz = cbuffer.data(dg_off + 79);

        auto g_zz_xxzz = cbuffer.data(dg_off + 80);

        auto g_zz_xyyy = cbuffer.data(dg_off + 81);

        auto g_zz_xyyz = cbuffer.data(dg_off + 82);

        auto g_zz_xyzz = cbuffer.data(dg_off + 83);

        auto g_zz_xzzz = cbuffer.data(dg_off + 84);

        auto g_zz_yyyy = cbuffer.data(dg_off + 85);

        auto g_zz_yyyz = cbuffer.data(dg_off + 86);

        auto g_zz_yyzz = cbuffer.data(dg_off + 87);

        auto g_zz_yzzz = cbuffer.data(dg_off + 88);

        auto g_zz_zzzz = cbuffer.data(dg_off + 89);

        /// set up bra offset for contr_buffer_xff

        const auto ff_off = idx_xff + i * 100;

        /// Set up 0-10 components of targeted buffer : cbuffer.data(

        auto g_xxx_xxx = cbuffer.data(ff_off + 0);

        auto g_xxx_xxy = cbuffer.data(ff_off + 1);

        auto g_xxx_xxz = cbuffer.data(ff_off + 2);

        auto g_xxx_xyy = cbuffer.data(ff_off + 3);

        auto g_xxx_xyz = cbuffer.data(ff_off + 4);

        auto g_xxx_xzz = cbuffer.data(ff_off + 5);

        auto g_xxx_yyy = cbuffer.data(ff_off + 6);

        auto g_xxx_yyz = cbuffer.data(ff_off + 7);

        auto g_xxx_yzz = cbuffer.data(ff_off + 8);

        auto g_xxx_zzz = cbuffer.data(ff_off + 9);

        #pragma omp simd aligned(cd_x, g_xx_xxx, g_xx_xxxx, g_xx_xxxy, g_xx_xxxz, g_xx_xxy, g_xx_xxyy, g_xx_xxyz, g_xx_xxz, g_xx_xxzz, g_xx_xyy, g_xx_xyyy, g_xx_xyyz, g_xx_xyz, g_xx_xyzz, g_xx_xzz, g_xx_xzzz, g_xx_yyy, g_xx_yyz, g_xx_yzz, g_xx_zzz, g_xxx_xxx, g_xxx_xxy, g_xxx_xxz, g_xxx_xyy, g_xxx_xyz, g_xxx_xzz, g_xxx_yyy, g_xxx_yyz, g_xxx_yzz, g_xxx_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxx_xxx[k] = -g_xx_xxx[k] * cd_x[k] + g_xx_xxxx[k];

            g_xxx_xxy[k] = -g_xx_xxy[k] * cd_x[k] + g_xx_xxxy[k];

            g_xxx_xxz[k] = -g_xx_xxz[k] * cd_x[k] + g_xx_xxxz[k];

            g_xxx_xyy[k] = -g_xx_xyy[k] * cd_x[k] + g_xx_xxyy[k];

            g_xxx_xyz[k] = -g_xx_xyz[k] * cd_x[k] + g_xx_xxyz[k];

            g_xxx_xzz[k] = -g_xx_xzz[k] * cd_x[k] + g_xx_xxzz[k];

            g_xxx_yyy[k] = -g_xx_yyy[k] * cd_x[k] + g_xx_xyyy[k];

            g_xxx_yyz[k] = -g_xx_yyz[k] * cd_x[k] + g_xx_xyyz[k];

            g_xxx_yzz[k] = -g_xx_yzz[k] * cd_x[k] + g_xx_xyzz[k];

            g_xxx_zzz[k] = -g_xx_zzz[k] * cd_x[k] + g_xx_xzzz[k];
        }

        /// Set up 10-20 components of targeted buffer : cbuffer.data(

        auto g_xxy_xxx = cbuffer.data(ff_off + 10);

        auto g_xxy_xxy = cbuffer.data(ff_off + 11);

        auto g_xxy_xxz = cbuffer.data(ff_off + 12);

        auto g_xxy_xyy = cbuffer.data(ff_off + 13);

        auto g_xxy_xyz = cbuffer.data(ff_off + 14);

        auto g_xxy_xzz = cbuffer.data(ff_off + 15);

        auto g_xxy_yyy = cbuffer.data(ff_off + 16);

        auto g_xxy_yyz = cbuffer.data(ff_off + 17);

        auto g_xxy_yzz = cbuffer.data(ff_off + 18);

        auto g_xxy_zzz = cbuffer.data(ff_off + 19);

        #pragma omp simd aligned(cd_x, g_xxy_xxx, g_xxy_xxy, g_xxy_xxz, g_xxy_xyy, g_xxy_xyz, g_xxy_xzz, g_xxy_yyy, g_xxy_yyz, g_xxy_yzz, g_xxy_zzz, g_xy_xxx, g_xy_xxxx, g_xy_xxxy, g_xy_xxxz, g_xy_xxy, g_xy_xxyy, g_xy_xxyz, g_xy_xxz, g_xy_xxzz, g_xy_xyy, g_xy_xyyy, g_xy_xyyz, g_xy_xyz, g_xy_xyzz, g_xy_xzz, g_xy_xzzz, g_xy_yyy, g_xy_yyz, g_xy_yzz, g_xy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxy_xxx[k] = -g_xy_xxx[k] * cd_x[k] + g_xy_xxxx[k];

            g_xxy_xxy[k] = -g_xy_xxy[k] * cd_x[k] + g_xy_xxxy[k];

            g_xxy_xxz[k] = -g_xy_xxz[k] * cd_x[k] + g_xy_xxxz[k];

            g_xxy_xyy[k] = -g_xy_xyy[k] * cd_x[k] + g_xy_xxyy[k];

            g_xxy_xyz[k] = -g_xy_xyz[k] * cd_x[k] + g_xy_xxyz[k];

            g_xxy_xzz[k] = -g_xy_xzz[k] * cd_x[k] + g_xy_xxzz[k];

            g_xxy_yyy[k] = -g_xy_yyy[k] * cd_x[k] + g_xy_xyyy[k];

            g_xxy_yyz[k] = -g_xy_yyz[k] * cd_x[k] + g_xy_xyyz[k];

            g_xxy_yzz[k] = -g_xy_yzz[k] * cd_x[k] + g_xy_xyzz[k];

            g_xxy_zzz[k] = -g_xy_zzz[k] * cd_x[k] + g_xy_xzzz[k];
        }

        /// Set up 20-30 components of targeted buffer : cbuffer.data(

        auto g_xxz_xxx = cbuffer.data(ff_off + 20);

        auto g_xxz_xxy = cbuffer.data(ff_off + 21);

        auto g_xxz_xxz = cbuffer.data(ff_off + 22);

        auto g_xxz_xyy = cbuffer.data(ff_off + 23);

        auto g_xxz_xyz = cbuffer.data(ff_off + 24);

        auto g_xxz_xzz = cbuffer.data(ff_off + 25);

        auto g_xxz_yyy = cbuffer.data(ff_off + 26);

        auto g_xxz_yyz = cbuffer.data(ff_off + 27);

        auto g_xxz_yzz = cbuffer.data(ff_off + 28);

        auto g_xxz_zzz = cbuffer.data(ff_off + 29);

        #pragma omp simd aligned(cd_x, g_xxz_xxx, g_xxz_xxy, g_xxz_xxz, g_xxz_xyy, g_xxz_xyz, g_xxz_xzz, g_xxz_yyy, g_xxz_yyz, g_xxz_yzz, g_xxz_zzz, g_xz_xxx, g_xz_xxxx, g_xz_xxxy, g_xz_xxxz, g_xz_xxy, g_xz_xxyy, g_xz_xxyz, g_xz_xxz, g_xz_xxzz, g_xz_xyy, g_xz_xyyy, g_xz_xyyz, g_xz_xyz, g_xz_xyzz, g_xz_xzz, g_xz_xzzz, g_xz_yyy, g_xz_yyz, g_xz_yzz, g_xz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxz_xxx[k] = -g_xz_xxx[k] * cd_x[k] + g_xz_xxxx[k];

            g_xxz_xxy[k] = -g_xz_xxy[k] * cd_x[k] + g_xz_xxxy[k];

            g_xxz_xxz[k] = -g_xz_xxz[k] * cd_x[k] + g_xz_xxxz[k];

            g_xxz_xyy[k] = -g_xz_xyy[k] * cd_x[k] + g_xz_xxyy[k];

            g_xxz_xyz[k] = -g_xz_xyz[k] * cd_x[k] + g_xz_xxyz[k];

            g_xxz_xzz[k] = -g_xz_xzz[k] * cd_x[k] + g_xz_xxzz[k];

            g_xxz_yyy[k] = -g_xz_yyy[k] * cd_x[k] + g_xz_xyyy[k];

            g_xxz_yyz[k] = -g_xz_yyz[k] * cd_x[k] + g_xz_xyyz[k];

            g_xxz_yzz[k] = -g_xz_yzz[k] * cd_x[k] + g_xz_xyzz[k];

            g_xxz_zzz[k] = -g_xz_zzz[k] * cd_x[k] + g_xz_xzzz[k];
        }

        /// Set up 30-40 components of targeted buffer : cbuffer.data(

        auto g_xyy_xxx = cbuffer.data(ff_off + 30);

        auto g_xyy_xxy = cbuffer.data(ff_off + 31);

        auto g_xyy_xxz = cbuffer.data(ff_off + 32);

        auto g_xyy_xyy = cbuffer.data(ff_off + 33);

        auto g_xyy_xyz = cbuffer.data(ff_off + 34);

        auto g_xyy_xzz = cbuffer.data(ff_off + 35);

        auto g_xyy_yyy = cbuffer.data(ff_off + 36);

        auto g_xyy_yyz = cbuffer.data(ff_off + 37);

        auto g_xyy_yzz = cbuffer.data(ff_off + 38);

        auto g_xyy_zzz = cbuffer.data(ff_off + 39);

        #pragma omp simd aligned(cd_x, g_xyy_xxx, g_xyy_xxy, g_xyy_xxz, g_xyy_xyy, g_xyy_xyz, g_xyy_xzz, g_xyy_yyy, g_xyy_yyz, g_xyy_yzz, g_xyy_zzz, g_yy_xxx, g_yy_xxxx, g_yy_xxxy, g_yy_xxxz, g_yy_xxy, g_yy_xxyy, g_yy_xxyz, g_yy_xxz, g_yy_xxzz, g_yy_xyy, g_yy_xyyy, g_yy_xyyz, g_yy_xyz, g_yy_xyzz, g_yy_xzz, g_yy_xzzz, g_yy_yyy, g_yy_yyz, g_yy_yzz, g_yy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyy_xxx[k] = -g_yy_xxx[k] * cd_x[k] + g_yy_xxxx[k];

            g_xyy_xxy[k] = -g_yy_xxy[k] * cd_x[k] + g_yy_xxxy[k];

            g_xyy_xxz[k] = -g_yy_xxz[k] * cd_x[k] + g_yy_xxxz[k];

            g_xyy_xyy[k] = -g_yy_xyy[k] * cd_x[k] + g_yy_xxyy[k];

            g_xyy_xyz[k] = -g_yy_xyz[k] * cd_x[k] + g_yy_xxyz[k];

            g_xyy_xzz[k] = -g_yy_xzz[k] * cd_x[k] + g_yy_xxzz[k];

            g_xyy_yyy[k] = -g_yy_yyy[k] * cd_x[k] + g_yy_xyyy[k];

            g_xyy_yyz[k] = -g_yy_yyz[k] * cd_x[k] + g_yy_xyyz[k];

            g_xyy_yzz[k] = -g_yy_yzz[k] * cd_x[k] + g_yy_xyzz[k];

            g_xyy_zzz[k] = -g_yy_zzz[k] * cd_x[k] + g_yy_xzzz[k];
        }

        /// Set up 40-50 components of targeted buffer : cbuffer.data(

        auto g_xyz_xxx = cbuffer.data(ff_off + 40);

        auto g_xyz_xxy = cbuffer.data(ff_off + 41);

        auto g_xyz_xxz = cbuffer.data(ff_off + 42);

        auto g_xyz_xyy = cbuffer.data(ff_off + 43);

        auto g_xyz_xyz = cbuffer.data(ff_off + 44);

        auto g_xyz_xzz = cbuffer.data(ff_off + 45);

        auto g_xyz_yyy = cbuffer.data(ff_off + 46);

        auto g_xyz_yyz = cbuffer.data(ff_off + 47);

        auto g_xyz_yzz = cbuffer.data(ff_off + 48);

        auto g_xyz_zzz = cbuffer.data(ff_off + 49);

        #pragma omp simd aligned(cd_x, g_xyz_xxx, g_xyz_xxy, g_xyz_xxz, g_xyz_xyy, g_xyz_xyz, g_xyz_xzz, g_xyz_yyy, g_xyz_yyz, g_xyz_yzz, g_xyz_zzz, g_yz_xxx, g_yz_xxxx, g_yz_xxxy, g_yz_xxxz, g_yz_xxy, g_yz_xxyy, g_yz_xxyz, g_yz_xxz, g_yz_xxzz, g_yz_xyy, g_yz_xyyy, g_yz_xyyz, g_yz_xyz, g_yz_xyzz, g_yz_xzz, g_yz_xzzz, g_yz_yyy, g_yz_yyz, g_yz_yzz, g_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyz_xxx[k] = -g_yz_xxx[k] * cd_x[k] + g_yz_xxxx[k];

            g_xyz_xxy[k] = -g_yz_xxy[k] * cd_x[k] + g_yz_xxxy[k];

            g_xyz_xxz[k] = -g_yz_xxz[k] * cd_x[k] + g_yz_xxxz[k];

            g_xyz_xyy[k] = -g_yz_xyy[k] * cd_x[k] + g_yz_xxyy[k];

            g_xyz_xyz[k] = -g_yz_xyz[k] * cd_x[k] + g_yz_xxyz[k];

            g_xyz_xzz[k] = -g_yz_xzz[k] * cd_x[k] + g_yz_xxzz[k];

            g_xyz_yyy[k] = -g_yz_yyy[k] * cd_x[k] + g_yz_xyyy[k];

            g_xyz_yyz[k] = -g_yz_yyz[k] * cd_x[k] + g_yz_xyyz[k];

            g_xyz_yzz[k] = -g_yz_yzz[k] * cd_x[k] + g_yz_xyzz[k];

            g_xyz_zzz[k] = -g_yz_zzz[k] * cd_x[k] + g_yz_xzzz[k];
        }

        /// Set up 50-60 components of targeted buffer : cbuffer.data(

        auto g_xzz_xxx = cbuffer.data(ff_off + 50);

        auto g_xzz_xxy = cbuffer.data(ff_off + 51);

        auto g_xzz_xxz = cbuffer.data(ff_off + 52);

        auto g_xzz_xyy = cbuffer.data(ff_off + 53);

        auto g_xzz_xyz = cbuffer.data(ff_off + 54);

        auto g_xzz_xzz = cbuffer.data(ff_off + 55);

        auto g_xzz_yyy = cbuffer.data(ff_off + 56);

        auto g_xzz_yyz = cbuffer.data(ff_off + 57);

        auto g_xzz_yzz = cbuffer.data(ff_off + 58);

        auto g_xzz_zzz = cbuffer.data(ff_off + 59);

        #pragma omp simd aligned(cd_x, g_xzz_xxx, g_xzz_xxy, g_xzz_xxz, g_xzz_xyy, g_xzz_xyz, g_xzz_xzz, g_xzz_yyy, g_xzz_yyz, g_xzz_yzz, g_xzz_zzz, g_zz_xxx, g_zz_xxxx, g_zz_xxxy, g_zz_xxxz, g_zz_xxy, g_zz_xxyy, g_zz_xxyz, g_zz_xxz, g_zz_xxzz, g_zz_xyy, g_zz_xyyy, g_zz_xyyz, g_zz_xyz, g_zz_xyzz, g_zz_xzz, g_zz_xzzz, g_zz_yyy, g_zz_yyz, g_zz_yzz, g_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzz_xxx[k] = -g_zz_xxx[k] * cd_x[k] + g_zz_xxxx[k];

            g_xzz_xxy[k] = -g_zz_xxy[k] * cd_x[k] + g_zz_xxxy[k];

            g_xzz_xxz[k] = -g_zz_xxz[k] * cd_x[k] + g_zz_xxxz[k];

            g_xzz_xyy[k] = -g_zz_xyy[k] * cd_x[k] + g_zz_xxyy[k];

            g_xzz_xyz[k] = -g_zz_xyz[k] * cd_x[k] + g_zz_xxyz[k];

            g_xzz_xzz[k] = -g_zz_xzz[k] * cd_x[k] + g_zz_xxzz[k];

            g_xzz_yyy[k] = -g_zz_yyy[k] * cd_x[k] + g_zz_xyyy[k];

            g_xzz_yyz[k] = -g_zz_yyz[k] * cd_x[k] + g_zz_xyyz[k];

            g_xzz_yzz[k] = -g_zz_yzz[k] * cd_x[k] + g_zz_xyzz[k];

            g_xzz_zzz[k] = -g_zz_zzz[k] * cd_x[k] + g_zz_xzzz[k];
        }

        /// Set up 60-70 components of targeted buffer : cbuffer.data(

        auto g_yyy_xxx = cbuffer.data(ff_off + 60);

        auto g_yyy_xxy = cbuffer.data(ff_off + 61);

        auto g_yyy_xxz = cbuffer.data(ff_off + 62);

        auto g_yyy_xyy = cbuffer.data(ff_off + 63);

        auto g_yyy_xyz = cbuffer.data(ff_off + 64);

        auto g_yyy_xzz = cbuffer.data(ff_off + 65);

        auto g_yyy_yyy = cbuffer.data(ff_off + 66);

        auto g_yyy_yyz = cbuffer.data(ff_off + 67);

        auto g_yyy_yzz = cbuffer.data(ff_off + 68);

        auto g_yyy_zzz = cbuffer.data(ff_off + 69);

        #pragma omp simd aligned(cd_y, g_yy_xxx, g_yy_xxxy, g_yy_xxy, g_yy_xxyy, g_yy_xxyz, g_yy_xxz, g_yy_xyy, g_yy_xyyy, g_yy_xyyz, g_yy_xyz, g_yy_xyzz, g_yy_xzz, g_yy_yyy, g_yy_yyyy, g_yy_yyyz, g_yy_yyz, g_yy_yyzz, g_yy_yzz, g_yy_yzzz, g_yy_zzz, g_yyy_xxx, g_yyy_xxy, g_yyy_xxz, g_yyy_xyy, g_yyy_xyz, g_yyy_xzz, g_yyy_yyy, g_yyy_yyz, g_yyy_yzz, g_yyy_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyy_xxx[k] = -g_yy_xxx[k] * cd_y[k] + g_yy_xxxy[k];

            g_yyy_xxy[k] = -g_yy_xxy[k] * cd_y[k] + g_yy_xxyy[k];

            g_yyy_xxz[k] = -g_yy_xxz[k] * cd_y[k] + g_yy_xxyz[k];

            g_yyy_xyy[k] = -g_yy_xyy[k] * cd_y[k] + g_yy_xyyy[k];

            g_yyy_xyz[k] = -g_yy_xyz[k] * cd_y[k] + g_yy_xyyz[k];

            g_yyy_xzz[k] = -g_yy_xzz[k] * cd_y[k] + g_yy_xyzz[k];

            g_yyy_yyy[k] = -g_yy_yyy[k] * cd_y[k] + g_yy_yyyy[k];

            g_yyy_yyz[k] = -g_yy_yyz[k] * cd_y[k] + g_yy_yyyz[k];

            g_yyy_yzz[k] = -g_yy_yzz[k] * cd_y[k] + g_yy_yyzz[k];

            g_yyy_zzz[k] = -g_yy_zzz[k] * cd_y[k] + g_yy_yzzz[k];
        }

        /// Set up 70-80 components of targeted buffer : cbuffer.data(

        auto g_yyz_xxx = cbuffer.data(ff_off + 70);

        auto g_yyz_xxy = cbuffer.data(ff_off + 71);

        auto g_yyz_xxz = cbuffer.data(ff_off + 72);

        auto g_yyz_xyy = cbuffer.data(ff_off + 73);

        auto g_yyz_xyz = cbuffer.data(ff_off + 74);

        auto g_yyz_xzz = cbuffer.data(ff_off + 75);

        auto g_yyz_yyy = cbuffer.data(ff_off + 76);

        auto g_yyz_yyz = cbuffer.data(ff_off + 77);

        auto g_yyz_yzz = cbuffer.data(ff_off + 78);

        auto g_yyz_zzz = cbuffer.data(ff_off + 79);

        #pragma omp simd aligned(cd_y, g_yyz_xxx, g_yyz_xxy, g_yyz_xxz, g_yyz_xyy, g_yyz_xyz, g_yyz_xzz, g_yyz_yyy, g_yyz_yyz, g_yyz_yzz, g_yyz_zzz, g_yz_xxx, g_yz_xxxy, g_yz_xxy, g_yz_xxyy, g_yz_xxyz, g_yz_xxz, g_yz_xyy, g_yz_xyyy, g_yz_xyyz, g_yz_xyz, g_yz_xyzz, g_yz_xzz, g_yz_yyy, g_yz_yyyy, g_yz_yyyz, g_yz_yyz, g_yz_yyzz, g_yz_yzz, g_yz_yzzz, g_yz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyz_xxx[k] = -g_yz_xxx[k] * cd_y[k] + g_yz_xxxy[k];

            g_yyz_xxy[k] = -g_yz_xxy[k] * cd_y[k] + g_yz_xxyy[k];

            g_yyz_xxz[k] = -g_yz_xxz[k] * cd_y[k] + g_yz_xxyz[k];

            g_yyz_xyy[k] = -g_yz_xyy[k] * cd_y[k] + g_yz_xyyy[k];

            g_yyz_xyz[k] = -g_yz_xyz[k] * cd_y[k] + g_yz_xyyz[k];

            g_yyz_xzz[k] = -g_yz_xzz[k] * cd_y[k] + g_yz_xyzz[k];

            g_yyz_yyy[k] = -g_yz_yyy[k] * cd_y[k] + g_yz_yyyy[k];

            g_yyz_yyz[k] = -g_yz_yyz[k] * cd_y[k] + g_yz_yyyz[k];

            g_yyz_yzz[k] = -g_yz_yzz[k] * cd_y[k] + g_yz_yyzz[k];

            g_yyz_zzz[k] = -g_yz_zzz[k] * cd_y[k] + g_yz_yzzz[k];
        }

        /// Set up 80-90 components of targeted buffer : cbuffer.data(

        auto g_yzz_xxx = cbuffer.data(ff_off + 80);

        auto g_yzz_xxy = cbuffer.data(ff_off + 81);

        auto g_yzz_xxz = cbuffer.data(ff_off + 82);

        auto g_yzz_xyy = cbuffer.data(ff_off + 83);

        auto g_yzz_xyz = cbuffer.data(ff_off + 84);

        auto g_yzz_xzz = cbuffer.data(ff_off + 85);

        auto g_yzz_yyy = cbuffer.data(ff_off + 86);

        auto g_yzz_yyz = cbuffer.data(ff_off + 87);

        auto g_yzz_yzz = cbuffer.data(ff_off + 88);

        auto g_yzz_zzz = cbuffer.data(ff_off + 89);

        #pragma omp simd aligned(cd_y, g_yzz_xxx, g_yzz_xxy, g_yzz_xxz, g_yzz_xyy, g_yzz_xyz, g_yzz_xzz, g_yzz_yyy, g_yzz_yyz, g_yzz_yzz, g_yzz_zzz, g_zz_xxx, g_zz_xxxy, g_zz_xxy, g_zz_xxyy, g_zz_xxyz, g_zz_xxz, g_zz_xyy, g_zz_xyyy, g_zz_xyyz, g_zz_xyz, g_zz_xyzz, g_zz_xzz, g_zz_yyy, g_zz_yyyy, g_zz_yyyz, g_zz_yyz, g_zz_yyzz, g_zz_yzz, g_zz_yzzz, g_zz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzz_xxx[k] = -g_zz_xxx[k] * cd_y[k] + g_zz_xxxy[k];

            g_yzz_xxy[k] = -g_zz_xxy[k] * cd_y[k] + g_zz_xxyy[k];

            g_yzz_xxz[k] = -g_zz_xxz[k] * cd_y[k] + g_zz_xxyz[k];

            g_yzz_xyy[k] = -g_zz_xyy[k] * cd_y[k] + g_zz_xyyy[k];

            g_yzz_xyz[k] = -g_zz_xyz[k] * cd_y[k] + g_zz_xyyz[k];

            g_yzz_xzz[k] = -g_zz_xzz[k] * cd_y[k] + g_zz_xyzz[k];

            g_yzz_yyy[k] = -g_zz_yyy[k] * cd_y[k] + g_zz_yyyy[k];

            g_yzz_yyz[k] = -g_zz_yyz[k] * cd_y[k] + g_zz_yyyz[k];

            g_yzz_yzz[k] = -g_zz_yzz[k] * cd_y[k] + g_zz_yyzz[k];

            g_yzz_zzz[k] = -g_zz_zzz[k] * cd_y[k] + g_zz_yzzz[k];
        }

        /// Set up 90-100 components of targeted buffer : cbuffer.data(

        auto g_zzz_xxx = cbuffer.data(ff_off + 90);

        auto g_zzz_xxy = cbuffer.data(ff_off + 91);

        auto g_zzz_xxz = cbuffer.data(ff_off + 92);

        auto g_zzz_xyy = cbuffer.data(ff_off + 93);

        auto g_zzz_xyz = cbuffer.data(ff_off + 94);

        auto g_zzz_xzz = cbuffer.data(ff_off + 95);

        auto g_zzz_yyy = cbuffer.data(ff_off + 96);

        auto g_zzz_yyz = cbuffer.data(ff_off + 97);

        auto g_zzz_yzz = cbuffer.data(ff_off + 98);

        auto g_zzz_zzz = cbuffer.data(ff_off + 99);

        #pragma omp simd aligned(cd_z, g_zz_xxx, g_zz_xxxz, g_zz_xxy, g_zz_xxyz, g_zz_xxz, g_zz_xxzz, g_zz_xyy, g_zz_xyyz, g_zz_xyz, g_zz_xyzz, g_zz_xzz, g_zz_xzzz, g_zz_yyy, g_zz_yyyz, g_zz_yyz, g_zz_yyzz, g_zz_yzz, g_zz_yzzz, g_zz_zzz, g_zz_zzzz, g_zzz_xxx, g_zzz_xxy, g_zzz_xxz, g_zzz_xyy, g_zzz_xyz, g_zzz_xzz, g_zzz_yyy, g_zzz_yyz, g_zzz_yzz, g_zzz_zzz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzz_xxx[k] = -g_zz_xxx[k] * cd_z[k] + g_zz_xxxz[k];

            g_zzz_xxy[k] = -g_zz_xxy[k] * cd_z[k] + g_zz_xxyz[k];

            g_zzz_xxz[k] = -g_zz_xxz[k] * cd_z[k] + g_zz_xxzz[k];

            g_zzz_xyy[k] = -g_zz_xyy[k] * cd_z[k] + g_zz_xyyz[k];

            g_zzz_xyz[k] = -g_zz_xyz[k] * cd_z[k] + g_zz_xyzz[k];

            g_zzz_xzz[k] = -g_zz_xzz[k] * cd_z[k] + g_zz_xzzz[k];

            g_zzz_yyy[k] = -g_zz_yyy[k] * cd_z[k] + g_zz_yyyz[k];

            g_zzz_yyz[k] = -g_zz_yyz[k] * cd_z[k] + g_zz_yyzz[k];

            g_zzz_yzz[k] = -g_zz_yzz[k] * cd_z[k] + g_zz_yzzz[k];

            g_zzz_zzz[k] = -g_zz_zzz[k] * cd_z[k] + g_zz_zzzz[k];
        }
    }
}

} // t3ceri namespace

