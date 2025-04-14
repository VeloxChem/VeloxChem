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

#include "ThreeCenterElectronRepulsionContrRecXGD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xgd(CSimdArray<double>& cbuffer,
                                const size_t idx_xgd,
                                const size_t idx_xfd,
                                const size_t idx_xff,
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
        /// Set up components of auxilary buffer : SFD

        const auto fd_off = idx_xfd + i * 60;

        auto g_xxx_xx = cbuffer.data(fd_off + 0);

        auto g_xxx_xy = cbuffer.data(fd_off + 1);

        auto g_xxx_xz = cbuffer.data(fd_off + 2);

        auto g_xxx_yy = cbuffer.data(fd_off + 3);

        auto g_xxx_yz = cbuffer.data(fd_off + 4);

        auto g_xxx_zz = cbuffer.data(fd_off + 5);

        auto g_xxy_xx = cbuffer.data(fd_off + 6);

        auto g_xxy_xy = cbuffer.data(fd_off + 7);

        auto g_xxy_xz = cbuffer.data(fd_off + 8);

        auto g_xxy_yy = cbuffer.data(fd_off + 9);

        auto g_xxy_yz = cbuffer.data(fd_off + 10);

        auto g_xxy_zz = cbuffer.data(fd_off + 11);

        auto g_xxz_xx = cbuffer.data(fd_off + 12);

        auto g_xxz_xy = cbuffer.data(fd_off + 13);

        auto g_xxz_xz = cbuffer.data(fd_off + 14);

        auto g_xxz_yy = cbuffer.data(fd_off + 15);

        auto g_xxz_yz = cbuffer.data(fd_off + 16);

        auto g_xxz_zz = cbuffer.data(fd_off + 17);

        auto g_xyy_xx = cbuffer.data(fd_off + 18);

        auto g_xyy_xy = cbuffer.data(fd_off + 19);

        auto g_xyy_xz = cbuffer.data(fd_off + 20);

        auto g_xyy_yy = cbuffer.data(fd_off + 21);

        auto g_xyy_yz = cbuffer.data(fd_off + 22);

        auto g_xyy_zz = cbuffer.data(fd_off + 23);

        auto g_xyz_xx = cbuffer.data(fd_off + 24);

        auto g_xyz_xy = cbuffer.data(fd_off + 25);

        auto g_xyz_xz = cbuffer.data(fd_off + 26);

        auto g_xyz_yy = cbuffer.data(fd_off + 27);

        auto g_xyz_yz = cbuffer.data(fd_off + 28);

        auto g_xyz_zz = cbuffer.data(fd_off + 29);

        auto g_xzz_xx = cbuffer.data(fd_off + 30);

        auto g_xzz_xy = cbuffer.data(fd_off + 31);

        auto g_xzz_xz = cbuffer.data(fd_off + 32);

        auto g_xzz_yy = cbuffer.data(fd_off + 33);

        auto g_xzz_yz = cbuffer.data(fd_off + 34);

        auto g_xzz_zz = cbuffer.data(fd_off + 35);

        auto g_yyy_xx = cbuffer.data(fd_off + 36);

        auto g_yyy_xy = cbuffer.data(fd_off + 37);

        auto g_yyy_xz = cbuffer.data(fd_off + 38);

        auto g_yyy_yy = cbuffer.data(fd_off + 39);

        auto g_yyy_yz = cbuffer.data(fd_off + 40);

        auto g_yyy_zz = cbuffer.data(fd_off + 41);

        auto g_yyz_xx = cbuffer.data(fd_off + 42);

        auto g_yyz_xy = cbuffer.data(fd_off + 43);

        auto g_yyz_xz = cbuffer.data(fd_off + 44);

        auto g_yyz_yy = cbuffer.data(fd_off + 45);

        auto g_yyz_yz = cbuffer.data(fd_off + 46);

        auto g_yyz_zz = cbuffer.data(fd_off + 47);

        auto g_yzz_xx = cbuffer.data(fd_off + 48);

        auto g_yzz_xy = cbuffer.data(fd_off + 49);

        auto g_yzz_xz = cbuffer.data(fd_off + 50);

        auto g_yzz_yy = cbuffer.data(fd_off + 51);

        auto g_yzz_yz = cbuffer.data(fd_off + 52);

        auto g_yzz_zz = cbuffer.data(fd_off + 53);

        auto g_zzz_xx = cbuffer.data(fd_off + 54);

        auto g_zzz_xy = cbuffer.data(fd_off + 55);

        auto g_zzz_xz = cbuffer.data(fd_off + 56);

        auto g_zzz_yy = cbuffer.data(fd_off + 57);

        auto g_zzz_yz = cbuffer.data(fd_off + 58);

        auto g_zzz_zz = cbuffer.data(fd_off + 59);

        /// Set up components of auxilary buffer : SFF

        const auto ff_off = idx_xff + i * 100;

        auto g_xxx_xxx = cbuffer.data(ff_off + 0);

        auto g_xxx_xxy = cbuffer.data(ff_off + 1);

        auto g_xxx_xxz = cbuffer.data(ff_off + 2);

        auto g_xxx_xyy = cbuffer.data(ff_off + 3);

        auto g_xxx_xyz = cbuffer.data(ff_off + 4);

        auto g_xxx_xzz = cbuffer.data(ff_off + 5);

        auto g_xxy_xxx = cbuffer.data(ff_off + 10);

        auto g_xxy_xxy = cbuffer.data(ff_off + 11);

        auto g_xxy_xxz = cbuffer.data(ff_off + 12);

        auto g_xxy_xyy = cbuffer.data(ff_off + 13);

        auto g_xxy_xyz = cbuffer.data(ff_off + 14);

        auto g_xxy_xzz = cbuffer.data(ff_off + 15);

        auto g_xxz_xxx = cbuffer.data(ff_off + 20);

        auto g_xxz_xxy = cbuffer.data(ff_off + 21);

        auto g_xxz_xxz = cbuffer.data(ff_off + 22);

        auto g_xxz_xyy = cbuffer.data(ff_off + 23);

        auto g_xxz_xyz = cbuffer.data(ff_off + 24);

        auto g_xxz_xzz = cbuffer.data(ff_off + 25);

        auto g_xyy_xxx = cbuffer.data(ff_off + 30);

        auto g_xyy_xxy = cbuffer.data(ff_off + 31);

        auto g_xyy_xxz = cbuffer.data(ff_off + 32);

        auto g_xyy_xyy = cbuffer.data(ff_off + 33);

        auto g_xyy_xyz = cbuffer.data(ff_off + 34);

        auto g_xyy_xzz = cbuffer.data(ff_off + 35);

        auto g_xyz_xxx = cbuffer.data(ff_off + 40);

        auto g_xyz_xxy = cbuffer.data(ff_off + 41);

        auto g_xyz_xxz = cbuffer.data(ff_off + 42);

        auto g_xyz_xyy = cbuffer.data(ff_off + 43);

        auto g_xyz_xyz = cbuffer.data(ff_off + 44);

        auto g_xyz_xzz = cbuffer.data(ff_off + 45);

        auto g_xzz_xxx = cbuffer.data(ff_off + 50);

        auto g_xzz_xxy = cbuffer.data(ff_off + 51);

        auto g_xzz_xxz = cbuffer.data(ff_off + 52);

        auto g_xzz_xyy = cbuffer.data(ff_off + 53);

        auto g_xzz_xyz = cbuffer.data(ff_off + 54);

        auto g_xzz_xzz = cbuffer.data(ff_off + 55);

        auto g_yyy_xxx = cbuffer.data(ff_off + 60);

        auto g_yyy_xxy = cbuffer.data(ff_off + 61);

        auto g_yyy_xxz = cbuffer.data(ff_off + 62);

        auto g_yyy_xyy = cbuffer.data(ff_off + 63);

        auto g_yyy_xyz = cbuffer.data(ff_off + 64);

        auto g_yyy_xzz = cbuffer.data(ff_off + 65);

        auto g_yyy_yyy = cbuffer.data(ff_off + 66);

        auto g_yyy_yyz = cbuffer.data(ff_off + 67);

        auto g_yyy_yzz = cbuffer.data(ff_off + 68);

        auto g_yyz_xxx = cbuffer.data(ff_off + 70);

        auto g_yyz_xxy = cbuffer.data(ff_off + 71);

        auto g_yyz_xxz = cbuffer.data(ff_off + 72);

        auto g_yyz_xyy = cbuffer.data(ff_off + 73);

        auto g_yyz_xyz = cbuffer.data(ff_off + 74);

        auto g_yyz_xzz = cbuffer.data(ff_off + 75);

        auto g_yyz_yyy = cbuffer.data(ff_off + 76);

        auto g_yyz_yyz = cbuffer.data(ff_off + 77);

        auto g_yyz_yzz = cbuffer.data(ff_off + 78);

        auto g_yzz_xxx = cbuffer.data(ff_off + 80);

        auto g_yzz_xxy = cbuffer.data(ff_off + 81);

        auto g_yzz_xxz = cbuffer.data(ff_off + 82);

        auto g_yzz_xyy = cbuffer.data(ff_off + 83);

        auto g_yzz_xyz = cbuffer.data(ff_off + 84);

        auto g_yzz_xzz = cbuffer.data(ff_off + 85);

        auto g_yzz_yyy = cbuffer.data(ff_off + 86);

        auto g_yzz_yyz = cbuffer.data(ff_off + 87);

        auto g_yzz_yzz = cbuffer.data(ff_off + 88);

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

        /// set up bra offset for contr_buffer_xgd

        const auto gd_off = idx_xgd + i * 90;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_xxxx_xx = cbuffer.data(gd_off + 0);

        auto g_xxxx_xy = cbuffer.data(gd_off + 1);

        auto g_xxxx_xz = cbuffer.data(gd_off + 2);

        auto g_xxxx_yy = cbuffer.data(gd_off + 3);

        auto g_xxxx_yz = cbuffer.data(gd_off + 4);

        auto g_xxxx_zz = cbuffer.data(gd_off + 5);

        #pragma omp simd aligned(cd_x, g_xxx_xx, g_xxx_xxx, g_xxx_xxy, g_xxx_xxz, g_xxx_xy, g_xxx_xyy, g_xxx_xyz, g_xxx_xz, g_xxx_xzz, g_xxx_yy, g_xxx_yz, g_xxx_zz, g_xxxx_xx, g_xxxx_xy, g_xxxx_xz, g_xxxx_yy, g_xxxx_yz, g_xxxx_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxx_xx[k] = -g_xxx_xx[k] * cd_x[k] + g_xxx_xxx[k];

            g_xxxx_xy[k] = -g_xxx_xy[k] * cd_x[k] + g_xxx_xxy[k];

            g_xxxx_xz[k] = -g_xxx_xz[k] * cd_x[k] + g_xxx_xxz[k];

            g_xxxx_yy[k] = -g_xxx_yy[k] * cd_x[k] + g_xxx_xyy[k];

            g_xxxx_yz[k] = -g_xxx_yz[k] * cd_x[k] + g_xxx_xyz[k];

            g_xxxx_zz[k] = -g_xxx_zz[k] * cd_x[k] + g_xxx_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_xxxy_xx = cbuffer.data(gd_off + 6);

        auto g_xxxy_xy = cbuffer.data(gd_off + 7);

        auto g_xxxy_xz = cbuffer.data(gd_off + 8);

        auto g_xxxy_yy = cbuffer.data(gd_off + 9);

        auto g_xxxy_yz = cbuffer.data(gd_off + 10);

        auto g_xxxy_zz = cbuffer.data(gd_off + 11);

        #pragma omp simd aligned(cd_x, g_xxxy_xx, g_xxxy_xy, g_xxxy_xz, g_xxxy_yy, g_xxxy_yz, g_xxxy_zz, g_xxy_xx, g_xxy_xxx, g_xxy_xxy, g_xxy_xxz, g_xxy_xy, g_xxy_xyy, g_xxy_xyz, g_xxy_xz, g_xxy_xzz, g_xxy_yy, g_xxy_yz, g_xxy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxy_xx[k] = -g_xxy_xx[k] * cd_x[k] + g_xxy_xxx[k];

            g_xxxy_xy[k] = -g_xxy_xy[k] * cd_x[k] + g_xxy_xxy[k];

            g_xxxy_xz[k] = -g_xxy_xz[k] * cd_x[k] + g_xxy_xxz[k];

            g_xxxy_yy[k] = -g_xxy_yy[k] * cd_x[k] + g_xxy_xyy[k];

            g_xxxy_yz[k] = -g_xxy_yz[k] * cd_x[k] + g_xxy_xyz[k];

            g_xxxy_zz[k] = -g_xxy_zz[k] * cd_x[k] + g_xxy_xzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_xxxz_xx = cbuffer.data(gd_off + 12);

        auto g_xxxz_xy = cbuffer.data(gd_off + 13);

        auto g_xxxz_xz = cbuffer.data(gd_off + 14);

        auto g_xxxz_yy = cbuffer.data(gd_off + 15);

        auto g_xxxz_yz = cbuffer.data(gd_off + 16);

        auto g_xxxz_zz = cbuffer.data(gd_off + 17);

        #pragma omp simd aligned(cd_x, g_xxxz_xx, g_xxxz_xy, g_xxxz_xz, g_xxxz_yy, g_xxxz_yz, g_xxxz_zz, g_xxz_xx, g_xxz_xxx, g_xxz_xxy, g_xxz_xxz, g_xxz_xy, g_xxz_xyy, g_xxz_xyz, g_xxz_xz, g_xxz_xzz, g_xxz_yy, g_xxz_yz, g_xxz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxxz_xx[k] = -g_xxz_xx[k] * cd_x[k] + g_xxz_xxx[k];

            g_xxxz_xy[k] = -g_xxz_xy[k] * cd_x[k] + g_xxz_xxy[k];

            g_xxxz_xz[k] = -g_xxz_xz[k] * cd_x[k] + g_xxz_xxz[k];

            g_xxxz_yy[k] = -g_xxz_yy[k] * cd_x[k] + g_xxz_xyy[k];

            g_xxxz_yz[k] = -g_xxz_yz[k] * cd_x[k] + g_xxz_xyz[k];

            g_xxxz_zz[k] = -g_xxz_zz[k] * cd_x[k] + g_xxz_xzz[k];
        }

        /// Set up 18-24 components of targeted buffer : cbuffer.data(

        auto g_xxyy_xx = cbuffer.data(gd_off + 18);

        auto g_xxyy_xy = cbuffer.data(gd_off + 19);

        auto g_xxyy_xz = cbuffer.data(gd_off + 20);

        auto g_xxyy_yy = cbuffer.data(gd_off + 21);

        auto g_xxyy_yz = cbuffer.data(gd_off + 22);

        auto g_xxyy_zz = cbuffer.data(gd_off + 23);

        #pragma omp simd aligned(cd_x, g_xxyy_xx, g_xxyy_xy, g_xxyy_xz, g_xxyy_yy, g_xxyy_yz, g_xxyy_zz, g_xyy_xx, g_xyy_xxx, g_xyy_xxy, g_xyy_xxz, g_xyy_xy, g_xyy_xyy, g_xyy_xyz, g_xyy_xz, g_xyy_xzz, g_xyy_yy, g_xyy_yz, g_xyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyy_xx[k] = -g_xyy_xx[k] * cd_x[k] + g_xyy_xxx[k];

            g_xxyy_xy[k] = -g_xyy_xy[k] * cd_x[k] + g_xyy_xxy[k];

            g_xxyy_xz[k] = -g_xyy_xz[k] * cd_x[k] + g_xyy_xxz[k];

            g_xxyy_yy[k] = -g_xyy_yy[k] * cd_x[k] + g_xyy_xyy[k];

            g_xxyy_yz[k] = -g_xyy_yz[k] * cd_x[k] + g_xyy_xyz[k];

            g_xxyy_zz[k] = -g_xyy_zz[k] * cd_x[k] + g_xyy_xzz[k];
        }

        /// Set up 24-30 components of targeted buffer : cbuffer.data(

        auto g_xxyz_xx = cbuffer.data(gd_off + 24);

        auto g_xxyz_xy = cbuffer.data(gd_off + 25);

        auto g_xxyz_xz = cbuffer.data(gd_off + 26);

        auto g_xxyz_yy = cbuffer.data(gd_off + 27);

        auto g_xxyz_yz = cbuffer.data(gd_off + 28);

        auto g_xxyz_zz = cbuffer.data(gd_off + 29);

        #pragma omp simd aligned(cd_x, g_xxyz_xx, g_xxyz_xy, g_xxyz_xz, g_xxyz_yy, g_xxyz_yz, g_xxyz_zz, g_xyz_xx, g_xyz_xxx, g_xyz_xxy, g_xyz_xxz, g_xyz_xy, g_xyz_xyy, g_xyz_xyz, g_xyz_xz, g_xyz_xzz, g_xyz_yy, g_xyz_yz, g_xyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxyz_xx[k] = -g_xyz_xx[k] * cd_x[k] + g_xyz_xxx[k];

            g_xxyz_xy[k] = -g_xyz_xy[k] * cd_x[k] + g_xyz_xxy[k];

            g_xxyz_xz[k] = -g_xyz_xz[k] * cd_x[k] + g_xyz_xxz[k];

            g_xxyz_yy[k] = -g_xyz_yy[k] * cd_x[k] + g_xyz_xyy[k];

            g_xxyz_yz[k] = -g_xyz_yz[k] * cd_x[k] + g_xyz_xyz[k];

            g_xxyz_zz[k] = -g_xyz_zz[k] * cd_x[k] + g_xyz_xzz[k];
        }

        /// Set up 30-36 components of targeted buffer : cbuffer.data(

        auto g_xxzz_xx = cbuffer.data(gd_off + 30);

        auto g_xxzz_xy = cbuffer.data(gd_off + 31);

        auto g_xxzz_xz = cbuffer.data(gd_off + 32);

        auto g_xxzz_yy = cbuffer.data(gd_off + 33);

        auto g_xxzz_yz = cbuffer.data(gd_off + 34);

        auto g_xxzz_zz = cbuffer.data(gd_off + 35);

        #pragma omp simd aligned(cd_x, g_xxzz_xx, g_xxzz_xy, g_xxzz_xz, g_xxzz_yy, g_xxzz_yz, g_xxzz_zz, g_xzz_xx, g_xzz_xxx, g_xzz_xxy, g_xzz_xxz, g_xzz_xy, g_xzz_xyy, g_xzz_xyz, g_xzz_xz, g_xzz_xzz, g_xzz_yy, g_xzz_yz, g_xzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxzz_xx[k] = -g_xzz_xx[k] * cd_x[k] + g_xzz_xxx[k];

            g_xxzz_xy[k] = -g_xzz_xy[k] * cd_x[k] + g_xzz_xxy[k];

            g_xxzz_xz[k] = -g_xzz_xz[k] * cd_x[k] + g_xzz_xxz[k];

            g_xxzz_yy[k] = -g_xzz_yy[k] * cd_x[k] + g_xzz_xyy[k];

            g_xxzz_yz[k] = -g_xzz_yz[k] * cd_x[k] + g_xzz_xyz[k];

            g_xxzz_zz[k] = -g_xzz_zz[k] * cd_x[k] + g_xzz_xzz[k];
        }

        /// Set up 36-42 components of targeted buffer : cbuffer.data(

        auto g_xyyy_xx = cbuffer.data(gd_off + 36);

        auto g_xyyy_xy = cbuffer.data(gd_off + 37);

        auto g_xyyy_xz = cbuffer.data(gd_off + 38);

        auto g_xyyy_yy = cbuffer.data(gd_off + 39);

        auto g_xyyy_yz = cbuffer.data(gd_off + 40);

        auto g_xyyy_zz = cbuffer.data(gd_off + 41);

        #pragma omp simd aligned(cd_x, g_xyyy_xx, g_xyyy_xy, g_xyyy_xz, g_xyyy_yy, g_xyyy_yz, g_xyyy_zz, g_yyy_xx, g_yyy_xxx, g_yyy_xxy, g_yyy_xxz, g_yyy_xy, g_yyy_xyy, g_yyy_xyz, g_yyy_xz, g_yyy_xzz, g_yyy_yy, g_yyy_yz, g_yyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyy_xx[k] = -g_yyy_xx[k] * cd_x[k] + g_yyy_xxx[k];

            g_xyyy_xy[k] = -g_yyy_xy[k] * cd_x[k] + g_yyy_xxy[k];

            g_xyyy_xz[k] = -g_yyy_xz[k] * cd_x[k] + g_yyy_xxz[k];

            g_xyyy_yy[k] = -g_yyy_yy[k] * cd_x[k] + g_yyy_xyy[k];

            g_xyyy_yz[k] = -g_yyy_yz[k] * cd_x[k] + g_yyy_xyz[k];

            g_xyyy_zz[k] = -g_yyy_zz[k] * cd_x[k] + g_yyy_xzz[k];
        }

        /// Set up 42-48 components of targeted buffer : cbuffer.data(

        auto g_xyyz_xx = cbuffer.data(gd_off + 42);

        auto g_xyyz_xy = cbuffer.data(gd_off + 43);

        auto g_xyyz_xz = cbuffer.data(gd_off + 44);

        auto g_xyyz_yy = cbuffer.data(gd_off + 45);

        auto g_xyyz_yz = cbuffer.data(gd_off + 46);

        auto g_xyyz_zz = cbuffer.data(gd_off + 47);

        #pragma omp simd aligned(cd_x, g_xyyz_xx, g_xyyz_xy, g_xyyz_xz, g_xyyz_yy, g_xyyz_yz, g_xyyz_zz, g_yyz_xx, g_yyz_xxx, g_yyz_xxy, g_yyz_xxz, g_yyz_xy, g_yyz_xyy, g_yyz_xyz, g_yyz_xz, g_yyz_xzz, g_yyz_yy, g_yyz_yz, g_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyyz_xx[k] = -g_yyz_xx[k] * cd_x[k] + g_yyz_xxx[k];

            g_xyyz_xy[k] = -g_yyz_xy[k] * cd_x[k] + g_yyz_xxy[k];

            g_xyyz_xz[k] = -g_yyz_xz[k] * cd_x[k] + g_yyz_xxz[k];

            g_xyyz_yy[k] = -g_yyz_yy[k] * cd_x[k] + g_yyz_xyy[k];

            g_xyyz_yz[k] = -g_yyz_yz[k] * cd_x[k] + g_yyz_xyz[k];

            g_xyyz_zz[k] = -g_yyz_zz[k] * cd_x[k] + g_yyz_xzz[k];
        }

        /// Set up 48-54 components of targeted buffer : cbuffer.data(

        auto g_xyzz_xx = cbuffer.data(gd_off + 48);

        auto g_xyzz_xy = cbuffer.data(gd_off + 49);

        auto g_xyzz_xz = cbuffer.data(gd_off + 50);

        auto g_xyzz_yy = cbuffer.data(gd_off + 51);

        auto g_xyzz_yz = cbuffer.data(gd_off + 52);

        auto g_xyzz_zz = cbuffer.data(gd_off + 53);

        #pragma omp simd aligned(cd_x, g_xyzz_xx, g_xyzz_xy, g_xyzz_xz, g_xyzz_yy, g_xyzz_yz, g_xyzz_zz, g_yzz_xx, g_yzz_xxx, g_yzz_xxy, g_yzz_xxz, g_yzz_xy, g_yzz_xyy, g_yzz_xyz, g_yzz_xz, g_yzz_xzz, g_yzz_yy, g_yzz_yz, g_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyzz_xx[k] = -g_yzz_xx[k] * cd_x[k] + g_yzz_xxx[k];

            g_xyzz_xy[k] = -g_yzz_xy[k] * cd_x[k] + g_yzz_xxy[k];

            g_xyzz_xz[k] = -g_yzz_xz[k] * cd_x[k] + g_yzz_xxz[k];

            g_xyzz_yy[k] = -g_yzz_yy[k] * cd_x[k] + g_yzz_xyy[k];

            g_xyzz_yz[k] = -g_yzz_yz[k] * cd_x[k] + g_yzz_xyz[k];

            g_xyzz_zz[k] = -g_yzz_zz[k] * cd_x[k] + g_yzz_xzz[k];
        }

        /// Set up 54-60 components of targeted buffer : cbuffer.data(

        auto g_xzzz_xx = cbuffer.data(gd_off + 54);

        auto g_xzzz_xy = cbuffer.data(gd_off + 55);

        auto g_xzzz_xz = cbuffer.data(gd_off + 56);

        auto g_xzzz_yy = cbuffer.data(gd_off + 57);

        auto g_xzzz_yz = cbuffer.data(gd_off + 58);

        auto g_xzzz_zz = cbuffer.data(gd_off + 59);

        #pragma omp simd aligned(cd_x, g_xzzz_xx, g_xzzz_xy, g_xzzz_xz, g_xzzz_yy, g_xzzz_yz, g_xzzz_zz, g_zzz_xx, g_zzz_xxx, g_zzz_xxy, g_zzz_xxz, g_zzz_xy, g_zzz_xyy, g_zzz_xyz, g_zzz_xz, g_zzz_xzz, g_zzz_yy, g_zzz_yz, g_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzzz_xx[k] = -g_zzz_xx[k] * cd_x[k] + g_zzz_xxx[k];

            g_xzzz_xy[k] = -g_zzz_xy[k] * cd_x[k] + g_zzz_xxy[k];

            g_xzzz_xz[k] = -g_zzz_xz[k] * cd_x[k] + g_zzz_xxz[k];

            g_xzzz_yy[k] = -g_zzz_yy[k] * cd_x[k] + g_zzz_xyy[k];

            g_xzzz_yz[k] = -g_zzz_yz[k] * cd_x[k] + g_zzz_xyz[k];

            g_xzzz_zz[k] = -g_zzz_zz[k] * cd_x[k] + g_zzz_xzz[k];
        }

        /// Set up 60-66 components of targeted buffer : cbuffer.data(

        auto g_yyyy_xx = cbuffer.data(gd_off + 60);

        auto g_yyyy_xy = cbuffer.data(gd_off + 61);

        auto g_yyyy_xz = cbuffer.data(gd_off + 62);

        auto g_yyyy_yy = cbuffer.data(gd_off + 63);

        auto g_yyyy_yz = cbuffer.data(gd_off + 64);

        auto g_yyyy_zz = cbuffer.data(gd_off + 65);

        #pragma omp simd aligned(cd_y, g_yyy_xx, g_yyy_xxy, g_yyy_xy, g_yyy_xyy, g_yyy_xyz, g_yyy_xz, g_yyy_yy, g_yyy_yyy, g_yyy_yyz, g_yyy_yz, g_yyy_yzz, g_yyy_zz, g_yyyy_xx, g_yyyy_xy, g_yyyy_xz, g_yyyy_yy, g_yyyy_yz, g_yyyy_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyy_xx[k] = -g_yyy_xx[k] * cd_y[k] + g_yyy_xxy[k];

            g_yyyy_xy[k] = -g_yyy_xy[k] * cd_y[k] + g_yyy_xyy[k];

            g_yyyy_xz[k] = -g_yyy_xz[k] * cd_y[k] + g_yyy_xyz[k];

            g_yyyy_yy[k] = -g_yyy_yy[k] * cd_y[k] + g_yyy_yyy[k];

            g_yyyy_yz[k] = -g_yyy_yz[k] * cd_y[k] + g_yyy_yyz[k];

            g_yyyy_zz[k] = -g_yyy_zz[k] * cd_y[k] + g_yyy_yzz[k];
        }

        /// Set up 66-72 components of targeted buffer : cbuffer.data(

        auto g_yyyz_xx = cbuffer.data(gd_off + 66);

        auto g_yyyz_xy = cbuffer.data(gd_off + 67);

        auto g_yyyz_xz = cbuffer.data(gd_off + 68);

        auto g_yyyz_yy = cbuffer.data(gd_off + 69);

        auto g_yyyz_yz = cbuffer.data(gd_off + 70);

        auto g_yyyz_zz = cbuffer.data(gd_off + 71);

        #pragma omp simd aligned(cd_y, g_yyyz_xx, g_yyyz_xy, g_yyyz_xz, g_yyyz_yy, g_yyyz_yz, g_yyyz_zz, g_yyz_xx, g_yyz_xxy, g_yyz_xy, g_yyz_xyy, g_yyz_xyz, g_yyz_xz, g_yyz_yy, g_yyz_yyy, g_yyz_yyz, g_yyz_yz, g_yyz_yzz, g_yyz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyyz_xx[k] = -g_yyz_xx[k] * cd_y[k] + g_yyz_xxy[k];

            g_yyyz_xy[k] = -g_yyz_xy[k] * cd_y[k] + g_yyz_xyy[k];

            g_yyyz_xz[k] = -g_yyz_xz[k] * cd_y[k] + g_yyz_xyz[k];

            g_yyyz_yy[k] = -g_yyz_yy[k] * cd_y[k] + g_yyz_yyy[k];

            g_yyyz_yz[k] = -g_yyz_yz[k] * cd_y[k] + g_yyz_yyz[k];

            g_yyyz_zz[k] = -g_yyz_zz[k] * cd_y[k] + g_yyz_yzz[k];
        }

        /// Set up 72-78 components of targeted buffer : cbuffer.data(

        auto g_yyzz_xx = cbuffer.data(gd_off + 72);

        auto g_yyzz_xy = cbuffer.data(gd_off + 73);

        auto g_yyzz_xz = cbuffer.data(gd_off + 74);

        auto g_yyzz_yy = cbuffer.data(gd_off + 75);

        auto g_yyzz_yz = cbuffer.data(gd_off + 76);

        auto g_yyzz_zz = cbuffer.data(gd_off + 77);

        #pragma omp simd aligned(cd_y, g_yyzz_xx, g_yyzz_xy, g_yyzz_xz, g_yyzz_yy, g_yyzz_yz, g_yyzz_zz, g_yzz_xx, g_yzz_xxy, g_yzz_xy, g_yzz_xyy, g_yzz_xyz, g_yzz_xz, g_yzz_yy, g_yzz_yyy, g_yzz_yyz, g_yzz_yz, g_yzz_yzz, g_yzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyzz_xx[k] = -g_yzz_xx[k] * cd_y[k] + g_yzz_xxy[k];

            g_yyzz_xy[k] = -g_yzz_xy[k] * cd_y[k] + g_yzz_xyy[k];

            g_yyzz_xz[k] = -g_yzz_xz[k] * cd_y[k] + g_yzz_xyz[k];

            g_yyzz_yy[k] = -g_yzz_yy[k] * cd_y[k] + g_yzz_yyy[k];

            g_yyzz_yz[k] = -g_yzz_yz[k] * cd_y[k] + g_yzz_yyz[k];

            g_yyzz_zz[k] = -g_yzz_zz[k] * cd_y[k] + g_yzz_yzz[k];
        }

        /// Set up 78-84 components of targeted buffer : cbuffer.data(

        auto g_yzzz_xx = cbuffer.data(gd_off + 78);

        auto g_yzzz_xy = cbuffer.data(gd_off + 79);

        auto g_yzzz_xz = cbuffer.data(gd_off + 80);

        auto g_yzzz_yy = cbuffer.data(gd_off + 81);

        auto g_yzzz_yz = cbuffer.data(gd_off + 82);

        auto g_yzzz_zz = cbuffer.data(gd_off + 83);

        #pragma omp simd aligned(cd_y, g_yzzz_xx, g_yzzz_xy, g_yzzz_xz, g_yzzz_yy, g_yzzz_yz, g_yzzz_zz, g_zzz_xx, g_zzz_xxy, g_zzz_xy, g_zzz_xyy, g_zzz_xyz, g_zzz_xz, g_zzz_yy, g_zzz_yyy, g_zzz_yyz, g_zzz_yz, g_zzz_yzz, g_zzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzzz_xx[k] = -g_zzz_xx[k] * cd_y[k] + g_zzz_xxy[k];

            g_yzzz_xy[k] = -g_zzz_xy[k] * cd_y[k] + g_zzz_xyy[k];

            g_yzzz_xz[k] = -g_zzz_xz[k] * cd_y[k] + g_zzz_xyz[k];

            g_yzzz_yy[k] = -g_zzz_yy[k] * cd_y[k] + g_zzz_yyy[k];

            g_yzzz_yz[k] = -g_zzz_yz[k] * cd_y[k] + g_zzz_yyz[k];

            g_yzzz_zz[k] = -g_zzz_zz[k] * cd_y[k] + g_zzz_yzz[k];
        }

        /// Set up 84-90 components of targeted buffer : cbuffer.data(

        auto g_zzzz_xx = cbuffer.data(gd_off + 84);

        auto g_zzzz_xy = cbuffer.data(gd_off + 85);

        auto g_zzzz_xz = cbuffer.data(gd_off + 86);

        auto g_zzzz_yy = cbuffer.data(gd_off + 87);

        auto g_zzzz_yz = cbuffer.data(gd_off + 88);

        auto g_zzzz_zz = cbuffer.data(gd_off + 89);

        #pragma omp simd aligned(cd_z, g_zzz_xx, g_zzz_xxz, g_zzz_xy, g_zzz_xyz, g_zzz_xz, g_zzz_xzz, g_zzz_yy, g_zzz_yyz, g_zzz_yz, g_zzz_yzz, g_zzz_zz, g_zzz_zzz, g_zzzz_xx, g_zzzz_xy, g_zzzz_xz, g_zzzz_yy, g_zzzz_yz, g_zzzz_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzzz_xx[k] = -g_zzz_xx[k] * cd_z[k] + g_zzz_xxz[k];

            g_zzzz_xy[k] = -g_zzz_xy[k] * cd_z[k] + g_zzz_xyz[k];

            g_zzzz_xz[k] = -g_zzz_xz[k] * cd_z[k] + g_zzz_xzz[k];

            g_zzzz_yy[k] = -g_zzz_yy[k] * cd_z[k] + g_zzz_yyz[k];

            g_zzzz_yz[k] = -g_zzz_yz[k] * cd_z[k] + g_zzz_yzz[k];

            g_zzzz_zz[k] = -g_zzz_zz[k] * cd_z[k] + g_zzz_zzz[k];
        }
    }
}

} // t3ceri namespace

