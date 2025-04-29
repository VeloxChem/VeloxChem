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

#include "ThreeCenterElectronRepulsionContrRecXFP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_hrr_electron_repulsion_xfp(CSimdArray<double>& cbuffer,
                                const size_t idx_xfp,
                                const size_t idx_xdp,
                                const size_t idx_xdd,
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
        /// Set up components of auxilary buffer : SDP

        const auto dp_off = idx_xdp + i * 18;

        auto g_xx_x = cbuffer.data(dp_off + 0);

        auto g_xx_y = cbuffer.data(dp_off + 1);

        auto g_xx_z = cbuffer.data(dp_off + 2);

        auto g_xy_x = cbuffer.data(dp_off + 3);

        auto g_xy_y = cbuffer.data(dp_off + 4);

        auto g_xy_z = cbuffer.data(dp_off + 5);

        auto g_xz_x = cbuffer.data(dp_off + 6);

        auto g_xz_y = cbuffer.data(dp_off + 7);

        auto g_xz_z = cbuffer.data(dp_off + 8);

        auto g_yy_x = cbuffer.data(dp_off + 9);

        auto g_yy_y = cbuffer.data(dp_off + 10);

        auto g_yy_z = cbuffer.data(dp_off + 11);

        auto g_yz_x = cbuffer.data(dp_off + 12);

        auto g_yz_y = cbuffer.data(dp_off + 13);

        auto g_yz_z = cbuffer.data(dp_off + 14);

        auto g_zz_x = cbuffer.data(dp_off + 15);

        auto g_zz_y = cbuffer.data(dp_off + 16);

        auto g_zz_z = cbuffer.data(dp_off + 17);

        /// Set up components of auxilary buffer : SDD

        const auto dd_off = idx_xdd + i * 36;

        auto g_xx_xx = cbuffer.data(dd_off + 0);

        auto g_xx_xy = cbuffer.data(dd_off + 1);

        auto g_xx_xz = cbuffer.data(dd_off + 2);

        auto g_xy_xx = cbuffer.data(dd_off + 6);

        auto g_xy_xy = cbuffer.data(dd_off + 7);

        auto g_xy_xz = cbuffer.data(dd_off + 8);

        auto g_xz_xx = cbuffer.data(dd_off + 12);

        auto g_xz_xy = cbuffer.data(dd_off + 13);

        auto g_xz_xz = cbuffer.data(dd_off + 14);

        auto g_yy_xx = cbuffer.data(dd_off + 18);

        auto g_yy_xy = cbuffer.data(dd_off + 19);

        auto g_yy_xz = cbuffer.data(dd_off + 20);

        auto g_yy_yy = cbuffer.data(dd_off + 21);

        auto g_yy_yz = cbuffer.data(dd_off + 22);

        auto g_yz_xx = cbuffer.data(dd_off + 24);

        auto g_yz_xy = cbuffer.data(dd_off + 25);

        auto g_yz_xz = cbuffer.data(dd_off + 26);

        auto g_yz_yy = cbuffer.data(dd_off + 27);

        auto g_yz_yz = cbuffer.data(dd_off + 28);

        auto g_zz_xx = cbuffer.data(dd_off + 30);

        auto g_zz_xy = cbuffer.data(dd_off + 31);

        auto g_zz_xz = cbuffer.data(dd_off + 32);

        auto g_zz_yy = cbuffer.data(dd_off + 33);

        auto g_zz_yz = cbuffer.data(dd_off + 34);

        auto g_zz_zz = cbuffer.data(dd_off + 35);

        /// set up bra offset for contr_buffer_xfp

        const auto fp_off = idx_xfp + i * 30;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_xxx_x = cbuffer.data(fp_off + 0);

        auto g_xxx_y = cbuffer.data(fp_off + 1);

        auto g_xxx_z = cbuffer.data(fp_off + 2);

        #pragma omp simd aligned(cd_x, g_xx_x, g_xx_xx, g_xx_xy, g_xx_xz, g_xx_y, g_xx_z, g_xxx_x, g_xxx_y, g_xxx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxx_x[k] = -g_xx_x[k] * cd_x[k] + g_xx_xx[k];

            g_xxx_y[k] = -g_xx_y[k] * cd_x[k] + g_xx_xy[k];

            g_xxx_z[k] = -g_xx_z[k] * cd_x[k] + g_xx_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_xxy_x = cbuffer.data(fp_off + 3);

        auto g_xxy_y = cbuffer.data(fp_off + 4);

        auto g_xxy_z = cbuffer.data(fp_off + 5);

        #pragma omp simd aligned(cd_x, g_xxy_x, g_xxy_y, g_xxy_z, g_xy_x, g_xy_xx, g_xy_xy, g_xy_xz, g_xy_y, g_xy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxy_x[k] = -g_xy_x[k] * cd_x[k] + g_xy_xx[k];

            g_xxy_y[k] = -g_xy_y[k] * cd_x[k] + g_xy_xy[k];

            g_xxy_z[k] = -g_xy_z[k] * cd_x[k] + g_xy_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_xxz_x = cbuffer.data(fp_off + 6);

        auto g_xxz_y = cbuffer.data(fp_off + 7);

        auto g_xxz_z = cbuffer.data(fp_off + 8);

        #pragma omp simd aligned(cd_x, g_xxz_x, g_xxz_y, g_xxz_z, g_xz_x, g_xz_xx, g_xz_xy, g_xz_xz, g_xz_y, g_xz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xxz_x[k] = -g_xz_x[k] * cd_x[k] + g_xz_xx[k];

            g_xxz_y[k] = -g_xz_y[k] * cd_x[k] + g_xz_xy[k];

            g_xxz_z[k] = -g_xz_z[k] * cd_x[k] + g_xz_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_xyy_x = cbuffer.data(fp_off + 9);

        auto g_xyy_y = cbuffer.data(fp_off + 10);

        auto g_xyy_z = cbuffer.data(fp_off + 11);

        #pragma omp simd aligned(cd_x, g_xyy_x, g_xyy_y, g_xyy_z, g_yy_x, g_yy_xx, g_yy_xy, g_yy_xz, g_yy_y, g_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyy_x[k] = -g_yy_x[k] * cd_x[k] + g_yy_xx[k];

            g_xyy_y[k] = -g_yy_y[k] * cd_x[k] + g_yy_xy[k];

            g_xyy_z[k] = -g_yy_z[k] * cd_x[k] + g_yy_xz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_xyz_x = cbuffer.data(fp_off + 12);

        auto g_xyz_y = cbuffer.data(fp_off + 13);

        auto g_xyz_z = cbuffer.data(fp_off + 14);

        #pragma omp simd aligned(cd_x, g_xyz_x, g_xyz_y, g_xyz_z, g_yz_x, g_yz_xx, g_yz_xy, g_yz_xz, g_yz_y, g_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xyz_x[k] = -g_yz_x[k] * cd_x[k] + g_yz_xx[k];

            g_xyz_y[k] = -g_yz_y[k] * cd_x[k] + g_yz_xy[k];

            g_xyz_z[k] = -g_yz_z[k] * cd_x[k] + g_yz_xz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_xzz_x = cbuffer.data(fp_off + 15);

        auto g_xzz_y = cbuffer.data(fp_off + 16);

        auto g_xzz_z = cbuffer.data(fp_off + 17);

        #pragma omp simd aligned(cd_x, g_xzz_x, g_xzz_y, g_xzz_z, g_zz_x, g_zz_xx, g_zz_xy, g_zz_xz, g_zz_y, g_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_xzz_x[k] = -g_zz_x[k] * cd_x[k] + g_zz_xx[k];

            g_xzz_y[k] = -g_zz_y[k] * cd_x[k] + g_zz_xy[k];

            g_xzz_z[k] = -g_zz_z[k] * cd_x[k] + g_zz_xz[k];
        }

        /// Set up 18-21 components of targeted buffer : cbuffer.data(

        auto g_yyy_x = cbuffer.data(fp_off + 18);

        auto g_yyy_y = cbuffer.data(fp_off + 19);

        auto g_yyy_z = cbuffer.data(fp_off + 20);

        #pragma omp simd aligned(cd_y, g_yy_x, g_yy_xy, g_yy_y, g_yy_yy, g_yy_yz, g_yy_z, g_yyy_x, g_yyy_y, g_yyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyy_x[k] = -g_yy_x[k] * cd_y[k] + g_yy_xy[k];

            g_yyy_y[k] = -g_yy_y[k] * cd_y[k] + g_yy_yy[k];

            g_yyy_z[k] = -g_yy_z[k] * cd_y[k] + g_yy_yz[k];
        }

        /// Set up 21-24 components of targeted buffer : cbuffer.data(

        auto g_yyz_x = cbuffer.data(fp_off + 21);

        auto g_yyz_y = cbuffer.data(fp_off + 22);

        auto g_yyz_z = cbuffer.data(fp_off + 23);

        #pragma omp simd aligned(cd_y, g_yyz_x, g_yyz_y, g_yyz_z, g_yz_x, g_yz_xy, g_yz_y, g_yz_yy, g_yz_yz, g_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yyz_x[k] = -g_yz_x[k] * cd_y[k] + g_yz_xy[k];

            g_yyz_y[k] = -g_yz_y[k] * cd_y[k] + g_yz_yy[k];

            g_yyz_z[k] = -g_yz_z[k] * cd_y[k] + g_yz_yz[k];
        }

        /// Set up 24-27 components of targeted buffer : cbuffer.data(

        auto g_yzz_x = cbuffer.data(fp_off + 24);

        auto g_yzz_y = cbuffer.data(fp_off + 25);

        auto g_yzz_z = cbuffer.data(fp_off + 26);

        #pragma omp simd aligned(cd_y, g_yzz_x, g_yzz_y, g_yzz_z, g_zz_x, g_zz_xy, g_zz_y, g_zz_yy, g_zz_yz, g_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_yzz_x[k] = -g_zz_x[k] * cd_y[k] + g_zz_xy[k];

            g_yzz_y[k] = -g_zz_y[k] * cd_y[k] + g_zz_yy[k];

            g_yzz_z[k] = -g_zz_z[k] * cd_y[k] + g_zz_yz[k];
        }

        /// Set up 27-30 components of targeted buffer : cbuffer.data(

        auto g_zzz_x = cbuffer.data(fp_off + 27);

        auto g_zzz_y = cbuffer.data(fp_off + 28);

        auto g_zzz_z = cbuffer.data(fp_off + 29);

        #pragma omp simd aligned(cd_z, g_zz_x, g_zz_xz, g_zz_y, g_zz_yz, g_zz_z, g_zz_zz, g_zzz_x, g_zzz_y, g_zzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_zzz_x[k] = -g_zz_x[k] * cd_z[k] + g_zz_xz[k];

            g_zzz_y[k] = -g_zz_y[k] * cd_z[k] + g_zz_yz[k];

            g_zzz_z[k] = -g_zz_z[k] * cd_z[k] + g_zz_zz[k];
        }
    }
}

} // t3ceri namespace

