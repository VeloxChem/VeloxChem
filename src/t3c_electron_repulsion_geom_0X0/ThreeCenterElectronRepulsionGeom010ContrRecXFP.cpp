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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXFP.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xfp(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xfp,
                                        const size_t idx_xdp,
                                        const size_t idx_geom_10_xdp,
                                        const size_t idx_geom_10_xdd,
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

        auto g_yy_x = cbuffer.data(dp_off + 9);

        auto g_yy_y = cbuffer.data(dp_off + 10);

        auto g_yy_z = cbuffer.data(dp_off + 11);

        auto g_zz_x = cbuffer.data(dp_off + 15);

        auto g_zz_y = cbuffer.data(dp_off + 16);

        auto g_zz_z = cbuffer.data(dp_off + 17);

        /// Set up components of auxilary buffer : SDP

        const auto dp_geom_10_off = idx_geom_10_xdp + i * 18;

        auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xy_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xy_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xz_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 6);

        auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 8);

        auto g_x_0_yy_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 9);

        auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 10);

        auto g_x_0_yy_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 11);

        auto g_x_0_yz_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 12);

        auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 13);

        auto g_x_0_yz_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 14);

        auto g_x_0_zz_x = cbuffer.data(dp_geom_10_off + 0 * acomps + 15);

        auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 0 * acomps + 16);

        auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 0 * acomps + 17);

        auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 0);

        auto g_y_0_xx_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 1);

        auto g_y_0_xx_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 2);

        auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 3);

        auto g_y_0_xy_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 4);

        auto g_y_0_xy_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 5);

        auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 6);

        auto g_y_0_xz_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 7);

        auto g_y_0_xz_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 8);

        auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 9);

        auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 10);

        auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 11);

        auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 12);

        auto g_y_0_yz_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 13);

        auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 14);

        auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 18 * acomps + 15);

        auto g_y_0_zz_y = cbuffer.data(dp_geom_10_off + 18 * acomps + 16);

        auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 18 * acomps + 17);

        auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 0);

        auto g_z_0_xx_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 1);

        auto g_z_0_xx_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 2);

        auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 3);

        auto g_z_0_xy_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 4);

        auto g_z_0_xy_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 5);

        auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 6);

        auto g_z_0_xz_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 7);

        auto g_z_0_xz_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 8);

        auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 9);

        auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 10);

        auto g_z_0_yy_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 11);

        auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 12);

        auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 13);

        auto g_z_0_yz_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 14);

        auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 36 * acomps + 15);

        auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 36 * acomps + 16);

        auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 36 * acomps + 17);

        /// Set up components of auxilary buffer : SDD

        const auto dd_geom_10_off = idx_geom_10_xdd + i * 36;

        auto g_x_0_xx_xx = cbuffer.data(dd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_xx_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_xx_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_xx_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_xx_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_xx_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 5);

        auto g_x_0_xy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 7);

        auto g_x_0_xy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 9);

        auto g_x_0_xy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 10);

        auto g_x_0_xz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 13);

        auto g_x_0_xz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 14);

        auto g_x_0_xz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 15);

        auto g_x_0_xz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 16);

        auto g_x_0_xz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 17);

        auto g_x_0_yy_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 19);

        auto g_x_0_yy_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 21);

        auto g_x_0_yy_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 22);

        auto g_x_0_yz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 25);

        auto g_x_0_yz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 27);

        auto g_x_0_yz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 28);

        auto g_x_0_zz_xy = cbuffer.data(dd_geom_10_off + 0 * acomps + 31);

        auto g_x_0_zz_xz = cbuffer.data(dd_geom_10_off + 0 * acomps + 32);

        auto g_x_0_zz_yy = cbuffer.data(dd_geom_10_off + 0 * acomps + 33);

        auto g_x_0_zz_yz = cbuffer.data(dd_geom_10_off + 0 * acomps + 34);

        auto g_x_0_zz_zz = cbuffer.data(dd_geom_10_off + 0 * acomps + 35);

        auto g_y_0_xx_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 0);

        auto g_y_0_xx_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 1);

        auto g_y_0_xx_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 2);

        auto g_y_0_xy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 6);

        auto g_y_0_xy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 7);

        auto g_y_0_xy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 8);

        auto g_y_0_xz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 12);

        auto g_y_0_xz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 13);

        auto g_y_0_xz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 14);

        auto g_y_0_yy_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 18);

        auto g_y_0_yy_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 19);

        auto g_y_0_yy_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 20);

        auto g_y_0_yy_yy = cbuffer.data(dd_geom_10_off + 36 * acomps + 21);

        auto g_y_0_yy_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 22);

        auto g_y_0_yy_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 23);

        auto g_y_0_yz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 24);

        auto g_y_0_yz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 25);

        auto g_y_0_yz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 26);

        auto g_y_0_yz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 28);

        auto g_y_0_yz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 29);

        auto g_y_0_zz_xx = cbuffer.data(dd_geom_10_off + 36 * acomps + 30);

        auto g_y_0_zz_xy = cbuffer.data(dd_geom_10_off + 36 * acomps + 31);

        auto g_y_0_zz_xz = cbuffer.data(dd_geom_10_off + 36 * acomps + 32);

        auto g_y_0_zz_yz = cbuffer.data(dd_geom_10_off + 36 * acomps + 34);

        auto g_y_0_zz_zz = cbuffer.data(dd_geom_10_off + 36 * acomps + 35);

        auto g_z_0_xx_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 0);

        auto g_z_0_xx_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 1);

        auto g_z_0_xx_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 2);

        auto g_z_0_xy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 6);

        auto g_z_0_xy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 7);

        auto g_z_0_xy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 8);

        auto g_z_0_xz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 12);

        auto g_z_0_xz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 13);

        auto g_z_0_xz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 14);

        auto g_z_0_yy_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 18);

        auto g_z_0_yy_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 19);

        auto g_z_0_yy_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 20);

        auto g_z_0_yy_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 21);

        auto g_z_0_yy_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 22);

        auto g_z_0_yz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 24);

        auto g_z_0_yz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 25);

        auto g_z_0_yz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 26);

        auto g_z_0_yz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 27);

        auto g_z_0_yz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 28);

        auto g_z_0_zz_xx = cbuffer.data(dd_geom_10_off + 72 * acomps + 30);

        auto g_z_0_zz_xy = cbuffer.data(dd_geom_10_off + 72 * acomps + 31);

        auto g_z_0_zz_xz = cbuffer.data(dd_geom_10_off + 72 * acomps + 32);

        auto g_z_0_zz_yy = cbuffer.data(dd_geom_10_off + 72 * acomps + 33);

        auto g_z_0_zz_yz = cbuffer.data(dd_geom_10_off + 72 * acomps + 34);

        auto g_z_0_zz_zz = cbuffer.data(dd_geom_10_off + 72 * acomps + 35);

        /// set up bra offset for contr_buffer_xxfp

        const auto fp_geom_10_off = idx_geom_10_xfp + i * 30;

        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_x_0_xx_x, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_y, g_x_0_xx_z, g_x_0_xxx_x, g_x_0_xxx_y, g_x_0_xxx_z, g_xx_x, g_xx_y, g_xx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxx_x[k] = -g_xx_x[k] - g_x_0_xx_x[k] * cd_x[k] + g_x_0_xx_xx[k];

            g_x_0_xxx_y[k] = -g_xx_y[k] - g_x_0_xx_y[k] * cd_x[k] + g_x_0_xx_xy[k];

            g_x_0_xxx_z[k] = -g_xx_z[k] - g_x_0_xx_z[k] * cd_x[k] + g_x_0_xx_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxy_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_xxy_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_y, g_x_0_xx_x, g_x_0_xx_xy, g_x_0_xx_y, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_z, g_x_0_xxy_x, g_x_0_xxy_y, g_x_0_xxy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxy_x[k] = -g_x_0_xx_x[k] * cd_y[k] + g_x_0_xx_xy[k];

            g_x_0_xxy_y[k] = -g_x_0_xx_y[k] * cd_y[k] + g_x_0_xx_yy[k];

            g_x_0_xxy_z[k] = -g_x_0_xx_z[k] * cd_y[k] + g_x_0_xx_yz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_x_0_xxz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 8);

        #pragma omp simd aligned(cd_z, g_x_0_xx_x, g_x_0_xx_xz, g_x_0_xx_y, g_x_0_xx_yz, g_x_0_xx_z, g_x_0_xx_zz, g_x_0_xxz_x, g_x_0_xxz_y, g_x_0_xxz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xxz_x[k] = -g_x_0_xx_x[k] * cd_z[k] + g_x_0_xx_xz[k];

            g_x_0_xxz_y[k] = -g_x_0_xx_y[k] * cd_z[k] + g_x_0_xx_yz[k];

            g_x_0_xxz_z[k] = -g_x_0_xx_z[k] * cd_z[k] + g_x_0_xx_zz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_xyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_xy_x, g_x_0_xy_xy, g_x_0_xy_y, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_z, g_x_0_xyy_x, g_x_0_xyy_y, g_x_0_xyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyy_x[k] = -g_x_0_xy_x[k] * cd_y[k] + g_x_0_xy_xy[k];

            g_x_0_xyy_y[k] = -g_x_0_xy_y[k] * cd_y[k] + g_x_0_xy_yy[k];

            g_x_0_xyy_z[k] = -g_x_0_xy_z[k] * cd_y[k] + g_x_0_xy_yz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_x_0_xyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_xyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 14);

        #pragma omp simd aligned(cd_y, g_x_0_xyz_x, g_x_0_xyz_y, g_x_0_xyz_z, g_x_0_xz_x, g_x_0_xz_xy, g_x_0_xz_y, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xyz_x[k] = -g_x_0_xz_x[k] * cd_y[k] + g_x_0_xz_xy[k];

            g_x_0_xyz_y[k] = -g_x_0_xz_y[k] * cd_y[k] + g_x_0_xz_yy[k];

            g_x_0_xyz_z[k] = -g_x_0_xz_z[k] * cd_y[k] + g_x_0_xz_yz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_xzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_xz_x, g_x_0_xz_xz, g_x_0_xz_y, g_x_0_xz_yz, g_x_0_xz_z, g_x_0_xz_zz, g_x_0_xzz_x, g_x_0_xzz_y, g_x_0_xzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_xzz_x[k] = -g_x_0_xz_x[k] * cd_z[k] + g_x_0_xz_xz[k];

            g_x_0_xzz_y[k] = -g_x_0_xz_y[k] * cd_z[k] + g_x_0_xz_yz[k];

            g_x_0_xzz_z[k] = -g_x_0_xz_z[k] * cd_z[k] + g_x_0_xz_zz[k];
        }

        /// Set up 18-21 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyy_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 18);

        auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 19);

        auto g_x_0_yyy_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 20);

        #pragma omp simd aligned(cd_y, g_x_0_yy_x, g_x_0_yy_xy, g_x_0_yy_y, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_z, g_x_0_yyy_x, g_x_0_yyy_y, g_x_0_yyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyy_x[k] = -g_x_0_yy_x[k] * cd_y[k] + g_x_0_yy_xy[k];

            g_x_0_yyy_y[k] = -g_x_0_yy_y[k] * cd_y[k] + g_x_0_yy_yy[k];

            g_x_0_yyy_z[k] = -g_x_0_yy_z[k] * cd_y[k] + g_x_0_yy_yz[k];
        }

        /// Set up 21-24 components of targeted buffer : cbuffer.data(

        auto g_x_0_yyz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 21);

        auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 22);

        auto g_x_0_yyz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_x_0_yyz_x, g_x_0_yyz_y, g_x_0_yyz_z, g_x_0_yz_x, g_x_0_yz_xy, g_x_0_yz_y, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yyz_x[k] = -g_x_0_yz_x[k] * cd_y[k] + g_x_0_yz_xy[k];

            g_x_0_yyz_y[k] = -g_x_0_yz_y[k] * cd_y[k] + g_x_0_yz_yy[k];

            g_x_0_yyz_z[k] = -g_x_0_yz_z[k] * cd_y[k] + g_x_0_yz_yz[k];
        }

        /// Set up 24-27 components of targeted buffer : cbuffer.data(

        auto g_x_0_yzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 24);

        auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 25);

        auto g_x_0_yzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 26);

        #pragma omp simd aligned(cd_y, g_x_0_yzz_x, g_x_0_yzz_y, g_x_0_yzz_z, g_x_0_zz_x, g_x_0_zz_xy, g_x_0_zz_y, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_yzz_x[k] = -g_x_0_zz_x[k] * cd_y[k] + g_x_0_zz_xy[k];

            g_x_0_yzz_y[k] = -g_x_0_zz_y[k] * cd_y[k] + g_x_0_zz_yy[k];

            g_x_0_yzz_z[k] = -g_x_0_zz_z[k] * cd_y[k] + g_x_0_zz_yz[k];
        }

        /// Set up 27-30 components of targeted buffer : cbuffer.data(

        auto g_x_0_zzz_x = cbuffer.data(fp_geom_10_off + 0 * acomps  + 27);

        auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 0 * acomps  + 28);

        auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 0 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_x_0_zz_x, g_x_0_zz_xz, g_x_0_zz_y, g_x_0_zz_yz, g_x_0_zz_z, g_x_0_zz_zz, g_x_0_zzz_x, g_x_0_zzz_y, g_x_0_zzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_zzz_x[k] = -g_x_0_zz_x[k] * cd_z[k] + g_x_0_zz_xz[k];

            g_x_0_zzz_y[k] = -g_x_0_zz_y[k] * cd_z[k] + g_x_0_zz_yz[k];

            g_x_0_zzz_z[k] = -g_x_0_zz_z[k] * cd_z[k] + g_x_0_zz_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 0);

        auto g_y_0_xxx_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 1);

        auto g_y_0_xxx_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_y_0_xx_x, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_y, g_y_0_xx_z, g_y_0_xxx_x, g_y_0_xxx_y, g_y_0_xxx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxx_x[k] = -g_y_0_xx_x[k] * cd_x[k] + g_y_0_xx_xx[k];

            g_y_0_xxx_y[k] = -g_y_0_xx_y[k] * cd_x[k] + g_y_0_xx_xy[k];

            g_y_0_xxx_z[k] = -g_y_0_xx_z[k] * cd_x[k] + g_y_0_xx_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 3);

        auto g_y_0_xxy_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 4);

        auto g_y_0_xxy_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_xxy_x, g_y_0_xxy_y, g_y_0_xxy_z, g_y_0_xy_x, g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_y, g_y_0_xy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxy_x[k] = -g_y_0_xy_x[k] * cd_x[k] + g_y_0_xy_xx[k];

            g_y_0_xxy_y[k] = -g_y_0_xy_y[k] * cd_x[k] + g_y_0_xy_xy[k];

            g_y_0_xxy_z[k] = -g_y_0_xy_z[k] * cd_x[k] + g_y_0_xy_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 6);

        auto g_y_0_xxz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 7);

        auto g_y_0_xxz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_y_0_xxz_x, g_y_0_xxz_y, g_y_0_xxz_z, g_y_0_xz_x, g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_y, g_y_0_xz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xxz_x[k] = -g_y_0_xz_x[k] * cd_x[k] + g_y_0_xz_xx[k];

            g_y_0_xxz_y[k] = -g_y_0_xz_y[k] * cd_x[k] + g_y_0_xz_xy[k];

            g_y_0_xxz_z[k] = -g_y_0_xz_z[k] * cd_x[k] + g_y_0_xz_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 9);

        auto g_y_0_xyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 10);

        auto g_y_0_xyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_y_0_xyy_x, g_y_0_xyy_y, g_y_0_xyy_z, g_y_0_yy_x, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_y, g_y_0_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyy_x[k] = -g_y_0_yy_x[k] * cd_x[k] + g_y_0_yy_xx[k];

            g_y_0_xyy_y[k] = -g_y_0_yy_y[k] * cd_x[k] + g_y_0_yy_xy[k];

            g_y_0_xyy_z[k] = -g_y_0_yy_z[k] * cd_x[k] + g_y_0_yy_xz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 12);

        auto g_y_0_xyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 13);

        auto g_y_0_xyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_y_0_xyz_x, g_y_0_xyz_y, g_y_0_xyz_z, g_y_0_yz_x, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_y, g_y_0_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xyz_x[k] = -g_y_0_yz_x[k] * cd_x[k] + g_y_0_yz_xx[k];

            g_y_0_xyz_y[k] = -g_y_0_yz_y[k] * cd_x[k] + g_y_0_yz_xy[k];

            g_y_0_xyz_z[k] = -g_y_0_yz_z[k] * cd_x[k] + g_y_0_yz_xz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 15);

        auto g_y_0_xzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 16);

        auto g_y_0_xzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_y_0_xzz_x, g_y_0_xzz_y, g_y_0_xzz_z, g_y_0_zz_x, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_y, g_y_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_xzz_x[k] = -g_y_0_zz_x[k] * cd_x[k] + g_y_0_zz_xx[k];

            g_y_0_xzz_y[k] = -g_y_0_zz_y[k] * cd_x[k] + g_y_0_zz_xy[k];

            g_y_0_xzz_z[k] = -g_y_0_zz_z[k] * cd_x[k] + g_y_0_zz_xz[k];
        }

        /// Set up 18-21 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 18);

        auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 19);

        auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 20);

        #pragma omp simd aligned(cd_y, g_y_0_yy_x, g_y_0_yy_xy, g_y_0_yy_y, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_z, g_y_0_yyy_x, g_y_0_yyy_y, g_y_0_yyy_z, g_yy_x, g_yy_y, g_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyy_x[k] = -g_yy_x[k] - g_y_0_yy_x[k] * cd_y[k] + g_y_0_yy_xy[k];

            g_y_0_yyy_y[k] = -g_yy_y[k] - g_y_0_yy_y[k] * cd_y[k] + g_y_0_yy_yy[k];

            g_y_0_yyy_z[k] = -g_yy_z[k] - g_y_0_yy_z[k] * cd_y[k] + g_y_0_yy_yz[k];
        }

        /// Set up 21-24 components of targeted buffer : cbuffer.data(

        auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 21);

        auto g_y_0_yyz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 22);

        auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 23);

        #pragma omp simd aligned(cd_z, g_y_0_yy_x, g_y_0_yy_xz, g_y_0_yy_y, g_y_0_yy_yz, g_y_0_yy_z, g_y_0_yy_zz, g_y_0_yyz_x, g_y_0_yyz_y, g_y_0_yyz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yyz_x[k] = -g_y_0_yy_x[k] * cd_z[k] + g_y_0_yy_xz[k];

            g_y_0_yyz_y[k] = -g_y_0_yy_y[k] * cd_z[k] + g_y_0_yy_yz[k];

            g_y_0_yyz_z[k] = -g_y_0_yy_z[k] * cd_z[k] + g_y_0_yy_zz[k];
        }

        /// Set up 24-27 components of targeted buffer : cbuffer.data(

        auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 24);

        auto g_y_0_yzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 25);

        auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 26);

        #pragma omp simd aligned(cd_z, g_y_0_yz_x, g_y_0_yz_xz, g_y_0_yz_y, g_y_0_yz_yz, g_y_0_yz_z, g_y_0_yz_zz, g_y_0_yzz_x, g_y_0_yzz_y, g_y_0_yzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_yzz_x[k] = -g_y_0_yz_x[k] * cd_z[k] + g_y_0_yz_xz[k];

            g_y_0_yzz_y[k] = -g_y_0_yz_y[k] * cd_z[k] + g_y_0_yz_yz[k];

            g_y_0_yzz_z[k] = -g_y_0_yz_z[k] * cd_z[k] + g_y_0_yz_zz[k];
        }

        /// Set up 27-30 components of targeted buffer : cbuffer.data(

        auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 30 * acomps  + 27);

        auto g_y_0_zzz_y = cbuffer.data(fp_geom_10_off + 30 * acomps  + 28);

        auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 30 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_y_0_zz_x, g_y_0_zz_xz, g_y_0_zz_y, g_y_0_zz_yz, g_y_0_zz_z, g_y_0_zz_zz, g_y_0_zzz_x, g_y_0_zzz_y, g_y_0_zzz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_zzz_x[k] = -g_y_0_zz_x[k] * cd_z[k] + g_y_0_zz_xz[k];

            g_y_0_zzz_y[k] = -g_y_0_zz_y[k] * cd_z[k] + g_y_0_zz_yz[k];

            g_y_0_zzz_z[k] = -g_y_0_zz_z[k] * cd_z[k] + g_y_0_zz_zz[k];
        }
        /// Set up 0-3 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 0);

        auto g_z_0_xxx_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 1);

        auto g_z_0_xxx_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 2);

        #pragma omp simd aligned(cd_x, g_z_0_xx_x, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_y, g_z_0_xx_z, g_z_0_xxx_x, g_z_0_xxx_y, g_z_0_xxx_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxx_x[k] = -g_z_0_xx_x[k] * cd_x[k] + g_z_0_xx_xx[k];

            g_z_0_xxx_y[k] = -g_z_0_xx_y[k] * cd_x[k] + g_z_0_xx_xy[k];

            g_z_0_xxx_z[k] = -g_z_0_xx_z[k] * cd_x[k] + g_z_0_xx_xz[k];
        }

        /// Set up 3-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 3);

        auto g_z_0_xxy_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 4);

        auto g_z_0_xxy_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_xxy_x, g_z_0_xxy_y, g_z_0_xxy_z, g_z_0_xy_x, g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_y, g_z_0_xy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxy_x[k] = -g_z_0_xy_x[k] * cd_x[k] + g_z_0_xy_xx[k];

            g_z_0_xxy_y[k] = -g_z_0_xy_y[k] * cd_x[k] + g_z_0_xy_xy[k];

            g_z_0_xxy_z[k] = -g_z_0_xy_z[k] * cd_x[k] + g_z_0_xy_xz[k];
        }

        /// Set up 6-9 components of targeted buffer : cbuffer.data(

        auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 6);

        auto g_z_0_xxz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 7);

        auto g_z_0_xxz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 8);

        #pragma omp simd aligned(cd_x, g_z_0_xxz_x, g_z_0_xxz_y, g_z_0_xxz_z, g_z_0_xz_x, g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_y, g_z_0_xz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xxz_x[k] = -g_z_0_xz_x[k] * cd_x[k] + g_z_0_xz_xx[k];

            g_z_0_xxz_y[k] = -g_z_0_xz_y[k] * cd_x[k] + g_z_0_xz_xy[k];

            g_z_0_xxz_z[k] = -g_z_0_xz_z[k] * cd_x[k] + g_z_0_xz_xz[k];
        }

        /// Set up 9-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 9);

        auto g_z_0_xyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 10);

        auto g_z_0_xyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 11);

        #pragma omp simd aligned(cd_x, g_z_0_xyy_x, g_z_0_xyy_y, g_z_0_xyy_z, g_z_0_yy_x, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_y, g_z_0_yy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyy_x[k] = -g_z_0_yy_x[k] * cd_x[k] + g_z_0_yy_xx[k];

            g_z_0_xyy_y[k] = -g_z_0_yy_y[k] * cd_x[k] + g_z_0_yy_xy[k];

            g_z_0_xyy_z[k] = -g_z_0_yy_z[k] * cd_x[k] + g_z_0_yy_xz[k];
        }

        /// Set up 12-15 components of targeted buffer : cbuffer.data(

        auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 12);

        auto g_z_0_xyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 13);

        auto g_z_0_xyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 14);

        #pragma omp simd aligned(cd_x, g_z_0_xyz_x, g_z_0_xyz_y, g_z_0_xyz_z, g_z_0_yz_x, g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_y, g_z_0_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xyz_x[k] = -g_z_0_yz_x[k] * cd_x[k] + g_z_0_yz_xx[k];

            g_z_0_xyz_y[k] = -g_z_0_yz_y[k] * cd_x[k] + g_z_0_yz_xy[k];

            g_z_0_xyz_z[k] = -g_z_0_yz_z[k] * cd_x[k] + g_z_0_yz_xz[k];
        }

        /// Set up 15-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 15);

        auto g_z_0_xzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 16);

        auto g_z_0_xzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 17);

        #pragma omp simd aligned(cd_x, g_z_0_xzz_x, g_z_0_xzz_y, g_z_0_xzz_z, g_z_0_zz_x, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_y, g_z_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_xzz_x[k] = -g_z_0_zz_x[k] * cd_x[k] + g_z_0_zz_xx[k];

            g_z_0_xzz_y[k] = -g_z_0_zz_y[k] * cd_x[k] + g_z_0_zz_xy[k];

            g_z_0_xzz_z[k] = -g_z_0_zz_z[k] * cd_x[k] + g_z_0_zz_xz[k];
        }

        /// Set up 18-21 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 18);

        auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 19);

        auto g_z_0_yyy_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 20);

        #pragma omp simd aligned(cd_y, g_z_0_yy_x, g_z_0_yy_xy, g_z_0_yy_y, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_z, g_z_0_yyy_x, g_z_0_yyy_y, g_z_0_yyy_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyy_x[k] = -g_z_0_yy_x[k] * cd_y[k] + g_z_0_yy_xy[k];

            g_z_0_yyy_y[k] = -g_z_0_yy_y[k] * cd_y[k] + g_z_0_yy_yy[k];

            g_z_0_yyy_z[k] = -g_z_0_yy_z[k] * cd_y[k] + g_z_0_yy_yz[k];
        }

        /// Set up 21-24 components of targeted buffer : cbuffer.data(

        auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 21);

        auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 22);

        auto g_z_0_yyz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 23);

        #pragma omp simd aligned(cd_y, g_z_0_yyz_x, g_z_0_yyz_y, g_z_0_yyz_z, g_z_0_yz_x, g_z_0_yz_xy, g_z_0_yz_y, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yyz_x[k] = -g_z_0_yz_x[k] * cd_y[k] + g_z_0_yz_xy[k];

            g_z_0_yyz_y[k] = -g_z_0_yz_y[k] * cd_y[k] + g_z_0_yz_yy[k];

            g_z_0_yyz_z[k] = -g_z_0_yz_z[k] * cd_y[k] + g_z_0_yz_yz[k];
        }

        /// Set up 24-27 components of targeted buffer : cbuffer.data(

        auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 24);

        auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 25);

        auto g_z_0_yzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 26);

        #pragma omp simd aligned(cd_y, g_z_0_yzz_x, g_z_0_yzz_y, g_z_0_yzz_z, g_z_0_zz_x, g_z_0_zz_xy, g_z_0_zz_y, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_yzz_x[k] = -g_z_0_zz_x[k] * cd_y[k] + g_z_0_zz_xy[k];

            g_z_0_yzz_y[k] = -g_z_0_zz_y[k] * cd_y[k] + g_z_0_zz_yy[k];

            g_z_0_yzz_z[k] = -g_z_0_zz_z[k] * cd_y[k] + g_z_0_zz_yz[k];
        }

        /// Set up 27-30 components of targeted buffer : cbuffer.data(

        auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 60 * acomps  + 27);

        auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 60 * acomps  + 28);

        auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 60 * acomps  + 29);

        #pragma omp simd aligned(cd_z, g_z_0_zz_x, g_z_0_zz_xz, g_z_0_zz_y, g_z_0_zz_yz, g_z_0_zz_z, g_z_0_zz_zz, g_z_0_zzz_x, g_z_0_zzz_y, g_z_0_zzz_z, g_zz_x, g_zz_y, g_zz_z  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_zzz_x[k] = -g_zz_x[k] - g_z_0_zz_x[k] * cd_z[k] + g_z_0_zz_xz[k];

            g_z_0_zzz_y[k] = -g_zz_y[k] - g_z_0_zz_y[k] * cd_z[k] + g_z_0_zz_yz[k];

            g_z_0_zzz_z[k] = -g_zz_z[k] - g_z_0_zz_z[k] * cd_z[k] + g_z_0_zz_zz[k];
        }
    }
}

} // t3ceri namespace

