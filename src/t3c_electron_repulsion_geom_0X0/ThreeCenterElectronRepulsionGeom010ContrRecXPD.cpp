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

#include "ThreeCenterElectronRepulsionGeom010ContrRecXPD.hpp"

#include "TensorComponents.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_ket_geom010_electron_repulsion_xpd(CSimdArray<double>& cbuffer,
                                        const size_t idx_geom_10_xpd,
                                        const size_t idx_xsd,
                                        const size_t idx_geom_10_xsd,
                                        const size_t idx_geom_10_xsf,
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
        /// Set up components of auxilary buffer : SSD

        const auto sd_off = idx_xsd + i * 6;

        auto g_0_xx = cbuffer.data(sd_off + 0);

        auto g_0_xy = cbuffer.data(sd_off + 1);

        auto g_0_xz = cbuffer.data(sd_off + 2);

        auto g_0_yy = cbuffer.data(sd_off + 3);

        auto g_0_yz = cbuffer.data(sd_off + 4);

        auto g_0_zz = cbuffer.data(sd_off + 5);

        /// Set up components of auxilary buffer : SSD

        const auto sd_geom_10_off = idx_geom_10_xsd + i * 6;

        auto g_x_0_0_xx = cbuffer.data(sd_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xy = cbuffer.data(sd_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xz = cbuffer.data(sd_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_yy = cbuffer.data(sd_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_yz = cbuffer.data(sd_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_zz = cbuffer.data(sd_geom_10_off + 0 * acomps + 5);

        auto g_y_0_0_xx = cbuffer.data(sd_geom_10_off + 6 * acomps + 0);

        auto g_y_0_0_xy = cbuffer.data(sd_geom_10_off + 6 * acomps + 1);

        auto g_y_0_0_xz = cbuffer.data(sd_geom_10_off + 6 * acomps + 2);

        auto g_y_0_0_yy = cbuffer.data(sd_geom_10_off + 6 * acomps + 3);

        auto g_y_0_0_yz = cbuffer.data(sd_geom_10_off + 6 * acomps + 4);

        auto g_y_0_0_zz = cbuffer.data(sd_geom_10_off + 6 * acomps + 5);

        auto g_z_0_0_xx = cbuffer.data(sd_geom_10_off + 12 * acomps + 0);

        auto g_z_0_0_xy = cbuffer.data(sd_geom_10_off + 12 * acomps + 1);

        auto g_z_0_0_xz = cbuffer.data(sd_geom_10_off + 12 * acomps + 2);

        auto g_z_0_0_yy = cbuffer.data(sd_geom_10_off + 12 * acomps + 3);

        auto g_z_0_0_yz = cbuffer.data(sd_geom_10_off + 12 * acomps + 4);

        auto g_z_0_0_zz = cbuffer.data(sd_geom_10_off + 12 * acomps + 5);

        /// Set up components of auxilary buffer : SSF

        const auto sf_geom_10_off = idx_geom_10_xsf + i * 10;

        auto g_x_0_0_xxx = cbuffer.data(sf_geom_10_off + 0 * acomps + 0);

        auto g_x_0_0_xxy = cbuffer.data(sf_geom_10_off + 0 * acomps + 1);

        auto g_x_0_0_xxz = cbuffer.data(sf_geom_10_off + 0 * acomps + 2);

        auto g_x_0_0_xyy = cbuffer.data(sf_geom_10_off + 0 * acomps + 3);

        auto g_x_0_0_xyz = cbuffer.data(sf_geom_10_off + 0 * acomps + 4);

        auto g_x_0_0_xzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 5);

        auto g_x_0_0_yyy = cbuffer.data(sf_geom_10_off + 0 * acomps + 6);

        auto g_x_0_0_yyz = cbuffer.data(sf_geom_10_off + 0 * acomps + 7);

        auto g_x_0_0_yzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 8);

        auto g_x_0_0_zzz = cbuffer.data(sf_geom_10_off + 0 * acomps + 9);

        auto g_y_0_0_xxx = cbuffer.data(sf_geom_10_off + 10 * acomps + 0);

        auto g_y_0_0_xxy = cbuffer.data(sf_geom_10_off + 10 * acomps + 1);

        auto g_y_0_0_xxz = cbuffer.data(sf_geom_10_off + 10 * acomps + 2);

        auto g_y_0_0_xyy = cbuffer.data(sf_geom_10_off + 10 * acomps + 3);

        auto g_y_0_0_xyz = cbuffer.data(sf_geom_10_off + 10 * acomps + 4);

        auto g_y_0_0_xzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 5);

        auto g_y_0_0_yyy = cbuffer.data(sf_geom_10_off + 10 * acomps + 6);

        auto g_y_0_0_yyz = cbuffer.data(sf_geom_10_off + 10 * acomps + 7);

        auto g_y_0_0_yzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 8);

        auto g_y_0_0_zzz = cbuffer.data(sf_geom_10_off + 10 * acomps + 9);

        auto g_z_0_0_xxx = cbuffer.data(sf_geom_10_off + 20 * acomps + 0);

        auto g_z_0_0_xxy = cbuffer.data(sf_geom_10_off + 20 * acomps + 1);

        auto g_z_0_0_xxz = cbuffer.data(sf_geom_10_off + 20 * acomps + 2);

        auto g_z_0_0_xyy = cbuffer.data(sf_geom_10_off + 20 * acomps + 3);

        auto g_z_0_0_xyz = cbuffer.data(sf_geom_10_off + 20 * acomps + 4);

        auto g_z_0_0_xzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 5);

        auto g_z_0_0_yyy = cbuffer.data(sf_geom_10_off + 20 * acomps + 6);

        auto g_z_0_0_yyz = cbuffer.data(sf_geom_10_off + 20 * acomps + 7);

        auto g_z_0_0_yzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 8);

        auto g_z_0_0_zzz = cbuffer.data(sf_geom_10_off + 20 * acomps + 9);

        /// set up bra offset for contr_buffer_xxpd

        const auto pd_geom_10_off = idx_geom_10_xpd + i * 18;

        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_x_0_x_xx = cbuffer.data(pd_geom_10_off + 0 * acomps  + 0);

        auto g_x_0_x_xy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 1);

        auto g_x_0_x_xz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 2);

        auto g_x_0_x_yy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 3);

        auto g_x_0_x_yz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 4);

        auto g_x_0_x_zz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz, g_0_zz, g_x_0_0_xx, g_x_0_0_xxx, g_x_0_0_xxy, g_x_0_0_xxz, g_x_0_0_xy, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xz, g_x_0_0_xzz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_x_0_x_xx, g_x_0_x_xy, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_x_xx[k] = -g_0_xx[k] - g_x_0_0_xx[k] * cd_x[k] + g_x_0_0_xxx[k];

            g_x_0_x_xy[k] = -g_0_xy[k] - g_x_0_0_xy[k] * cd_x[k] + g_x_0_0_xxy[k];

            g_x_0_x_xz[k] = -g_0_xz[k] - g_x_0_0_xz[k] * cd_x[k] + g_x_0_0_xxz[k];

            g_x_0_x_yy[k] = -g_0_yy[k] - g_x_0_0_yy[k] * cd_x[k] + g_x_0_0_xyy[k];

            g_x_0_x_yz[k] = -g_0_yz[k] - g_x_0_0_yz[k] * cd_x[k] + g_x_0_0_xyz[k];

            g_x_0_x_zz[k] = -g_0_zz[k] - g_x_0_0_zz[k] * cd_x[k] + g_x_0_0_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_x_0_y_xx = cbuffer.data(pd_geom_10_off + 0 * acomps  + 6);

        auto g_x_0_y_xy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 7);

        auto g_x_0_y_xz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 8);

        auto g_x_0_y_yy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 9);

        auto g_x_0_y_yz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 10);

        auto g_x_0_y_zz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_x_0_0_xx, g_x_0_0_xxy, g_x_0_0_xy, g_x_0_0_xyy, g_x_0_0_xyz, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yyy, g_x_0_0_yyz, g_x_0_0_yz, g_x_0_0_yzz, g_x_0_0_zz, g_x_0_y_xx, g_x_0_y_xy, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yz, g_x_0_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_y_xx[k] = -g_x_0_0_xx[k] * cd_y[k] + g_x_0_0_xxy[k];

            g_x_0_y_xy[k] = -g_x_0_0_xy[k] * cd_y[k] + g_x_0_0_xyy[k];

            g_x_0_y_xz[k] = -g_x_0_0_xz[k] * cd_y[k] + g_x_0_0_xyz[k];

            g_x_0_y_yy[k] = -g_x_0_0_yy[k] * cd_y[k] + g_x_0_0_yyy[k];

            g_x_0_y_yz[k] = -g_x_0_0_yz[k] * cd_y[k] + g_x_0_0_yyz[k];

            g_x_0_y_zz[k] = -g_x_0_0_zz[k] * cd_y[k] + g_x_0_0_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_x_0_z_xx = cbuffer.data(pd_geom_10_off + 0 * acomps  + 12);

        auto g_x_0_z_xy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 13);

        auto g_x_0_z_xz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 14);

        auto g_x_0_z_yy = cbuffer.data(pd_geom_10_off + 0 * acomps  + 15);

        auto g_x_0_z_yz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 16);

        auto g_x_0_z_zz = cbuffer.data(pd_geom_10_off + 0 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_x_0_0_xx, g_x_0_0_xxz, g_x_0_0_xy, g_x_0_0_xyz, g_x_0_0_xz, g_x_0_0_xzz, g_x_0_0_yy, g_x_0_0_yyz, g_x_0_0_yz, g_x_0_0_yzz, g_x_0_0_zz, g_x_0_0_zzz, g_x_0_z_xx, g_x_0_z_xy, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yz, g_x_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_x_0_z_xx[k] = -g_x_0_0_xx[k] * cd_z[k] + g_x_0_0_xxz[k];

            g_x_0_z_xy[k] = -g_x_0_0_xy[k] * cd_z[k] + g_x_0_0_xyz[k];

            g_x_0_z_xz[k] = -g_x_0_0_xz[k] * cd_z[k] + g_x_0_0_xzz[k];

            g_x_0_z_yy[k] = -g_x_0_0_yy[k] * cd_z[k] + g_x_0_0_yyz[k];

            g_x_0_z_yz[k] = -g_x_0_0_yz[k] * cd_z[k] + g_x_0_0_yzz[k];

            g_x_0_z_zz[k] = -g_x_0_0_zz[k] * cd_z[k] + g_x_0_0_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_y_0_x_xx = cbuffer.data(pd_geom_10_off + 18 * acomps  + 0);

        auto g_y_0_x_xy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 1);

        auto g_y_0_x_xz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 2);

        auto g_y_0_x_yy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 3);

        auto g_y_0_x_yz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 4);

        auto g_y_0_x_zz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_y_0_0_xx, g_y_0_0_xxx, g_y_0_0_xxy, g_y_0_0_xxz, g_y_0_0_xy, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xz, g_y_0_0_xzz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_y_0_x_xx, g_y_0_x_xy, g_y_0_x_xz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_x_xx[k] = -g_y_0_0_xx[k] * cd_x[k] + g_y_0_0_xxx[k];

            g_y_0_x_xy[k] = -g_y_0_0_xy[k] * cd_x[k] + g_y_0_0_xxy[k];

            g_y_0_x_xz[k] = -g_y_0_0_xz[k] * cd_x[k] + g_y_0_0_xxz[k];

            g_y_0_x_yy[k] = -g_y_0_0_yy[k] * cd_x[k] + g_y_0_0_xyy[k];

            g_y_0_x_yz[k] = -g_y_0_0_yz[k] * cd_x[k] + g_y_0_0_xyz[k];

            g_y_0_x_zz[k] = -g_y_0_0_zz[k] * cd_x[k] + g_y_0_0_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_y_0_y_xx = cbuffer.data(pd_geom_10_off + 18 * acomps  + 6);

        auto g_y_0_y_xy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 7);

        auto g_y_0_y_xz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 8);

        auto g_y_0_y_yy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 9);

        auto g_y_0_y_yz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 10);

        auto g_y_0_y_zz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz, g_0_zz, g_y_0_0_xx, g_y_0_0_xxy, g_y_0_0_xy, g_y_0_0_xyy, g_y_0_0_xyz, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yyy, g_y_0_0_yyz, g_y_0_0_yz, g_y_0_0_yzz, g_y_0_0_zz, g_y_0_y_xx, g_y_0_y_xy, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_y_xx[k] = -g_0_xx[k] - g_y_0_0_xx[k] * cd_y[k] + g_y_0_0_xxy[k];

            g_y_0_y_xy[k] = -g_0_xy[k] - g_y_0_0_xy[k] * cd_y[k] + g_y_0_0_xyy[k];

            g_y_0_y_xz[k] = -g_0_xz[k] - g_y_0_0_xz[k] * cd_y[k] + g_y_0_0_xyz[k];

            g_y_0_y_yy[k] = -g_0_yy[k] - g_y_0_0_yy[k] * cd_y[k] + g_y_0_0_yyy[k];

            g_y_0_y_yz[k] = -g_0_yz[k] - g_y_0_0_yz[k] * cd_y[k] + g_y_0_0_yyz[k];

            g_y_0_y_zz[k] = -g_0_zz[k] - g_y_0_0_zz[k] * cd_y[k] + g_y_0_0_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_y_0_z_xx = cbuffer.data(pd_geom_10_off + 18 * acomps  + 12);

        auto g_y_0_z_xy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 13);

        auto g_y_0_z_xz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 14);

        auto g_y_0_z_yy = cbuffer.data(pd_geom_10_off + 18 * acomps  + 15);

        auto g_y_0_z_yz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 16);

        auto g_y_0_z_zz = cbuffer.data(pd_geom_10_off + 18 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_y_0_0_xx, g_y_0_0_xxz, g_y_0_0_xy, g_y_0_0_xyz, g_y_0_0_xz, g_y_0_0_xzz, g_y_0_0_yy, g_y_0_0_yyz, g_y_0_0_yz, g_y_0_0_yzz, g_y_0_0_zz, g_y_0_0_zzz, g_y_0_z_xx, g_y_0_z_xy, g_y_0_z_xz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_y_0_z_xx[k] = -g_y_0_0_xx[k] * cd_z[k] + g_y_0_0_xxz[k];

            g_y_0_z_xy[k] = -g_y_0_0_xy[k] * cd_z[k] + g_y_0_0_xyz[k];

            g_y_0_z_xz[k] = -g_y_0_0_xz[k] * cd_z[k] + g_y_0_0_xzz[k];

            g_y_0_z_yy[k] = -g_y_0_0_yy[k] * cd_z[k] + g_y_0_0_yyz[k];

            g_y_0_z_yz[k] = -g_y_0_0_yz[k] * cd_z[k] + g_y_0_0_yzz[k];

            g_y_0_z_zz[k] = -g_y_0_0_zz[k] * cd_z[k] + g_y_0_0_zzz[k];
        }
        /// Set up 0-6 components of targeted buffer : cbuffer.data(

        auto g_z_0_x_xx = cbuffer.data(pd_geom_10_off + 36 * acomps  + 0);

        auto g_z_0_x_xy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 1);

        auto g_z_0_x_xz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 2);

        auto g_z_0_x_yy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 3);

        auto g_z_0_x_yz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 4);

        auto g_z_0_x_zz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 5);

        #pragma omp simd aligned(cd_x, g_z_0_0_xx, g_z_0_0_xxx, g_z_0_0_xxy, g_z_0_0_xxz, g_z_0_0_xy, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xz, g_z_0_0_xzz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_z_0_x_xx, g_z_0_x_xy, g_z_0_x_xz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_x_xx[k] = -g_z_0_0_xx[k] * cd_x[k] + g_z_0_0_xxx[k];

            g_z_0_x_xy[k] = -g_z_0_0_xy[k] * cd_x[k] + g_z_0_0_xxy[k];

            g_z_0_x_xz[k] = -g_z_0_0_xz[k] * cd_x[k] + g_z_0_0_xxz[k];

            g_z_0_x_yy[k] = -g_z_0_0_yy[k] * cd_x[k] + g_z_0_0_xyy[k];

            g_z_0_x_yz[k] = -g_z_0_0_yz[k] * cd_x[k] + g_z_0_0_xyz[k];

            g_z_0_x_zz[k] = -g_z_0_0_zz[k] * cd_x[k] + g_z_0_0_xzz[k];
        }

        /// Set up 6-12 components of targeted buffer : cbuffer.data(

        auto g_z_0_y_xx = cbuffer.data(pd_geom_10_off + 36 * acomps  + 6);

        auto g_z_0_y_xy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 7);

        auto g_z_0_y_xz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 8);

        auto g_z_0_y_yy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 9);

        auto g_z_0_y_yz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 10);

        auto g_z_0_y_zz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 11);

        #pragma omp simd aligned(cd_y, g_z_0_0_xx, g_z_0_0_xxy, g_z_0_0_xy, g_z_0_0_xyy, g_z_0_0_xyz, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yyy, g_z_0_0_yyz, g_z_0_0_yz, g_z_0_0_yzz, g_z_0_0_zz, g_z_0_y_xx, g_z_0_y_xy, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_y_xx[k] = -g_z_0_0_xx[k] * cd_y[k] + g_z_0_0_xxy[k];

            g_z_0_y_xy[k] = -g_z_0_0_xy[k] * cd_y[k] + g_z_0_0_xyy[k];

            g_z_0_y_xz[k] = -g_z_0_0_xz[k] * cd_y[k] + g_z_0_0_xyz[k];

            g_z_0_y_yy[k] = -g_z_0_0_yy[k] * cd_y[k] + g_z_0_0_yyy[k];

            g_z_0_y_yz[k] = -g_z_0_0_yz[k] * cd_y[k] + g_z_0_0_yyz[k];

            g_z_0_y_zz[k] = -g_z_0_0_zz[k] * cd_y[k] + g_z_0_0_yzz[k];
        }

        /// Set up 12-18 components of targeted buffer : cbuffer.data(

        auto g_z_0_z_xx = cbuffer.data(pd_geom_10_off + 36 * acomps  + 12);

        auto g_z_0_z_xy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 13);

        auto g_z_0_z_xz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 14);

        auto g_z_0_z_yy = cbuffer.data(pd_geom_10_off + 36 * acomps  + 15);

        auto g_z_0_z_yz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 16);

        auto g_z_0_z_zz = cbuffer.data(pd_geom_10_off + 36 * acomps  + 17);

        #pragma omp simd aligned(cd_z, g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz, g_0_zz, g_z_0_0_xx, g_z_0_0_xxz, g_z_0_0_xy, g_z_0_0_xyz, g_z_0_0_xz, g_z_0_0_xzz, g_z_0_0_yy, g_z_0_0_yyz, g_z_0_0_yz, g_z_0_0_yzz, g_z_0_0_zz, g_z_0_0_zzz, g_z_0_z_xx, g_z_0_z_xy, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz  : 64)
        for (size_t k = 0; k < nelems; k++)
        {
            g_z_0_z_xx[k] = -g_0_xx[k] - g_z_0_0_xx[k] * cd_z[k] + g_z_0_0_xxz[k];

            g_z_0_z_xy[k] = -g_0_xy[k] - g_z_0_0_xy[k] * cd_z[k] + g_z_0_0_xyz[k];

            g_z_0_z_xz[k] = -g_0_xz[k] - g_z_0_0_xz[k] * cd_z[k] + g_z_0_0_xzz[k];

            g_z_0_z_yy[k] = -g_0_yy[k] - g_z_0_0_yy[k] * cd_z[k] + g_z_0_0_yyz[k];

            g_z_0_z_yz[k] = -g_0_yz[k] - g_z_0_0_yz[k] * cd_z[k] + g_z_0_0_yzz[k];

            g_z_0_z_zz[k] = -g_0_zz[k] - g_z_0_0_zz[k] * cd_z[k] + g_z_0_0_zzz[k];
        }
    }
}

} // t3ceri namespace

