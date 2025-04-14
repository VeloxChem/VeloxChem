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

#include "ElectronRepulsionGeom2000ContrRecPDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_pdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_pdxx,
                                            const size_t idx_geom_10_sdxx,
                                            const size_t idx_geom_20_sdxx,
                                            const size_t idx_geom_20_sfxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_10_off = idx_geom_10_sdxx + i * dcomps + j;

            auto g_x_0_0_xx = cbuffer.data(sd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xy = cbuffer.data(sd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xz = cbuffer.data(sd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_yy = cbuffer.data(sd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_yz = cbuffer.data(sd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_zz = cbuffer.data(sd_geom_10_off + 5 * ccomps * dcomps);

            auto g_y_0_0_xx = cbuffer.data(sd_geom_10_off + 6 * ccomps * dcomps);

            auto g_y_0_0_xy = cbuffer.data(sd_geom_10_off + 7 * ccomps * dcomps);

            auto g_y_0_0_xz = cbuffer.data(sd_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_0_yy = cbuffer.data(sd_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_0_yz = cbuffer.data(sd_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_0_zz = cbuffer.data(sd_geom_10_off + 11 * ccomps * dcomps);

            auto g_z_0_0_xx = cbuffer.data(sd_geom_10_off + 12 * ccomps * dcomps);

            auto g_z_0_0_xy = cbuffer.data(sd_geom_10_off + 13 * ccomps * dcomps);

            auto g_z_0_0_xz = cbuffer.data(sd_geom_10_off + 14 * ccomps * dcomps);

            auto g_z_0_0_yy = cbuffer.data(sd_geom_10_off + 15 * ccomps * dcomps);

            auto g_z_0_0_yz = cbuffer.data(sd_geom_10_off + 16 * ccomps * dcomps);

            auto g_z_0_0_zz = cbuffer.data(sd_geom_10_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_20_off = idx_geom_20_sdxx + i * dcomps + j;

            auto g_xx_0_0_xx = cbuffer.data(sd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xy = cbuffer.data(sd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xz = cbuffer.data(sd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_yy = cbuffer.data(sd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_yz = cbuffer.data(sd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_zz = cbuffer.data(sd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xy_0_0_xx = cbuffer.data(sd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xy_0_0_xy = cbuffer.data(sd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xy_0_0_xz = cbuffer.data(sd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xy_0_0_yy = cbuffer.data(sd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xy_0_0_yz = cbuffer.data(sd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xy_0_0_zz = cbuffer.data(sd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xz_0_0_xx = cbuffer.data(sd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xz_0_0_xy = cbuffer.data(sd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xz_0_0_xz = cbuffer.data(sd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xz_0_0_yy = cbuffer.data(sd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xz_0_0_yz = cbuffer.data(sd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xz_0_0_zz = cbuffer.data(sd_geom_20_off + 17 * ccomps * dcomps);

            auto g_yy_0_0_xx = cbuffer.data(sd_geom_20_off + 18 * ccomps * dcomps);

            auto g_yy_0_0_xy = cbuffer.data(sd_geom_20_off + 19 * ccomps * dcomps);

            auto g_yy_0_0_xz = cbuffer.data(sd_geom_20_off + 20 * ccomps * dcomps);

            auto g_yy_0_0_yy = cbuffer.data(sd_geom_20_off + 21 * ccomps * dcomps);

            auto g_yy_0_0_yz = cbuffer.data(sd_geom_20_off + 22 * ccomps * dcomps);

            auto g_yy_0_0_zz = cbuffer.data(sd_geom_20_off + 23 * ccomps * dcomps);

            auto g_yz_0_0_xx = cbuffer.data(sd_geom_20_off + 24 * ccomps * dcomps);

            auto g_yz_0_0_xy = cbuffer.data(sd_geom_20_off + 25 * ccomps * dcomps);

            auto g_yz_0_0_xz = cbuffer.data(sd_geom_20_off + 26 * ccomps * dcomps);

            auto g_yz_0_0_yy = cbuffer.data(sd_geom_20_off + 27 * ccomps * dcomps);

            auto g_yz_0_0_yz = cbuffer.data(sd_geom_20_off + 28 * ccomps * dcomps);

            auto g_yz_0_0_zz = cbuffer.data(sd_geom_20_off + 29 * ccomps * dcomps);

            auto g_zz_0_0_xx = cbuffer.data(sd_geom_20_off + 30 * ccomps * dcomps);

            auto g_zz_0_0_xy = cbuffer.data(sd_geom_20_off + 31 * ccomps * dcomps);

            auto g_zz_0_0_xz = cbuffer.data(sd_geom_20_off + 32 * ccomps * dcomps);

            auto g_zz_0_0_yy = cbuffer.data(sd_geom_20_off + 33 * ccomps * dcomps);

            auto g_zz_0_0_yz = cbuffer.data(sd_geom_20_off + 34 * ccomps * dcomps);

            auto g_zz_0_0_zz = cbuffer.data(sd_geom_20_off + 35 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_20_off = idx_geom_20_sfxx + i * dcomps + j;

            auto g_xx_0_0_xxx = cbuffer.data(sf_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxy = cbuffer.data(sf_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxz = cbuffer.data(sf_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xyy = cbuffer.data(sf_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xyz = cbuffer.data(sf_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xzz = cbuffer.data(sf_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_yyy = cbuffer.data(sf_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_yyz = cbuffer.data(sf_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_yzz = cbuffer.data(sf_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_zzz = cbuffer.data(sf_geom_20_off + 9 * ccomps * dcomps);

            auto g_xy_0_0_xxx = cbuffer.data(sf_geom_20_off + 10 * ccomps * dcomps);

            auto g_xy_0_0_xxy = cbuffer.data(sf_geom_20_off + 11 * ccomps * dcomps);

            auto g_xy_0_0_xxz = cbuffer.data(sf_geom_20_off + 12 * ccomps * dcomps);

            auto g_xy_0_0_xyy = cbuffer.data(sf_geom_20_off + 13 * ccomps * dcomps);

            auto g_xy_0_0_xyz = cbuffer.data(sf_geom_20_off + 14 * ccomps * dcomps);

            auto g_xy_0_0_xzz = cbuffer.data(sf_geom_20_off + 15 * ccomps * dcomps);

            auto g_xy_0_0_yyy = cbuffer.data(sf_geom_20_off + 16 * ccomps * dcomps);

            auto g_xy_0_0_yyz = cbuffer.data(sf_geom_20_off + 17 * ccomps * dcomps);

            auto g_xy_0_0_yzz = cbuffer.data(sf_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_0_zzz = cbuffer.data(sf_geom_20_off + 19 * ccomps * dcomps);

            auto g_xz_0_0_xxx = cbuffer.data(sf_geom_20_off + 20 * ccomps * dcomps);

            auto g_xz_0_0_xxy = cbuffer.data(sf_geom_20_off + 21 * ccomps * dcomps);

            auto g_xz_0_0_xxz = cbuffer.data(sf_geom_20_off + 22 * ccomps * dcomps);

            auto g_xz_0_0_xyy = cbuffer.data(sf_geom_20_off + 23 * ccomps * dcomps);

            auto g_xz_0_0_xyz = cbuffer.data(sf_geom_20_off + 24 * ccomps * dcomps);

            auto g_xz_0_0_xzz = cbuffer.data(sf_geom_20_off + 25 * ccomps * dcomps);

            auto g_xz_0_0_yyy = cbuffer.data(sf_geom_20_off + 26 * ccomps * dcomps);

            auto g_xz_0_0_yyz = cbuffer.data(sf_geom_20_off + 27 * ccomps * dcomps);

            auto g_xz_0_0_yzz = cbuffer.data(sf_geom_20_off + 28 * ccomps * dcomps);

            auto g_xz_0_0_zzz = cbuffer.data(sf_geom_20_off + 29 * ccomps * dcomps);

            auto g_yy_0_0_xxx = cbuffer.data(sf_geom_20_off + 30 * ccomps * dcomps);

            auto g_yy_0_0_xxy = cbuffer.data(sf_geom_20_off + 31 * ccomps * dcomps);

            auto g_yy_0_0_xxz = cbuffer.data(sf_geom_20_off + 32 * ccomps * dcomps);

            auto g_yy_0_0_xyy = cbuffer.data(sf_geom_20_off + 33 * ccomps * dcomps);

            auto g_yy_0_0_xyz = cbuffer.data(sf_geom_20_off + 34 * ccomps * dcomps);

            auto g_yy_0_0_xzz = cbuffer.data(sf_geom_20_off + 35 * ccomps * dcomps);

            auto g_yy_0_0_yyy = cbuffer.data(sf_geom_20_off + 36 * ccomps * dcomps);

            auto g_yy_0_0_yyz = cbuffer.data(sf_geom_20_off + 37 * ccomps * dcomps);

            auto g_yy_0_0_yzz = cbuffer.data(sf_geom_20_off + 38 * ccomps * dcomps);

            auto g_yy_0_0_zzz = cbuffer.data(sf_geom_20_off + 39 * ccomps * dcomps);

            auto g_yz_0_0_xxx = cbuffer.data(sf_geom_20_off + 40 * ccomps * dcomps);

            auto g_yz_0_0_xxy = cbuffer.data(sf_geom_20_off + 41 * ccomps * dcomps);

            auto g_yz_0_0_xxz = cbuffer.data(sf_geom_20_off + 42 * ccomps * dcomps);

            auto g_yz_0_0_xyy = cbuffer.data(sf_geom_20_off + 43 * ccomps * dcomps);

            auto g_yz_0_0_xyz = cbuffer.data(sf_geom_20_off + 44 * ccomps * dcomps);

            auto g_yz_0_0_xzz = cbuffer.data(sf_geom_20_off + 45 * ccomps * dcomps);

            auto g_yz_0_0_yyy = cbuffer.data(sf_geom_20_off + 46 * ccomps * dcomps);

            auto g_yz_0_0_yyz = cbuffer.data(sf_geom_20_off + 47 * ccomps * dcomps);

            auto g_yz_0_0_yzz = cbuffer.data(sf_geom_20_off + 48 * ccomps * dcomps);

            auto g_yz_0_0_zzz = cbuffer.data(sf_geom_20_off + 49 * ccomps * dcomps);

            auto g_zz_0_0_xxx = cbuffer.data(sf_geom_20_off + 50 * ccomps * dcomps);

            auto g_zz_0_0_xxy = cbuffer.data(sf_geom_20_off + 51 * ccomps * dcomps);

            auto g_zz_0_0_xxz = cbuffer.data(sf_geom_20_off + 52 * ccomps * dcomps);

            auto g_zz_0_0_xyy = cbuffer.data(sf_geom_20_off + 53 * ccomps * dcomps);

            auto g_zz_0_0_xyz = cbuffer.data(sf_geom_20_off + 54 * ccomps * dcomps);

            auto g_zz_0_0_xzz = cbuffer.data(sf_geom_20_off + 55 * ccomps * dcomps);

            auto g_zz_0_0_yyy = cbuffer.data(sf_geom_20_off + 56 * ccomps * dcomps);

            auto g_zz_0_0_yyz = cbuffer.data(sf_geom_20_off + 57 * ccomps * dcomps);

            auto g_zz_0_0_yzz = cbuffer.data(sf_geom_20_off + 58 * ccomps * dcomps);

            auto g_zz_0_0_zzz = cbuffer.data(sf_geom_20_off + 59 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pdxx

            const auto pd_geom_20_off = idx_geom_20_pdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_x_xx = cbuffer.data(pd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xy = cbuffer.data(pd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xz = cbuffer.data(pd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_yy = cbuffer.data(pd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_yz = cbuffer.data(pd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_zz = cbuffer.data(pd_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_xx_0_0_xx, g_xx_0_0_xxx, g_xx_0_0_xxy, g_xx_0_0_xxz, g_xx_0_0_xy, g_xx_0_0_xyy, g_xx_0_0_xyz, g_xx_0_0_xz, g_xx_0_0_xzz, g_xx_0_0_yy, g_xx_0_0_yz, g_xx_0_0_zz, g_xx_0_x_xx, g_xx_0_x_xy, g_xx_0_x_xz, g_xx_0_x_yy, g_xx_0_x_yz, g_xx_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_x_xx[k] = -2.0 * g_x_0_0_xx[k] - g_xx_0_0_xx[k] * ab_x + g_xx_0_0_xxx[k];

                g_xx_0_x_xy[k] = -2.0 * g_x_0_0_xy[k] - g_xx_0_0_xy[k] * ab_x + g_xx_0_0_xxy[k];

                g_xx_0_x_xz[k] = -2.0 * g_x_0_0_xz[k] - g_xx_0_0_xz[k] * ab_x + g_xx_0_0_xxz[k];

                g_xx_0_x_yy[k] = -2.0 * g_x_0_0_yy[k] - g_xx_0_0_yy[k] * ab_x + g_xx_0_0_xyy[k];

                g_xx_0_x_yz[k] = -2.0 * g_x_0_0_yz[k] - g_xx_0_0_yz[k] * ab_x + g_xx_0_0_xyz[k];

                g_xx_0_x_zz[k] = -2.0 * g_x_0_0_zz[k] - g_xx_0_0_zz[k] * ab_x + g_xx_0_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_y_xx = cbuffer.data(pd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_y_xy = cbuffer.data(pd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_y_xz = cbuffer.data(pd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_y_yy = cbuffer.data(pd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_y_yz = cbuffer.data(pd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_y_zz = cbuffer.data(pd_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xx, g_xx_0_0_xxy, g_xx_0_0_xy, g_xx_0_0_xyy, g_xx_0_0_xyz, g_xx_0_0_xz, g_xx_0_0_yy, g_xx_0_0_yyy, g_xx_0_0_yyz, g_xx_0_0_yz, g_xx_0_0_yzz, g_xx_0_0_zz, g_xx_0_y_xx, g_xx_0_y_xy, g_xx_0_y_xz, g_xx_0_y_yy, g_xx_0_y_yz, g_xx_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_y_xx[k] = -g_xx_0_0_xx[k] * ab_y + g_xx_0_0_xxy[k];

                g_xx_0_y_xy[k] = -g_xx_0_0_xy[k] * ab_y + g_xx_0_0_xyy[k];

                g_xx_0_y_xz[k] = -g_xx_0_0_xz[k] * ab_y + g_xx_0_0_xyz[k];

                g_xx_0_y_yy[k] = -g_xx_0_0_yy[k] * ab_y + g_xx_0_0_yyy[k];

                g_xx_0_y_yz[k] = -g_xx_0_0_yz[k] * ab_y + g_xx_0_0_yyz[k];

                g_xx_0_y_zz[k] = -g_xx_0_0_zz[k] * ab_y + g_xx_0_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_z_xx = cbuffer.data(pd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_z_xy = cbuffer.data(pd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_z_xz = cbuffer.data(pd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_z_yy = cbuffer.data(pd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_z_yz = cbuffer.data(pd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_z_zz = cbuffer.data(pd_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xx, g_xx_0_0_xxz, g_xx_0_0_xy, g_xx_0_0_xyz, g_xx_0_0_xz, g_xx_0_0_xzz, g_xx_0_0_yy, g_xx_0_0_yyz, g_xx_0_0_yz, g_xx_0_0_yzz, g_xx_0_0_zz, g_xx_0_0_zzz, g_xx_0_z_xx, g_xx_0_z_xy, g_xx_0_z_xz, g_xx_0_z_yy, g_xx_0_z_yz, g_xx_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_z_xx[k] = -g_xx_0_0_xx[k] * ab_z + g_xx_0_0_xxz[k];

                g_xx_0_z_xy[k] = -g_xx_0_0_xy[k] * ab_z + g_xx_0_0_xyz[k];

                g_xx_0_z_xz[k] = -g_xx_0_0_xz[k] * ab_z + g_xx_0_0_xzz[k];

                g_xx_0_z_yy[k] = -g_xx_0_0_yy[k] * ab_z + g_xx_0_0_yyz[k];

                g_xx_0_z_yz[k] = -g_xx_0_0_yz[k] * ab_z + g_xx_0_0_yzz[k];

                g_xx_0_z_zz[k] = -g_xx_0_0_zz[k] * ab_z + g_xx_0_0_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_xy_0_x_xx = cbuffer.data(pd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xy_0_x_xy = cbuffer.data(pd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xy_0_x_xz = cbuffer.data(pd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xy_0_x_yy = cbuffer.data(pd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xy_0_x_yz = cbuffer.data(pd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xy_0_x_zz = cbuffer.data(pd_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xx, g_xy_0_0_xxx, g_xy_0_0_xxy, g_xy_0_0_xxz, g_xy_0_0_xy, g_xy_0_0_xyy, g_xy_0_0_xyz, g_xy_0_0_xz, g_xy_0_0_xzz, g_xy_0_0_yy, g_xy_0_0_yz, g_xy_0_0_zz, g_xy_0_x_xx, g_xy_0_x_xy, g_xy_0_x_xz, g_xy_0_x_yy, g_xy_0_x_yz, g_xy_0_x_zz, g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_x_xx[k] = -g_y_0_0_xx[k] - g_xy_0_0_xx[k] * ab_x + g_xy_0_0_xxx[k];

                g_xy_0_x_xy[k] = -g_y_0_0_xy[k] - g_xy_0_0_xy[k] * ab_x + g_xy_0_0_xxy[k];

                g_xy_0_x_xz[k] = -g_y_0_0_xz[k] - g_xy_0_0_xz[k] * ab_x + g_xy_0_0_xxz[k];

                g_xy_0_x_yy[k] = -g_y_0_0_yy[k] - g_xy_0_0_yy[k] * ab_x + g_xy_0_0_xyy[k];

                g_xy_0_x_yz[k] = -g_y_0_0_yz[k] - g_xy_0_0_yz[k] * ab_x + g_xy_0_0_xyz[k];

                g_xy_0_x_zz[k] = -g_y_0_0_zz[k] - g_xy_0_0_zz[k] * ab_x + g_xy_0_0_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_xy_0_y_xx = cbuffer.data(pd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xy_0_y_xy = cbuffer.data(pd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xy_0_y_xz = cbuffer.data(pd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xy_0_y_yy = cbuffer.data(pd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_y_yz = cbuffer.data(pd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_y_zz = cbuffer.data(pd_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_xy_0_0_xx, g_xy_0_0_xxy, g_xy_0_0_xy, g_xy_0_0_xyy, g_xy_0_0_xyz, g_xy_0_0_xz, g_xy_0_0_yy, g_xy_0_0_yyy, g_xy_0_0_yyz, g_xy_0_0_yz, g_xy_0_0_yzz, g_xy_0_0_zz, g_xy_0_y_xx, g_xy_0_y_xy, g_xy_0_y_xz, g_xy_0_y_yy, g_xy_0_y_yz, g_xy_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_y_xx[k] = -g_x_0_0_xx[k] - g_xy_0_0_xx[k] * ab_y + g_xy_0_0_xxy[k];

                g_xy_0_y_xy[k] = -g_x_0_0_xy[k] - g_xy_0_0_xy[k] * ab_y + g_xy_0_0_xyy[k];

                g_xy_0_y_xz[k] = -g_x_0_0_xz[k] - g_xy_0_0_xz[k] * ab_y + g_xy_0_0_xyz[k];

                g_xy_0_y_yy[k] = -g_x_0_0_yy[k] - g_xy_0_0_yy[k] * ab_y + g_xy_0_0_yyy[k];

                g_xy_0_y_yz[k] = -g_x_0_0_yz[k] - g_xy_0_0_yz[k] * ab_y + g_xy_0_0_yyz[k];

                g_xy_0_y_zz[k] = -g_x_0_0_zz[k] - g_xy_0_0_zz[k] * ab_y + g_xy_0_0_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_xy_0_z_xx = cbuffer.data(pd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_z_xy = cbuffer.data(pd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_z_xz = cbuffer.data(pd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_z_yy = cbuffer.data(pd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_z_yz = cbuffer.data(pd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_z_zz = cbuffer.data(pd_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xx, g_xy_0_0_xxz, g_xy_0_0_xy, g_xy_0_0_xyz, g_xy_0_0_xz, g_xy_0_0_xzz, g_xy_0_0_yy, g_xy_0_0_yyz, g_xy_0_0_yz, g_xy_0_0_yzz, g_xy_0_0_zz, g_xy_0_0_zzz, g_xy_0_z_xx, g_xy_0_z_xy, g_xy_0_z_xz, g_xy_0_z_yy, g_xy_0_z_yz, g_xy_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_z_xx[k] = -g_xy_0_0_xx[k] * ab_z + g_xy_0_0_xxz[k];

                g_xy_0_z_xy[k] = -g_xy_0_0_xy[k] * ab_z + g_xy_0_0_xyz[k];

                g_xy_0_z_xz[k] = -g_xy_0_0_xz[k] * ab_z + g_xy_0_0_xzz[k];

                g_xy_0_z_yy[k] = -g_xy_0_0_yy[k] * ab_z + g_xy_0_0_yyz[k];

                g_xy_0_z_yz[k] = -g_xy_0_0_yz[k] * ab_z + g_xy_0_0_yzz[k];

                g_xy_0_z_zz[k] = -g_xy_0_0_zz[k] * ab_z + g_xy_0_0_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_xz_0_x_xx = cbuffer.data(pd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xz_0_x_xy = cbuffer.data(pd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xz_0_x_xz = cbuffer.data(pd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xz_0_x_yy = cbuffer.data(pd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xz_0_x_yz = cbuffer.data(pd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xz_0_x_zz = cbuffer.data(pd_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xx, g_xz_0_0_xxx, g_xz_0_0_xxy, g_xz_0_0_xxz, g_xz_0_0_xy, g_xz_0_0_xyy, g_xz_0_0_xyz, g_xz_0_0_xz, g_xz_0_0_xzz, g_xz_0_0_yy, g_xz_0_0_yz, g_xz_0_0_zz, g_xz_0_x_xx, g_xz_0_x_xy, g_xz_0_x_xz, g_xz_0_x_yy, g_xz_0_x_yz, g_xz_0_x_zz, g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_x_xx[k] = -g_z_0_0_xx[k] - g_xz_0_0_xx[k] * ab_x + g_xz_0_0_xxx[k];

                g_xz_0_x_xy[k] = -g_z_0_0_xy[k] - g_xz_0_0_xy[k] * ab_x + g_xz_0_0_xxy[k];

                g_xz_0_x_xz[k] = -g_z_0_0_xz[k] - g_xz_0_0_xz[k] * ab_x + g_xz_0_0_xxz[k];

                g_xz_0_x_yy[k] = -g_z_0_0_yy[k] - g_xz_0_0_yy[k] * ab_x + g_xz_0_0_xyy[k];

                g_xz_0_x_yz[k] = -g_z_0_0_yz[k] - g_xz_0_0_yz[k] * ab_x + g_xz_0_0_xyz[k];

                g_xz_0_x_zz[k] = -g_z_0_0_zz[k] - g_xz_0_0_zz[k] * ab_x + g_xz_0_0_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_xz_0_y_xx = cbuffer.data(pd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xz_0_y_xy = cbuffer.data(pd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xz_0_y_xz = cbuffer.data(pd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xz_0_y_yy = cbuffer.data(pd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xz_0_y_yz = cbuffer.data(pd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xz_0_y_zz = cbuffer.data(pd_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xx, g_xz_0_0_xxy, g_xz_0_0_xy, g_xz_0_0_xyy, g_xz_0_0_xyz, g_xz_0_0_xz, g_xz_0_0_yy, g_xz_0_0_yyy, g_xz_0_0_yyz, g_xz_0_0_yz, g_xz_0_0_yzz, g_xz_0_0_zz, g_xz_0_y_xx, g_xz_0_y_xy, g_xz_0_y_xz, g_xz_0_y_yy, g_xz_0_y_yz, g_xz_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_y_xx[k] = -g_xz_0_0_xx[k] * ab_y + g_xz_0_0_xxy[k];

                g_xz_0_y_xy[k] = -g_xz_0_0_xy[k] * ab_y + g_xz_0_0_xyy[k];

                g_xz_0_y_xz[k] = -g_xz_0_0_xz[k] * ab_y + g_xz_0_0_xyz[k];

                g_xz_0_y_yy[k] = -g_xz_0_0_yy[k] * ab_y + g_xz_0_0_yyy[k];

                g_xz_0_y_yz[k] = -g_xz_0_0_yz[k] * ab_y + g_xz_0_0_yyz[k];

                g_xz_0_y_zz[k] = -g_xz_0_0_zz[k] * ab_y + g_xz_0_0_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_xz_0_z_xx = cbuffer.data(pd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xz_0_z_xy = cbuffer.data(pd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xz_0_z_xz = cbuffer.data(pd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xz_0_z_yy = cbuffer.data(pd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xz_0_z_yz = cbuffer.data(pd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xz_0_z_zz = cbuffer.data(pd_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xx, g_x_0_0_xy, g_x_0_0_xz, g_x_0_0_yy, g_x_0_0_yz, g_x_0_0_zz, g_xz_0_0_xx, g_xz_0_0_xxz, g_xz_0_0_xy, g_xz_0_0_xyz, g_xz_0_0_xz, g_xz_0_0_xzz, g_xz_0_0_yy, g_xz_0_0_yyz, g_xz_0_0_yz, g_xz_0_0_yzz, g_xz_0_0_zz, g_xz_0_0_zzz, g_xz_0_z_xx, g_xz_0_z_xy, g_xz_0_z_xz, g_xz_0_z_yy, g_xz_0_z_yz, g_xz_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_z_xx[k] = -g_x_0_0_xx[k] - g_xz_0_0_xx[k] * ab_z + g_xz_0_0_xxz[k];

                g_xz_0_z_xy[k] = -g_x_0_0_xy[k] - g_xz_0_0_xy[k] * ab_z + g_xz_0_0_xyz[k];

                g_xz_0_z_xz[k] = -g_x_0_0_xz[k] - g_xz_0_0_xz[k] * ab_z + g_xz_0_0_xzz[k];

                g_xz_0_z_yy[k] = -g_x_0_0_yy[k] - g_xz_0_0_yy[k] * ab_z + g_xz_0_0_yyz[k];

                g_xz_0_z_yz[k] = -g_x_0_0_yz[k] - g_xz_0_0_yz[k] * ab_z + g_xz_0_0_yzz[k];

                g_xz_0_z_zz[k] = -g_x_0_0_zz[k] - g_xz_0_0_zz[k] * ab_z + g_xz_0_0_zzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_yy_0_x_xx = cbuffer.data(pd_geom_20_off + 54 * ccomps * dcomps);

            auto g_yy_0_x_xy = cbuffer.data(pd_geom_20_off + 55 * ccomps * dcomps);

            auto g_yy_0_x_xz = cbuffer.data(pd_geom_20_off + 56 * ccomps * dcomps);

            auto g_yy_0_x_yy = cbuffer.data(pd_geom_20_off + 57 * ccomps * dcomps);

            auto g_yy_0_x_yz = cbuffer.data(pd_geom_20_off + 58 * ccomps * dcomps);

            auto g_yy_0_x_zz = cbuffer.data(pd_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xx, g_yy_0_0_xxx, g_yy_0_0_xxy, g_yy_0_0_xxz, g_yy_0_0_xy, g_yy_0_0_xyy, g_yy_0_0_xyz, g_yy_0_0_xz, g_yy_0_0_xzz, g_yy_0_0_yy, g_yy_0_0_yz, g_yy_0_0_zz, g_yy_0_x_xx, g_yy_0_x_xy, g_yy_0_x_xz, g_yy_0_x_yy, g_yy_0_x_yz, g_yy_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_x_xx[k] = -g_yy_0_0_xx[k] * ab_x + g_yy_0_0_xxx[k];

                g_yy_0_x_xy[k] = -g_yy_0_0_xy[k] * ab_x + g_yy_0_0_xxy[k];

                g_yy_0_x_xz[k] = -g_yy_0_0_xz[k] * ab_x + g_yy_0_0_xxz[k];

                g_yy_0_x_yy[k] = -g_yy_0_0_yy[k] * ab_x + g_yy_0_0_xyy[k];

                g_yy_0_x_yz[k] = -g_yy_0_0_yz[k] * ab_x + g_yy_0_0_xyz[k];

                g_yy_0_x_zz[k] = -g_yy_0_0_zz[k] * ab_x + g_yy_0_0_xzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_yy_0_y_xx = cbuffer.data(pd_geom_20_off + 60 * ccomps * dcomps);

            auto g_yy_0_y_xy = cbuffer.data(pd_geom_20_off + 61 * ccomps * dcomps);

            auto g_yy_0_y_xz = cbuffer.data(pd_geom_20_off + 62 * ccomps * dcomps);

            auto g_yy_0_y_yy = cbuffer.data(pd_geom_20_off + 63 * ccomps * dcomps);

            auto g_yy_0_y_yz = cbuffer.data(pd_geom_20_off + 64 * ccomps * dcomps);

            auto g_yy_0_y_zz = cbuffer.data(pd_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_yy_0_0_xx, g_yy_0_0_xxy, g_yy_0_0_xy, g_yy_0_0_xyy, g_yy_0_0_xyz, g_yy_0_0_xz, g_yy_0_0_yy, g_yy_0_0_yyy, g_yy_0_0_yyz, g_yy_0_0_yz, g_yy_0_0_yzz, g_yy_0_0_zz, g_yy_0_y_xx, g_yy_0_y_xy, g_yy_0_y_xz, g_yy_0_y_yy, g_yy_0_y_yz, g_yy_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_y_xx[k] = -2.0 * g_y_0_0_xx[k] - g_yy_0_0_xx[k] * ab_y + g_yy_0_0_xxy[k];

                g_yy_0_y_xy[k] = -2.0 * g_y_0_0_xy[k] - g_yy_0_0_xy[k] * ab_y + g_yy_0_0_xyy[k];

                g_yy_0_y_xz[k] = -2.0 * g_y_0_0_xz[k] - g_yy_0_0_xz[k] * ab_y + g_yy_0_0_xyz[k];

                g_yy_0_y_yy[k] = -2.0 * g_y_0_0_yy[k] - g_yy_0_0_yy[k] * ab_y + g_yy_0_0_yyy[k];

                g_yy_0_y_yz[k] = -2.0 * g_y_0_0_yz[k] - g_yy_0_0_yz[k] * ab_y + g_yy_0_0_yyz[k];

                g_yy_0_y_zz[k] = -2.0 * g_y_0_0_zz[k] - g_yy_0_0_zz[k] * ab_y + g_yy_0_0_yzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_yy_0_z_xx = cbuffer.data(pd_geom_20_off + 66 * ccomps * dcomps);

            auto g_yy_0_z_xy = cbuffer.data(pd_geom_20_off + 67 * ccomps * dcomps);

            auto g_yy_0_z_xz = cbuffer.data(pd_geom_20_off + 68 * ccomps * dcomps);

            auto g_yy_0_z_yy = cbuffer.data(pd_geom_20_off + 69 * ccomps * dcomps);

            auto g_yy_0_z_yz = cbuffer.data(pd_geom_20_off + 70 * ccomps * dcomps);

            auto g_yy_0_z_zz = cbuffer.data(pd_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xx, g_yy_0_0_xxz, g_yy_0_0_xy, g_yy_0_0_xyz, g_yy_0_0_xz, g_yy_0_0_xzz, g_yy_0_0_yy, g_yy_0_0_yyz, g_yy_0_0_yz, g_yy_0_0_yzz, g_yy_0_0_zz, g_yy_0_0_zzz, g_yy_0_z_xx, g_yy_0_z_xy, g_yy_0_z_xz, g_yy_0_z_yy, g_yy_0_z_yz, g_yy_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_z_xx[k] = -g_yy_0_0_xx[k] * ab_z + g_yy_0_0_xxz[k];

                g_yy_0_z_xy[k] = -g_yy_0_0_xy[k] * ab_z + g_yy_0_0_xyz[k];

                g_yy_0_z_xz[k] = -g_yy_0_0_xz[k] * ab_z + g_yy_0_0_xzz[k];

                g_yy_0_z_yy[k] = -g_yy_0_0_yy[k] * ab_z + g_yy_0_0_yyz[k];

                g_yy_0_z_yz[k] = -g_yy_0_0_yz[k] * ab_z + g_yy_0_0_yzz[k];

                g_yy_0_z_zz[k] = -g_yy_0_0_zz[k] * ab_z + g_yy_0_0_zzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_yz_0_x_xx = cbuffer.data(pd_geom_20_off + 72 * ccomps * dcomps);

            auto g_yz_0_x_xy = cbuffer.data(pd_geom_20_off + 73 * ccomps * dcomps);

            auto g_yz_0_x_xz = cbuffer.data(pd_geom_20_off + 74 * ccomps * dcomps);

            auto g_yz_0_x_yy = cbuffer.data(pd_geom_20_off + 75 * ccomps * dcomps);

            auto g_yz_0_x_yz = cbuffer.data(pd_geom_20_off + 76 * ccomps * dcomps);

            auto g_yz_0_x_zz = cbuffer.data(pd_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xx, g_yz_0_0_xxx, g_yz_0_0_xxy, g_yz_0_0_xxz, g_yz_0_0_xy, g_yz_0_0_xyy, g_yz_0_0_xyz, g_yz_0_0_xz, g_yz_0_0_xzz, g_yz_0_0_yy, g_yz_0_0_yz, g_yz_0_0_zz, g_yz_0_x_xx, g_yz_0_x_xy, g_yz_0_x_xz, g_yz_0_x_yy, g_yz_0_x_yz, g_yz_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_x_xx[k] = -g_yz_0_0_xx[k] * ab_x + g_yz_0_0_xxx[k];

                g_yz_0_x_xy[k] = -g_yz_0_0_xy[k] * ab_x + g_yz_0_0_xxy[k];

                g_yz_0_x_xz[k] = -g_yz_0_0_xz[k] * ab_x + g_yz_0_0_xxz[k];

                g_yz_0_x_yy[k] = -g_yz_0_0_yy[k] * ab_x + g_yz_0_0_xyy[k];

                g_yz_0_x_yz[k] = -g_yz_0_0_yz[k] * ab_x + g_yz_0_0_xyz[k];

                g_yz_0_x_zz[k] = -g_yz_0_0_zz[k] * ab_x + g_yz_0_0_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_yz_0_y_xx = cbuffer.data(pd_geom_20_off + 78 * ccomps * dcomps);

            auto g_yz_0_y_xy = cbuffer.data(pd_geom_20_off + 79 * ccomps * dcomps);

            auto g_yz_0_y_xz = cbuffer.data(pd_geom_20_off + 80 * ccomps * dcomps);

            auto g_yz_0_y_yy = cbuffer.data(pd_geom_20_off + 81 * ccomps * dcomps);

            auto g_yz_0_y_yz = cbuffer.data(pd_geom_20_off + 82 * ccomps * dcomps);

            auto g_yz_0_y_zz = cbuffer.data(pd_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xx, g_yz_0_0_xxy, g_yz_0_0_xy, g_yz_0_0_xyy, g_yz_0_0_xyz, g_yz_0_0_xz, g_yz_0_0_yy, g_yz_0_0_yyy, g_yz_0_0_yyz, g_yz_0_0_yz, g_yz_0_0_yzz, g_yz_0_0_zz, g_yz_0_y_xx, g_yz_0_y_xy, g_yz_0_y_xz, g_yz_0_y_yy, g_yz_0_y_yz, g_yz_0_y_zz, g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_y_xx[k] = -g_z_0_0_xx[k] - g_yz_0_0_xx[k] * ab_y + g_yz_0_0_xxy[k];

                g_yz_0_y_xy[k] = -g_z_0_0_xy[k] - g_yz_0_0_xy[k] * ab_y + g_yz_0_0_xyy[k];

                g_yz_0_y_xz[k] = -g_z_0_0_xz[k] - g_yz_0_0_xz[k] * ab_y + g_yz_0_0_xyz[k];

                g_yz_0_y_yy[k] = -g_z_0_0_yy[k] - g_yz_0_0_yy[k] * ab_y + g_yz_0_0_yyy[k];

                g_yz_0_y_yz[k] = -g_z_0_0_yz[k] - g_yz_0_0_yz[k] * ab_y + g_yz_0_0_yyz[k];

                g_yz_0_y_zz[k] = -g_z_0_0_zz[k] - g_yz_0_0_zz[k] * ab_y + g_yz_0_0_yzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_yz_0_z_xx = cbuffer.data(pd_geom_20_off + 84 * ccomps * dcomps);

            auto g_yz_0_z_xy = cbuffer.data(pd_geom_20_off + 85 * ccomps * dcomps);

            auto g_yz_0_z_xz = cbuffer.data(pd_geom_20_off + 86 * ccomps * dcomps);

            auto g_yz_0_z_yy = cbuffer.data(pd_geom_20_off + 87 * ccomps * dcomps);

            auto g_yz_0_z_yz = cbuffer.data(pd_geom_20_off + 88 * ccomps * dcomps);

            auto g_yz_0_z_zz = cbuffer.data(pd_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xx, g_y_0_0_xy, g_y_0_0_xz, g_y_0_0_yy, g_y_0_0_yz, g_y_0_0_zz, g_yz_0_0_xx, g_yz_0_0_xxz, g_yz_0_0_xy, g_yz_0_0_xyz, g_yz_0_0_xz, g_yz_0_0_xzz, g_yz_0_0_yy, g_yz_0_0_yyz, g_yz_0_0_yz, g_yz_0_0_yzz, g_yz_0_0_zz, g_yz_0_0_zzz, g_yz_0_z_xx, g_yz_0_z_xy, g_yz_0_z_xz, g_yz_0_z_yy, g_yz_0_z_yz, g_yz_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_z_xx[k] = -g_y_0_0_xx[k] - g_yz_0_0_xx[k] * ab_z + g_yz_0_0_xxz[k];

                g_yz_0_z_xy[k] = -g_y_0_0_xy[k] - g_yz_0_0_xy[k] * ab_z + g_yz_0_0_xyz[k];

                g_yz_0_z_xz[k] = -g_y_0_0_xz[k] - g_yz_0_0_xz[k] * ab_z + g_yz_0_0_xzz[k];

                g_yz_0_z_yy[k] = -g_y_0_0_yy[k] - g_yz_0_0_yy[k] * ab_z + g_yz_0_0_yyz[k];

                g_yz_0_z_yz[k] = -g_y_0_0_yz[k] - g_yz_0_0_yz[k] * ab_z + g_yz_0_0_yzz[k];

                g_yz_0_z_zz[k] = -g_y_0_0_zz[k] - g_yz_0_0_zz[k] * ab_z + g_yz_0_0_zzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_zz_0_x_xx = cbuffer.data(pd_geom_20_off + 90 * ccomps * dcomps);

            auto g_zz_0_x_xy = cbuffer.data(pd_geom_20_off + 91 * ccomps * dcomps);

            auto g_zz_0_x_xz = cbuffer.data(pd_geom_20_off + 92 * ccomps * dcomps);

            auto g_zz_0_x_yy = cbuffer.data(pd_geom_20_off + 93 * ccomps * dcomps);

            auto g_zz_0_x_yz = cbuffer.data(pd_geom_20_off + 94 * ccomps * dcomps);

            auto g_zz_0_x_zz = cbuffer.data(pd_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xx, g_zz_0_0_xxx, g_zz_0_0_xxy, g_zz_0_0_xxz, g_zz_0_0_xy, g_zz_0_0_xyy, g_zz_0_0_xyz, g_zz_0_0_xz, g_zz_0_0_xzz, g_zz_0_0_yy, g_zz_0_0_yz, g_zz_0_0_zz, g_zz_0_x_xx, g_zz_0_x_xy, g_zz_0_x_xz, g_zz_0_x_yy, g_zz_0_x_yz, g_zz_0_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_x_xx[k] = -g_zz_0_0_xx[k] * ab_x + g_zz_0_0_xxx[k];

                g_zz_0_x_xy[k] = -g_zz_0_0_xy[k] * ab_x + g_zz_0_0_xxy[k];

                g_zz_0_x_xz[k] = -g_zz_0_0_xz[k] * ab_x + g_zz_0_0_xxz[k];

                g_zz_0_x_yy[k] = -g_zz_0_0_yy[k] * ab_x + g_zz_0_0_xyy[k];

                g_zz_0_x_yz[k] = -g_zz_0_0_yz[k] * ab_x + g_zz_0_0_xyz[k];

                g_zz_0_x_zz[k] = -g_zz_0_0_zz[k] * ab_x + g_zz_0_0_xzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_zz_0_y_xx = cbuffer.data(pd_geom_20_off + 96 * ccomps * dcomps);

            auto g_zz_0_y_xy = cbuffer.data(pd_geom_20_off + 97 * ccomps * dcomps);

            auto g_zz_0_y_xz = cbuffer.data(pd_geom_20_off + 98 * ccomps * dcomps);

            auto g_zz_0_y_yy = cbuffer.data(pd_geom_20_off + 99 * ccomps * dcomps);

            auto g_zz_0_y_yz = cbuffer.data(pd_geom_20_off + 100 * ccomps * dcomps);

            auto g_zz_0_y_zz = cbuffer.data(pd_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xx, g_zz_0_0_xxy, g_zz_0_0_xy, g_zz_0_0_xyy, g_zz_0_0_xyz, g_zz_0_0_xz, g_zz_0_0_yy, g_zz_0_0_yyy, g_zz_0_0_yyz, g_zz_0_0_yz, g_zz_0_0_yzz, g_zz_0_0_zz, g_zz_0_y_xx, g_zz_0_y_xy, g_zz_0_y_xz, g_zz_0_y_yy, g_zz_0_y_yz, g_zz_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_y_xx[k] = -g_zz_0_0_xx[k] * ab_y + g_zz_0_0_xxy[k];

                g_zz_0_y_xy[k] = -g_zz_0_0_xy[k] * ab_y + g_zz_0_0_xyy[k];

                g_zz_0_y_xz[k] = -g_zz_0_0_xz[k] * ab_y + g_zz_0_0_xyz[k];

                g_zz_0_y_yy[k] = -g_zz_0_0_yy[k] * ab_y + g_zz_0_0_yyy[k];

                g_zz_0_y_yz[k] = -g_zz_0_0_yz[k] * ab_y + g_zz_0_0_yyz[k];

                g_zz_0_y_zz[k] = -g_zz_0_0_zz[k] * ab_y + g_zz_0_0_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_zz_0_z_xx = cbuffer.data(pd_geom_20_off + 102 * ccomps * dcomps);

            auto g_zz_0_z_xy = cbuffer.data(pd_geom_20_off + 103 * ccomps * dcomps);

            auto g_zz_0_z_xz = cbuffer.data(pd_geom_20_off + 104 * ccomps * dcomps);

            auto g_zz_0_z_yy = cbuffer.data(pd_geom_20_off + 105 * ccomps * dcomps);

            auto g_zz_0_z_yz = cbuffer.data(pd_geom_20_off + 106 * ccomps * dcomps);

            auto g_zz_0_z_zz = cbuffer.data(pd_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xx, g_z_0_0_xy, g_z_0_0_xz, g_z_0_0_yy, g_z_0_0_yz, g_z_0_0_zz, g_zz_0_0_xx, g_zz_0_0_xxz, g_zz_0_0_xy, g_zz_0_0_xyz, g_zz_0_0_xz, g_zz_0_0_xzz, g_zz_0_0_yy, g_zz_0_0_yyz, g_zz_0_0_yz, g_zz_0_0_yzz, g_zz_0_0_zz, g_zz_0_0_zzz, g_zz_0_z_xx, g_zz_0_z_xy, g_zz_0_z_xz, g_zz_0_z_yy, g_zz_0_z_yz, g_zz_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_z_xx[k] = -2.0 * g_z_0_0_xx[k] - g_zz_0_0_xx[k] * ab_z + g_zz_0_0_xxz[k];

                g_zz_0_z_xy[k] = -2.0 * g_z_0_0_xy[k] - g_zz_0_0_xy[k] * ab_z + g_zz_0_0_xyz[k];

                g_zz_0_z_xz[k] = -2.0 * g_z_0_0_xz[k] - g_zz_0_0_xz[k] * ab_z + g_zz_0_0_xzz[k];

                g_zz_0_z_yy[k] = -2.0 * g_z_0_0_yy[k] - g_zz_0_0_yy[k] * ab_z + g_zz_0_0_yyz[k];

                g_zz_0_z_yz[k] = -2.0 * g_z_0_0_yz[k] - g_zz_0_0_yz[k] * ab_z + g_zz_0_0_yzz[k];

                g_zz_0_z_zz[k] = -2.0 * g_z_0_0_zz[k] - g_zz_0_0_zz[k] * ab_z + g_zz_0_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

