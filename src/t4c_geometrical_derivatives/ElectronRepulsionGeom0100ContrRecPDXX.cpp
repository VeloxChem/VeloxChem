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

#include "ElectronRepulsionGeom0100ContrRecPDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_pdxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_pdxx,
                                            const size_t idx_sdxx,
                                            const size_t idx_geom_01_sdxx,
                                            const size_t idx_geom_01_sfxx,
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

            const auto sd_off = idx_sdxx + i * dcomps + j;

            auto g_0_xx = cbuffer.data(sd_off + 0 * ccomps * dcomps);

            auto g_0_xy = cbuffer.data(sd_off + 1 * ccomps * dcomps);

            auto g_0_xz = cbuffer.data(sd_off + 2 * ccomps * dcomps);

            auto g_0_yy = cbuffer.data(sd_off + 3 * ccomps * dcomps);

            auto g_0_yz = cbuffer.data(sd_off + 4 * ccomps * dcomps);

            auto g_0_zz = cbuffer.data(sd_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SDSS

            const auto sd_geom_01_off = idx_geom_01_sdxx + i * dcomps + j;

            auto g_0_x_0_xx = cbuffer.data(sd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xy = cbuffer.data(sd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xz = cbuffer.data(sd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_yy = cbuffer.data(sd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_yz = cbuffer.data(sd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_zz = cbuffer.data(sd_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_y_0_xx = cbuffer.data(sd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_y_0_xy = cbuffer.data(sd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_y_0_xz = cbuffer.data(sd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_y_0_yy = cbuffer.data(sd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_0_yz = cbuffer.data(sd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_0_zz = cbuffer.data(sd_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_z_0_xx = cbuffer.data(sd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_z_0_xy = cbuffer.data(sd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_z_0_xz = cbuffer.data(sd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_z_0_yy = cbuffer.data(sd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_z_0_yz = cbuffer.data(sd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_z_0_zz = cbuffer.data(sd_geom_01_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SFSS

            const auto sf_geom_01_off = idx_geom_01_sfxx + i * dcomps + j;

            auto g_0_x_0_xxx = cbuffer.data(sf_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxy = cbuffer.data(sf_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxz = cbuffer.data(sf_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xyy = cbuffer.data(sf_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xyz = cbuffer.data(sf_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xzz = cbuffer.data(sf_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_yyy = cbuffer.data(sf_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_yyz = cbuffer.data(sf_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_yzz = cbuffer.data(sf_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_zzz = cbuffer.data(sf_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_y_0_xxx = cbuffer.data(sf_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_y_0_xxy = cbuffer.data(sf_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_y_0_xxz = cbuffer.data(sf_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_y_0_xyy = cbuffer.data(sf_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_y_0_xyz = cbuffer.data(sf_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_y_0_xzz = cbuffer.data(sf_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_y_0_yyy = cbuffer.data(sf_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_y_0_yyz = cbuffer.data(sf_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_y_0_yzz = cbuffer.data(sf_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_0_zzz = cbuffer.data(sf_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_z_0_xxx = cbuffer.data(sf_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_z_0_xxy = cbuffer.data(sf_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_z_0_xxz = cbuffer.data(sf_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_z_0_xyy = cbuffer.data(sf_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_z_0_xyz = cbuffer.data(sf_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_z_0_xzz = cbuffer.data(sf_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_z_0_yyy = cbuffer.data(sf_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_z_0_yyz = cbuffer.data(sf_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_z_0_yzz = cbuffer.data(sf_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_z_0_zzz = cbuffer.data(sf_geom_01_off + 29 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pdxx

            const auto pd_geom_01_off = idx_geom_01_pdxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_x_xx = cbuffer.data(pd_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xy = cbuffer.data(pd_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xz = cbuffer.data(pd_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_yy = cbuffer.data(pd_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_yz = cbuffer.data(pd_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_zz = cbuffer.data(pd_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xxx, g_0_x_0_xxy, g_0_x_0_xxz, g_0_x_0_xy, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xz, g_0_x_0_xzz, g_0_x_0_yy, g_0_x_0_yz, g_0_x_0_zz, g_0_x_x_xx, g_0_x_x_xy, g_0_x_x_xz, g_0_x_x_yy, g_0_x_x_yz, g_0_x_x_zz, g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_x_xx[k] = g_0_xx[k] - g_0_x_0_xx[k] * ab_x + g_0_x_0_xxx[k];

                g_0_x_x_xy[k] = g_0_xy[k] - g_0_x_0_xy[k] * ab_x + g_0_x_0_xxy[k];

                g_0_x_x_xz[k] = g_0_xz[k] - g_0_x_0_xz[k] * ab_x + g_0_x_0_xxz[k];

                g_0_x_x_yy[k] = g_0_yy[k] - g_0_x_0_yy[k] * ab_x + g_0_x_0_xyy[k];

                g_0_x_x_yz[k] = g_0_yz[k] - g_0_x_0_yz[k] * ab_x + g_0_x_0_xyz[k];

                g_0_x_x_zz[k] = g_0_zz[k] - g_0_x_0_zz[k] * ab_x + g_0_x_0_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_y_xx = cbuffer.data(pd_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_y_xy = cbuffer.data(pd_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_y_xz = cbuffer.data(pd_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_y_yy = cbuffer.data(pd_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_y_yz = cbuffer.data(pd_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_y_zz = cbuffer.data(pd_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xxy, g_0_x_0_xy, g_0_x_0_xyy, g_0_x_0_xyz, g_0_x_0_xz, g_0_x_0_yy, g_0_x_0_yyy, g_0_x_0_yyz, g_0_x_0_yz, g_0_x_0_yzz, g_0_x_0_zz, g_0_x_y_xx, g_0_x_y_xy, g_0_x_y_xz, g_0_x_y_yy, g_0_x_y_yz, g_0_x_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_y_xx[k] = -g_0_x_0_xx[k] * ab_y + g_0_x_0_xxy[k];

                g_0_x_y_xy[k] = -g_0_x_0_xy[k] * ab_y + g_0_x_0_xyy[k];

                g_0_x_y_xz[k] = -g_0_x_0_xz[k] * ab_y + g_0_x_0_xyz[k];

                g_0_x_y_yy[k] = -g_0_x_0_yy[k] * ab_y + g_0_x_0_yyy[k];

                g_0_x_y_yz[k] = -g_0_x_0_yz[k] * ab_y + g_0_x_0_yyz[k];

                g_0_x_y_zz[k] = -g_0_x_0_zz[k] * ab_y + g_0_x_0_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_z_xx = cbuffer.data(pd_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_z_xy = cbuffer.data(pd_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_z_xz = cbuffer.data(pd_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_z_yy = cbuffer.data(pd_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_z_yz = cbuffer.data(pd_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_z_zz = cbuffer.data(pd_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xx, g_0_x_0_xxz, g_0_x_0_xy, g_0_x_0_xyz, g_0_x_0_xz, g_0_x_0_xzz, g_0_x_0_yy, g_0_x_0_yyz, g_0_x_0_yz, g_0_x_0_yzz, g_0_x_0_zz, g_0_x_0_zzz, g_0_x_z_xx, g_0_x_z_xy, g_0_x_z_xz, g_0_x_z_yy, g_0_x_z_yz, g_0_x_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_z_xx[k] = -g_0_x_0_xx[k] * ab_z + g_0_x_0_xxz[k];

                g_0_x_z_xy[k] = -g_0_x_0_xy[k] * ab_z + g_0_x_0_xyz[k];

                g_0_x_z_xz[k] = -g_0_x_0_xz[k] * ab_z + g_0_x_0_xzz[k];

                g_0_x_z_yy[k] = -g_0_x_0_yy[k] * ab_z + g_0_x_0_yyz[k];

                g_0_x_z_yz[k] = -g_0_x_0_yz[k] * ab_z + g_0_x_0_yzz[k];

                g_0_x_z_zz[k] = -g_0_x_0_zz[k] * ab_z + g_0_x_0_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_0_y_x_xx = cbuffer.data(pd_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_y_x_xy = cbuffer.data(pd_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_y_x_xz = cbuffer.data(pd_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_y_x_yy = cbuffer.data(pd_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_y_x_yz = cbuffer.data(pd_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_y_x_zz = cbuffer.data(pd_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xxx, g_0_y_0_xxy, g_0_y_0_xxz, g_0_y_0_xy, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xz, g_0_y_0_xzz, g_0_y_0_yy, g_0_y_0_yz, g_0_y_0_zz, g_0_y_x_xx, g_0_y_x_xy, g_0_y_x_xz, g_0_y_x_yy, g_0_y_x_yz, g_0_y_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_x_xx[k] = -g_0_y_0_xx[k] * ab_x + g_0_y_0_xxx[k];

                g_0_y_x_xy[k] = -g_0_y_0_xy[k] * ab_x + g_0_y_0_xxy[k];

                g_0_y_x_xz[k] = -g_0_y_0_xz[k] * ab_x + g_0_y_0_xxz[k];

                g_0_y_x_yy[k] = -g_0_y_0_yy[k] * ab_x + g_0_y_0_xyy[k];

                g_0_y_x_yz[k] = -g_0_y_0_yz[k] * ab_x + g_0_y_0_xyz[k];

                g_0_y_x_zz[k] = -g_0_y_0_zz[k] * ab_x + g_0_y_0_xzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_0_y_y_xx = cbuffer.data(pd_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_y_y_xy = cbuffer.data(pd_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_y_y_xz = cbuffer.data(pd_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_y_y_yy = cbuffer.data(pd_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_y_yz = cbuffer.data(pd_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_y_zz = cbuffer.data(pd_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xx, g_0_xy, g_0_xz, g_0_y_0_xx, g_0_y_0_xxy, g_0_y_0_xy, g_0_y_0_xyy, g_0_y_0_xyz, g_0_y_0_xz, g_0_y_0_yy, g_0_y_0_yyy, g_0_y_0_yyz, g_0_y_0_yz, g_0_y_0_yzz, g_0_y_0_zz, g_0_y_y_xx, g_0_y_y_xy, g_0_y_y_xz, g_0_y_y_yy, g_0_y_y_yz, g_0_y_y_zz, g_0_yy, g_0_yz, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_y_xx[k] = g_0_xx[k] - g_0_y_0_xx[k] * ab_y + g_0_y_0_xxy[k];

                g_0_y_y_xy[k] = g_0_xy[k] - g_0_y_0_xy[k] * ab_y + g_0_y_0_xyy[k];

                g_0_y_y_xz[k] = g_0_xz[k] - g_0_y_0_xz[k] * ab_y + g_0_y_0_xyz[k];

                g_0_y_y_yy[k] = g_0_yy[k] - g_0_y_0_yy[k] * ab_y + g_0_y_0_yyy[k];

                g_0_y_y_yz[k] = g_0_yz[k] - g_0_y_0_yz[k] * ab_y + g_0_y_0_yyz[k];

                g_0_y_y_zz[k] = g_0_zz[k] - g_0_y_0_zz[k] * ab_y + g_0_y_0_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_0_y_z_xx = cbuffer.data(pd_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_z_xy = cbuffer.data(pd_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_z_xz = cbuffer.data(pd_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_z_yy = cbuffer.data(pd_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_z_yz = cbuffer.data(pd_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_z_zz = cbuffer.data(pd_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xx, g_0_y_0_xxz, g_0_y_0_xy, g_0_y_0_xyz, g_0_y_0_xz, g_0_y_0_xzz, g_0_y_0_yy, g_0_y_0_yyz, g_0_y_0_yz, g_0_y_0_yzz, g_0_y_0_zz, g_0_y_0_zzz, g_0_y_z_xx, g_0_y_z_xy, g_0_y_z_xz, g_0_y_z_yy, g_0_y_z_yz, g_0_y_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_z_xx[k] = -g_0_y_0_xx[k] * ab_z + g_0_y_0_xxz[k];

                g_0_y_z_xy[k] = -g_0_y_0_xy[k] * ab_z + g_0_y_0_xyz[k];

                g_0_y_z_xz[k] = -g_0_y_0_xz[k] * ab_z + g_0_y_0_xzz[k];

                g_0_y_z_yy[k] = -g_0_y_0_yy[k] * ab_z + g_0_y_0_yyz[k];

                g_0_y_z_yz[k] = -g_0_y_0_yz[k] * ab_z + g_0_y_0_yzz[k];

                g_0_y_z_zz[k] = -g_0_y_0_zz[k] * ab_z + g_0_y_0_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_0_z_x_xx = cbuffer.data(pd_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_z_x_xy = cbuffer.data(pd_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_z_x_xz = cbuffer.data(pd_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_z_x_yy = cbuffer.data(pd_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_z_x_yz = cbuffer.data(pd_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_z_x_zz = cbuffer.data(pd_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xxx, g_0_z_0_xxy, g_0_z_0_xxz, g_0_z_0_xy, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xz, g_0_z_0_xzz, g_0_z_0_yy, g_0_z_0_yz, g_0_z_0_zz, g_0_z_x_xx, g_0_z_x_xy, g_0_z_x_xz, g_0_z_x_yy, g_0_z_x_yz, g_0_z_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_x_xx[k] = -g_0_z_0_xx[k] * ab_x + g_0_z_0_xxx[k];

                g_0_z_x_xy[k] = -g_0_z_0_xy[k] * ab_x + g_0_z_0_xxy[k];

                g_0_z_x_xz[k] = -g_0_z_0_xz[k] * ab_x + g_0_z_0_xxz[k];

                g_0_z_x_yy[k] = -g_0_z_0_yy[k] * ab_x + g_0_z_0_xyy[k];

                g_0_z_x_yz[k] = -g_0_z_0_yz[k] * ab_x + g_0_z_0_xyz[k];

                g_0_z_x_zz[k] = -g_0_z_0_zz[k] * ab_x + g_0_z_0_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_0_z_y_xx = cbuffer.data(pd_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_z_y_xy = cbuffer.data(pd_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_z_y_xz = cbuffer.data(pd_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_z_y_yy = cbuffer.data(pd_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_z_y_yz = cbuffer.data(pd_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_z_y_zz = cbuffer.data(pd_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xx, g_0_z_0_xxy, g_0_z_0_xy, g_0_z_0_xyy, g_0_z_0_xyz, g_0_z_0_xz, g_0_z_0_yy, g_0_z_0_yyy, g_0_z_0_yyz, g_0_z_0_yz, g_0_z_0_yzz, g_0_z_0_zz, g_0_z_y_xx, g_0_z_y_xy, g_0_z_y_xz, g_0_z_y_yy, g_0_z_y_yz, g_0_z_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_y_xx[k] = -g_0_z_0_xx[k] * ab_y + g_0_z_0_xxy[k];

                g_0_z_y_xy[k] = -g_0_z_0_xy[k] * ab_y + g_0_z_0_xyy[k];

                g_0_z_y_xz[k] = -g_0_z_0_xz[k] * ab_y + g_0_z_0_xyz[k];

                g_0_z_y_yy[k] = -g_0_z_0_yy[k] * ab_y + g_0_z_0_yyy[k];

                g_0_z_y_yz[k] = -g_0_z_0_yz[k] * ab_y + g_0_z_0_yyz[k];

                g_0_z_y_zz[k] = -g_0_z_0_zz[k] * ab_y + g_0_z_0_yzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_0_z_z_xx = cbuffer.data(pd_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_z_z_xy = cbuffer.data(pd_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_z_z_xz = cbuffer.data(pd_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_z_z_yy = cbuffer.data(pd_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_z_z_yz = cbuffer.data(pd_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_z_z_zz = cbuffer.data(pd_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xx, g_0_xy, g_0_xz, g_0_yy, g_0_yz, g_0_z_0_xx, g_0_z_0_xxz, g_0_z_0_xy, g_0_z_0_xyz, g_0_z_0_xz, g_0_z_0_xzz, g_0_z_0_yy, g_0_z_0_yyz, g_0_z_0_yz, g_0_z_0_yzz, g_0_z_0_zz, g_0_z_0_zzz, g_0_z_z_xx, g_0_z_z_xy, g_0_z_z_xz, g_0_z_z_yy, g_0_z_z_yz, g_0_z_z_zz, g_0_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_z_xx[k] = g_0_xx[k] - g_0_z_0_xx[k] * ab_z + g_0_z_0_xxz[k];

                g_0_z_z_xy[k] = g_0_xy[k] - g_0_z_0_xy[k] * ab_z + g_0_z_0_xyz[k];

                g_0_z_z_xz[k] = g_0_xz[k] - g_0_z_0_xz[k] * ab_z + g_0_z_0_xzz[k];

                g_0_z_z_yy[k] = g_0_yy[k] - g_0_z_0_yy[k] * ab_z + g_0_z_0_yyz[k];

                g_0_z_z_yz[k] = g_0_yz[k] - g_0_z_0_yz[k] * ab_z + g_0_z_0_yzz[k];

                g_0_z_z_zz[k] = g_0_zz[k] - g_0_z_0_zz[k] * ab_z + g_0_z_0_zzz[k];
            }
        }
    }
}

} // erirec namespace

