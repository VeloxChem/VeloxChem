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

#include "ElectronRepulsionGeom1000ContrRecDDXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_ddxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ddxx,
                                            const size_t idx_pdxx,
                                            const size_t idx_geom_10_pdxx,
                                            const size_t idx_geom_10_pfxx,
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
            /// Set up components of auxilary buffer : PDSS

            const auto pd_off = idx_pdxx + i * dcomps + j;

            auto g_x_xx = cbuffer.data(pd_off + 0 * ccomps * dcomps);

            auto g_x_xy = cbuffer.data(pd_off + 1 * ccomps * dcomps);

            auto g_x_xz = cbuffer.data(pd_off + 2 * ccomps * dcomps);

            auto g_x_yy = cbuffer.data(pd_off + 3 * ccomps * dcomps);

            auto g_x_yz = cbuffer.data(pd_off + 4 * ccomps * dcomps);

            auto g_x_zz = cbuffer.data(pd_off + 5 * ccomps * dcomps);

            auto g_y_xx = cbuffer.data(pd_off + 6 * ccomps * dcomps);

            auto g_y_xy = cbuffer.data(pd_off + 7 * ccomps * dcomps);

            auto g_y_xz = cbuffer.data(pd_off + 8 * ccomps * dcomps);

            auto g_y_yy = cbuffer.data(pd_off + 9 * ccomps * dcomps);

            auto g_y_yz = cbuffer.data(pd_off + 10 * ccomps * dcomps);

            auto g_y_zz = cbuffer.data(pd_off + 11 * ccomps * dcomps);

            auto g_z_xx = cbuffer.data(pd_off + 12 * ccomps * dcomps);

            auto g_z_xy = cbuffer.data(pd_off + 13 * ccomps * dcomps);

            auto g_z_xz = cbuffer.data(pd_off + 14 * ccomps * dcomps);

            auto g_z_yy = cbuffer.data(pd_off + 15 * ccomps * dcomps);

            auto g_z_yz = cbuffer.data(pd_off + 16 * ccomps * dcomps);

            auto g_z_zz = cbuffer.data(pd_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PDSS

            const auto pd_geom_10_off = idx_geom_10_pdxx + i * dcomps + j;

            auto g_x_0_x_xx = cbuffer.data(pd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xy = cbuffer.data(pd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xz = cbuffer.data(pd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_yy = cbuffer.data(pd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_yz = cbuffer.data(pd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_zz = cbuffer.data(pd_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_y_xx = cbuffer.data(pd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_y_xy = cbuffer.data(pd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_y_xz = cbuffer.data(pd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_y_yy = cbuffer.data(pd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_y_yz = cbuffer.data(pd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_y_zz = cbuffer.data(pd_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_z_xx = cbuffer.data(pd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_z_xy = cbuffer.data(pd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_z_xz = cbuffer.data(pd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_z_yy = cbuffer.data(pd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_z_yz = cbuffer.data(pd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_z_zz = cbuffer.data(pd_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_x_xx = cbuffer.data(pd_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_x_xy = cbuffer.data(pd_geom_10_off + 19 * ccomps * dcomps);

            auto g_y_0_x_xz = cbuffer.data(pd_geom_10_off + 20 * ccomps * dcomps);

            auto g_y_0_x_yy = cbuffer.data(pd_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_x_yz = cbuffer.data(pd_geom_10_off + 22 * ccomps * dcomps);

            auto g_y_0_x_zz = cbuffer.data(pd_geom_10_off + 23 * ccomps * dcomps);

            auto g_y_0_y_xx = cbuffer.data(pd_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_y_xy = cbuffer.data(pd_geom_10_off + 25 * ccomps * dcomps);

            auto g_y_0_y_xz = cbuffer.data(pd_geom_10_off + 26 * ccomps * dcomps);

            auto g_y_0_y_yy = cbuffer.data(pd_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_y_yz = cbuffer.data(pd_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_y_zz = cbuffer.data(pd_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_z_xx = cbuffer.data(pd_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_z_xy = cbuffer.data(pd_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_z_xz = cbuffer.data(pd_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_z_yy = cbuffer.data(pd_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_z_yz = cbuffer.data(pd_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_z_zz = cbuffer.data(pd_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_x_xx = cbuffer.data(pd_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_x_xy = cbuffer.data(pd_geom_10_off + 37 * ccomps * dcomps);

            auto g_z_0_x_xz = cbuffer.data(pd_geom_10_off + 38 * ccomps * dcomps);

            auto g_z_0_x_yy = cbuffer.data(pd_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_x_yz = cbuffer.data(pd_geom_10_off + 40 * ccomps * dcomps);

            auto g_z_0_x_zz = cbuffer.data(pd_geom_10_off + 41 * ccomps * dcomps);

            auto g_z_0_y_xx = cbuffer.data(pd_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_y_xy = cbuffer.data(pd_geom_10_off + 43 * ccomps * dcomps);

            auto g_z_0_y_xz = cbuffer.data(pd_geom_10_off + 44 * ccomps * dcomps);

            auto g_z_0_y_yy = cbuffer.data(pd_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_y_yz = cbuffer.data(pd_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_y_zz = cbuffer.data(pd_geom_10_off + 47 * ccomps * dcomps);

            auto g_z_0_z_xx = cbuffer.data(pd_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_z_xy = cbuffer.data(pd_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_z_xz = cbuffer.data(pd_geom_10_off + 50 * ccomps * dcomps);

            auto g_z_0_z_yy = cbuffer.data(pd_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_z_yz = cbuffer.data(pd_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_z_zz = cbuffer.data(pd_geom_10_off + 53 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PFSS

            const auto pf_geom_10_off = idx_geom_10_pfxx + i * dcomps + j;

            auto g_x_0_x_xxx = cbuffer.data(pf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxy = cbuffer.data(pf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxz = cbuffer.data(pf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xyy = cbuffer.data(pf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xyz = cbuffer.data(pf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xzz = cbuffer.data(pf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_yyy = cbuffer.data(pf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_yyz = cbuffer.data(pf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_yzz = cbuffer.data(pf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_zzz = cbuffer.data(pf_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_y_xxy = cbuffer.data(pf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_y_xyy = cbuffer.data(pf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_y_xyz = cbuffer.data(pf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_y_yyy = cbuffer.data(pf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_y_yyz = cbuffer.data(pf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_y_yzz = cbuffer.data(pf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_z_xxy = cbuffer.data(pf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_z_xxz = cbuffer.data(pf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_z_xyy = cbuffer.data(pf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_z_xyz = cbuffer.data(pf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_z_xzz = cbuffer.data(pf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_z_yyy = cbuffer.data(pf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_z_yyz = cbuffer.data(pf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_z_yzz = cbuffer.data(pf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_z_zzz = cbuffer.data(pf_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_x_xxx = cbuffer.data(pf_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_x_xxy = cbuffer.data(pf_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_x_xxz = cbuffer.data(pf_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_x_xyy = cbuffer.data(pf_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_x_xyz = cbuffer.data(pf_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_x_xzz = cbuffer.data(pf_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_y_xxx = cbuffer.data(pf_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_y_xxy = cbuffer.data(pf_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_y_xxz = cbuffer.data(pf_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_y_xyy = cbuffer.data(pf_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_y_xyz = cbuffer.data(pf_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_y_xzz = cbuffer.data(pf_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_y_yyy = cbuffer.data(pf_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_y_yyz = cbuffer.data(pf_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_y_yzz = cbuffer.data(pf_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_y_zzz = cbuffer.data(pf_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_z_xxx = cbuffer.data(pf_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_z_xxy = cbuffer.data(pf_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_z_xxz = cbuffer.data(pf_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_z_xyy = cbuffer.data(pf_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_z_xyz = cbuffer.data(pf_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_z_xzz = cbuffer.data(pf_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_z_yyz = cbuffer.data(pf_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_z_yzz = cbuffer.data(pf_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_z_zzz = cbuffer.data(pf_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_x_xxx = cbuffer.data(pf_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_x_xxy = cbuffer.data(pf_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_x_xxz = cbuffer.data(pf_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_x_xyy = cbuffer.data(pf_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_x_xyz = cbuffer.data(pf_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_x_xzz = cbuffer.data(pf_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_y_xxx = cbuffer.data(pf_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_y_xxy = cbuffer.data(pf_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_y_xxz = cbuffer.data(pf_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_y_xyy = cbuffer.data(pf_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_y_xyz = cbuffer.data(pf_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_y_xzz = cbuffer.data(pf_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_y_yyy = cbuffer.data(pf_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_y_yyz = cbuffer.data(pf_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_y_yzz = cbuffer.data(pf_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_z_xxx = cbuffer.data(pf_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_z_xxy = cbuffer.data(pf_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_z_xxz = cbuffer.data(pf_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_z_xyy = cbuffer.data(pf_geom_10_off + 83 * ccomps * dcomps);

            auto g_z_0_z_xyz = cbuffer.data(pf_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_z_xzz = cbuffer.data(pf_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_z_yyy = cbuffer.data(pf_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_z_yyz = cbuffer.data(pf_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_z_yzz = cbuffer.data(pf_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_z_zzz = cbuffer.data(pf_geom_10_off + 89 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ddxx

            const auto dd_geom_10_off = idx_geom_10_ddxx + i * dcomps + j;

            /// Set up 0-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_xx = cbuffer.data(dd_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xy = cbuffer.data(dd_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xz = cbuffer.data(dd_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_yy = cbuffer.data(dd_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_yz = cbuffer.data(dd_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_zz = cbuffer.data(dd_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xxx, g_x_0_x_xxy, g_x_0_x_xxz, g_x_0_x_xy, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_xzz, g_x_0_x_yy, g_x_0_x_yz, g_x_0_x_zz, g_x_0_xx_xx, g_x_0_xx_xy, g_x_0_xx_xz, g_x_0_xx_yy, g_x_0_xx_yz, g_x_0_xx_zz, g_x_xx, g_x_xy, g_x_xz, g_x_yy, g_x_yz, g_x_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_xx[k] = -g_x_xx[k] - g_x_0_x_xx[k] * ab_x + g_x_0_x_xxx[k];

                g_x_0_xx_xy[k] = -g_x_xy[k] - g_x_0_x_xy[k] * ab_x + g_x_0_x_xxy[k];

                g_x_0_xx_xz[k] = -g_x_xz[k] - g_x_0_x_xz[k] * ab_x + g_x_0_x_xxz[k];

                g_x_0_xx_yy[k] = -g_x_yy[k] - g_x_0_x_yy[k] * ab_x + g_x_0_x_xyy[k];

                g_x_0_xx_yz[k] = -g_x_yz[k] - g_x_0_x_yz[k] * ab_x + g_x_0_x_xyz[k];

                g_x_0_xx_zz[k] = -g_x_zz[k] - g_x_0_x_zz[k] * ab_x + g_x_0_x_xzz[k];
            }

            /// Set up 6-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_xx = cbuffer.data(dd_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xy_xy = cbuffer.data(dd_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xy_xz = cbuffer.data(dd_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xy_yy = cbuffer.data(dd_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xy_yz = cbuffer.data(dd_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xy_zz = cbuffer.data(dd_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xxy, g_x_0_x_xy, g_x_0_x_xyy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_yy, g_x_0_x_yyy, g_x_0_x_yyz, g_x_0_x_yz, g_x_0_x_yzz, g_x_0_x_zz, g_x_0_xy_xx, g_x_0_xy_xy, g_x_0_xy_xz, g_x_0_xy_yy, g_x_0_xy_yz, g_x_0_xy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_xx[k] = -g_x_0_x_xx[k] * ab_y + g_x_0_x_xxy[k];

                g_x_0_xy_xy[k] = -g_x_0_x_xy[k] * ab_y + g_x_0_x_xyy[k];

                g_x_0_xy_xz[k] = -g_x_0_x_xz[k] * ab_y + g_x_0_x_xyz[k];

                g_x_0_xy_yy[k] = -g_x_0_x_yy[k] * ab_y + g_x_0_x_yyy[k];

                g_x_0_xy_yz[k] = -g_x_0_x_yz[k] * ab_y + g_x_0_x_yyz[k];

                g_x_0_xy_zz[k] = -g_x_0_x_zz[k] * ab_y + g_x_0_x_yzz[k];
            }

            /// Set up 12-18 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_xx = cbuffer.data(dd_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xz_xy = cbuffer.data(dd_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xz_xz = cbuffer.data(dd_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xz_yy = cbuffer.data(dd_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xz_yz = cbuffer.data(dd_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xz_zz = cbuffer.data(dd_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xx, g_x_0_x_xxz, g_x_0_x_xy, g_x_0_x_xyz, g_x_0_x_xz, g_x_0_x_xzz, g_x_0_x_yy, g_x_0_x_yyz, g_x_0_x_yz, g_x_0_x_yzz, g_x_0_x_zz, g_x_0_x_zzz, g_x_0_xz_xx, g_x_0_xz_xy, g_x_0_xz_xz, g_x_0_xz_yy, g_x_0_xz_yz, g_x_0_xz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_xx[k] = -g_x_0_x_xx[k] * ab_z + g_x_0_x_xxz[k];

                g_x_0_xz_xy[k] = -g_x_0_x_xy[k] * ab_z + g_x_0_x_xyz[k];

                g_x_0_xz_xz[k] = -g_x_0_x_xz[k] * ab_z + g_x_0_x_xzz[k];

                g_x_0_xz_yy[k] = -g_x_0_x_yy[k] * ab_z + g_x_0_x_yyz[k];

                g_x_0_xz_yz[k] = -g_x_0_x_yz[k] * ab_z + g_x_0_x_yzz[k];

                g_x_0_xz_zz[k] = -g_x_0_x_zz[k] * ab_z + g_x_0_x_zzz[k];
            }

            /// Set up 18-24 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_xx = cbuffer.data(dd_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_yy_xy = cbuffer.data(dd_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_yy_xz = cbuffer.data(dd_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_yy_yy = cbuffer.data(dd_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_yy_yz = cbuffer.data(dd_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_yy_zz = cbuffer.data(dd_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xx, g_x_0_y_xxy, g_x_0_y_xy, g_x_0_y_xyy, g_x_0_y_xyz, g_x_0_y_xz, g_x_0_y_yy, g_x_0_y_yyy, g_x_0_y_yyz, g_x_0_y_yz, g_x_0_y_yzz, g_x_0_y_zz, g_x_0_yy_xx, g_x_0_yy_xy, g_x_0_yy_xz, g_x_0_yy_yy, g_x_0_yy_yz, g_x_0_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_xx[k] = -g_x_0_y_xx[k] * ab_y + g_x_0_y_xxy[k];

                g_x_0_yy_xy[k] = -g_x_0_y_xy[k] * ab_y + g_x_0_y_xyy[k];

                g_x_0_yy_xz[k] = -g_x_0_y_xz[k] * ab_y + g_x_0_y_xyz[k];

                g_x_0_yy_yy[k] = -g_x_0_y_yy[k] * ab_y + g_x_0_y_yyy[k];

                g_x_0_yy_yz[k] = -g_x_0_y_yz[k] * ab_y + g_x_0_y_yyz[k];

                g_x_0_yy_zz[k] = -g_x_0_y_zz[k] * ab_y + g_x_0_y_yzz[k];
            }

            /// Set up 24-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_xx = cbuffer.data(dd_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_yz_xy = cbuffer.data(dd_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_yz_xz = cbuffer.data(dd_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_yz_yy = cbuffer.data(dd_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_yz_yz = cbuffer.data(dd_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_yz_zz = cbuffer.data(dd_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yz_xx, g_x_0_yz_xy, g_x_0_yz_xz, g_x_0_yz_yy, g_x_0_yz_yz, g_x_0_yz_zz, g_x_0_z_xx, g_x_0_z_xxy, g_x_0_z_xy, g_x_0_z_xyy, g_x_0_z_xyz, g_x_0_z_xz, g_x_0_z_yy, g_x_0_z_yyy, g_x_0_z_yyz, g_x_0_z_yz, g_x_0_z_yzz, g_x_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_xx[k] = -g_x_0_z_xx[k] * ab_y + g_x_0_z_xxy[k];

                g_x_0_yz_xy[k] = -g_x_0_z_xy[k] * ab_y + g_x_0_z_xyy[k];

                g_x_0_yz_xz[k] = -g_x_0_z_xz[k] * ab_y + g_x_0_z_xyz[k];

                g_x_0_yz_yy[k] = -g_x_0_z_yy[k] * ab_y + g_x_0_z_yyy[k];

                g_x_0_yz_yz[k] = -g_x_0_z_yz[k] * ab_y + g_x_0_z_yyz[k];

                g_x_0_yz_zz[k] = -g_x_0_z_zz[k] * ab_y + g_x_0_z_yzz[k];
            }

            /// Set up 30-36 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_xx = cbuffer.data(dd_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_zz_xy = cbuffer.data(dd_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_zz_xz = cbuffer.data(dd_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_zz_yy = cbuffer.data(dd_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_zz_yz = cbuffer.data(dd_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_zz_zz = cbuffer.data(dd_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xx, g_x_0_z_xxz, g_x_0_z_xy, g_x_0_z_xyz, g_x_0_z_xz, g_x_0_z_xzz, g_x_0_z_yy, g_x_0_z_yyz, g_x_0_z_yz, g_x_0_z_yzz, g_x_0_z_zz, g_x_0_z_zzz, g_x_0_zz_xx, g_x_0_zz_xy, g_x_0_zz_xz, g_x_0_zz_yy, g_x_0_zz_yz, g_x_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_xx[k] = -g_x_0_z_xx[k] * ab_z + g_x_0_z_xxz[k];

                g_x_0_zz_xy[k] = -g_x_0_z_xy[k] * ab_z + g_x_0_z_xyz[k];

                g_x_0_zz_xz[k] = -g_x_0_z_xz[k] * ab_z + g_x_0_z_xzz[k];

                g_x_0_zz_yy[k] = -g_x_0_z_yy[k] * ab_z + g_x_0_z_yyz[k];

                g_x_0_zz_yz[k] = -g_x_0_z_yz[k] * ab_z + g_x_0_z_yzz[k];

                g_x_0_zz_zz[k] = -g_x_0_z_zz[k] * ab_z + g_x_0_z_zzz[k];
            }

            /// Set up 36-42 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_xx = cbuffer.data(dd_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_xx_xy = cbuffer.data(dd_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_xx_xz = cbuffer.data(dd_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_xx_yy = cbuffer.data(dd_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_xx_yz = cbuffer.data(dd_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_xx_zz = cbuffer.data(dd_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xx, g_y_0_x_xxx, g_y_0_x_xxy, g_y_0_x_xxz, g_y_0_x_xy, g_y_0_x_xyy, g_y_0_x_xyz, g_y_0_x_xz, g_y_0_x_xzz, g_y_0_x_yy, g_y_0_x_yz, g_y_0_x_zz, g_y_0_xx_xx, g_y_0_xx_xy, g_y_0_xx_xz, g_y_0_xx_yy, g_y_0_xx_yz, g_y_0_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_xx[k] = -g_y_0_x_xx[k] * ab_x + g_y_0_x_xxx[k];

                g_y_0_xx_xy[k] = -g_y_0_x_xy[k] * ab_x + g_y_0_x_xxy[k];

                g_y_0_xx_xz[k] = -g_y_0_x_xz[k] * ab_x + g_y_0_x_xxz[k];

                g_y_0_xx_yy[k] = -g_y_0_x_yy[k] * ab_x + g_y_0_x_xyy[k];

                g_y_0_xx_yz[k] = -g_y_0_x_yz[k] * ab_x + g_y_0_x_xyz[k];

                g_y_0_xx_zz[k] = -g_y_0_x_zz[k] * ab_x + g_y_0_x_xzz[k];
            }

            /// Set up 42-48 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_xx = cbuffer.data(dd_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_xy_xy = cbuffer.data(dd_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_xy_xz = cbuffer.data(dd_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_xy_yy = cbuffer.data(dd_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_xy_yz = cbuffer.data(dd_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_xy_zz = cbuffer.data(dd_geom_10_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xx, g_y_0_xy_xy, g_y_0_xy_xz, g_y_0_xy_yy, g_y_0_xy_yz, g_y_0_xy_zz, g_y_0_y_xx, g_y_0_y_xxx, g_y_0_y_xxy, g_y_0_y_xxz, g_y_0_y_xy, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_xzz, g_y_0_y_yy, g_y_0_y_yz, g_y_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_xx[k] = -g_y_0_y_xx[k] * ab_x + g_y_0_y_xxx[k];

                g_y_0_xy_xy[k] = -g_y_0_y_xy[k] * ab_x + g_y_0_y_xxy[k];

                g_y_0_xy_xz[k] = -g_y_0_y_xz[k] * ab_x + g_y_0_y_xxz[k];

                g_y_0_xy_yy[k] = -g_y_0_y_yy[k] * ab_x + g_y_0_y_xyy[k];

                g_y_0_xy_yz[k] = -g_y_0_y_yz[k] * ab_x + g_y_0_y_xyz[k];

                g_y_0_xy_zz[k] = -g_y_0_y_zz[k] * ab_x + g_y_0_y_xzz[k];
            }

            /// Set up 48-54 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_xx = cbuffer.data(dd_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_xz_xy = cbuffer.data(dd_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_xz_xz = cbuffer.data(dd_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_xz_yy = cbuffer.data(dd_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_xz_yz = cbuffer.data(dd_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_xz_zz = cbuffer.data(dd_geom_10_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xz_xx, g_y_0_xz_xy, g_y_0_xz_xz, g_y_0_xz_yy, g_y_0_xz_yz, g_y_0_xz_zz, g_y_0_z_xx, g_y_0_z_xxx, g_y_0_z_xxy, g_y_0_z_xxz, g_y_0_z_xy, g_y_0_z_xyy, g_y_0_z_xyz, g_y_0_z_xz, g_y_0_z_xzz, g_y_0_z_yy, g_y_0_z_yz, g_y_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_xx[k] = -g_y_0_z_xx[k] * ab_x + g_y_0_z_xxx[k];

                g_y_0_xz_xy[k] = -g_y_0_z_xy[k] * ab_x + g_y_0_z_xxy[k];

                g_y_0_xz_xz[k] = -g_y_0_z_xz[k] * ab_x + g_y_0_z_xxz[k];

                g_y_0_xz_yy[k] = -g_y_0_z_yy[k] * ab_x + g_y_0_z_xyy[k];

                g_y_0_xz_yz[k] = -g_y_0_z_yz[k] * ab_x + g_y_0_z_xyz[k];

                g_y_0_xz_zz[k] = -g_y_0_z_zz[k] * ab_x + g_y_0_z_xzz[k];
            }

            /// Set up 54-60 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_xx = cbuffer.data(dd_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_yy_xy = cbuffer.data(dd_geom_10_off + 55 * ccomps * dcomps);

            auto g_y_0_yy_xz = cbuffer.data(dd_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_yy_yy = cbuffer.data(dd_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_yy_yz = cbuffer.data(dd_geom_10_off + 58 * ccomps * dcomps);

            auto g_y_0_yy_zz = cbuffer.data(dd_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xxy, g_y_0_y_xy, g_y_0_y_xyy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_yy, g_y_0_y_yyy, g_y_0_y_yyz, g_y_0_y_yz, g_y_0_y_yzz, g_y_0_y_zz, g_y_0_yy_xx, g_y_0_yy_xy, g_y_0_yy_xz, g_y_0_yy_yy, g_y_0_yy_yz, g_y_0_yy_zz, g_y_xx, g_y_xy, g_y_xz, g_y_yy, g_y_yz, g_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_xx[k] = -g_y_xx[k] - g_y_0_y_xx[k] * ab_y + g_y_0_y_xxy[k];

                g_y_0_yy_xy[k] = -g_y_xy[k] - g_y_0_y_xy[k] * ab_y + g_y_0_y_xyy[k];

                g_y_0_yy_xz[k] = -g_y_xz[k] - g_y_0_y_xz[k] * ab_y + g_y_0_y_xyz[k];

                g_y_0_yy_yy[k] = -g_y_yy[k] - g_y_0_y_yy[k] * ab_y + g_y_0_y_yyy[k];

                g_y_0_yy_yz[k] = -g_y_yz[k] - g_y_0_y_yz[k] * ab_y + g_y_0_y_yyz[k];

                g_y_0_yy_zz[k] = -g_y_zz[k] - g_y_0_y_zz[k] * ab_y + g_y_0_y_yzz[k];
            }

            /// Set up 60-66 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_xx = cbuffer.data(dd_geom_10_off + 60 * ccomps * dcomps);

            auto g_y_0_yz_xy = cbuffer.data(dd_geom_10_off + 61 * ccomps * dcomps);

            auto g_y_0_yz_xz = cbuffer.data(dd_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_yz_yy = cbuffer.data(dd_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_yz_yz = cbuffer.data(dd_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_yz_zz = cbuffer.data(dd_geom_10_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xx, g_y_0_y_xxz, g_y_0_y_xy, g_y_0_y_xyz, g_y_0_y_xz, g_y_0_y_xzz, g_y_0_y_yy, g_y_0_y_yyz, g_y_0_y_yz, g_y_0_y_yzz, g_y_0_y_zz, g_y_0_y_zzz, g_y_0_yz_xx, g_y_0_yz_xy, g_y_0_yz_xz, g_y_0_yz_yy, g_y_0_yz_yz, g_y_0_yz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_xx[k] = -g_y_0_y_xx[k] * ab_z + g_y_0_y_xxz[k];

                g_y_0_yz_xy[k] = -g_y_0_y_xy[k] * ab_z + g_y_0_y_xyz[k];

                g_y_0_yz_xz[k] = -g_y_0_y_xz[k] * ab_z + g_y_0_y_xzz[k];

                g_y_0_yz_yy[k] = -g_y_0_y_yy[k] * ab_z + g_y_0_y_yyz[k];

                g_y_0_yz_yz[k] = -g_y_0_y_yz[k] * ab_z + g_y_0_y_yzz[k];

                g_y_0_yz_zz[k] = -g_y_0_y_zz[k] * ab_z + g_y_0_y_zzz[k];
            }

            /// Set up 66-72 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_xx = cbuffer.data(dd_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_zz_xy = cbuffer.data(dd_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_zz_xz = cbuffer.data(dd_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_zz_yy = cbuffer.data(dd_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_zz_yz = cbuffer.data(dd_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_zz_zz = cbuffer.data(dd_geom_10_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xx, g_y_0_z_xxz, g_y_0_z_xy, g_y_0_z_xyz, g_y_0_z_xz, g_y_0_z_xzz, g_y_0_z_yy, g_y_0_z_yyz, g_y_0_z_yz, g_y_0_z_yzz, g_y_0_z_zz, g_y_0_z_zzz, g_y_0_zz_xx, g_y_0_zz_xy, g_y_0_zz_xz, g_y_0_zz_yy, g_y_0_zz_yz, g_y_0_zz_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_xx[k] = -g_y_0_z_xx[k] * ab_z + g_y_0_z_xxz[k];

                g_y_0_zz_xy[k] = -g_y_0_z_xy[k] * ab_z + g_y_0_z_xyz[k];

                g_y_0_zz_xz[k] = -g_y_0_z_xz[k] * ab_z + g_y_0_z_xzz[k];

                g_y_0_zz_yy[k] = -g_y_0_z_yy[k] * ab_z + g_y_0_z_yyz[k];

                g_y_0_zz_yz[k] = -g_y_0_z_yz[k] * ab_z + g_y_0_z_yzz[k];

                g_y_0_zz_zz[k] = -g_y_0_z_zz[k] * ab_z + g_y_0_z_zzz[k];
            }

            /// Set up 72-78 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_xx = cbuffer.data(dd_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_xx_xy = cbuffer.data(dd_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_xx_xz = cbuffer.data(dd_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_xx_yy = cbuffer.data(dd_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_xx_yz = cbuffer.data(dd_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_xx_zz = cbuffer.data(dd_geom_10_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xx, g_z_0_x_xxx, g_z_0_x_xxy, g_z_0_x_xxz, g_z_0_x_xy, g_z_0_x_xyy, g_z_0_x_xyz, g_z_0_x_xz, g_z_0_x_xzz, g_z_0_x_yy, g_z_0_x_yz, g_z_0_x_zz, g_z_0_xx_xx, g_z_0_xx_xy, g_z_0_xx_xz, g_z_0_xx_yy, g_z_0_xx_yz, g_z_0_xx_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_xx[k] = -g_z_0_x_xx[k] * ab_x + g_z_0_x_xxx[k];

                g_z_0_xx_xy[k] = -g_z_0_x_xy[k] * ab_x + g_z_0_x_xxy[k];

                g_z_0_xx_xz[k] = -g_z_0_x_xz[k] * ab_x + g_z_0_x_xxz[k];

                g_z_0_xx_yy[k] = -g_z_0_x_yy[k] * ab_x + g_z_0_x_xyy[k];

                g_z_0_xx_yz[k] = -g_z_0_x_yz[k] * ab_x + g_z_0_x_xyz[k];

                g_z_0_xx_zz[k] = -g_z_0_x_zz[k] * ab_x + g_z_0_x_xzz[k];
            }

            /// Set up 78-84 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_xx = cbuffer.data(dd_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_xy_xy = cbuffer.data(dd_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_xy_xz = cbuffer.data(dd_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_xy_yy = cbuffer.data(dd_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_xy_yz = cbuffer.data(dd_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_xy_zz = cbuffer.data(dd_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xy_xx, g_z_0_xy_xy, g_z_0_xy_xz, g_z_0_xy_yy, g_z_0_xy_yz, g_z_0_xy_zz, g_z_0_y_xx, g_z_0_y_xxx, g_z_0_y_xxy, g_z_0_y_xxz, g_z_0_y_xy, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xz, g_z_0_y_xzz, g_z_0_y_yy, g_z_0_y_yz, g_z_0_y_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_xx[k] = -g_z_0_y_xx[k] * ab_x + g_z_0_y_xxx[k];

                g_z_0_xy_xy[k] = -g_z_0_y_xy[k] * ab_x + g_z_0_y_xxy[k];

                g_z_0_xy_xz[k] = -g_z_0_y_xz[k] * ab_x + g_z_0_y_xxz[k];

                g_z_0_xy_yy[k] = -g_z_0_y_yy[k] * ab_x + g_z_0_y_xyy[k];

                g_z_0_xy_yz[k] = -g_z_0_y_yz[k] * ab_x + g_z_0_y_xyz[k];

                g_z_0_xy_zz[k] = -g_z_0_y_zz[k] * ab_x + g_z_0_y_xzz[k];
            }

            /// Set up 84-90 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_xx = cbuffer.data(dd_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_xz_xy = cbuffer.data(dd_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_xz_xz = cbuffer.data(dd_geom_10_off + 86 * ccomps * dcomps);

            auto g_z_0_xz_yy = cbuffer.data(dd_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_xz_yz = cbuffer.data(dd_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_xz_zz = cbuffer.data(dd_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xx, g_z_0_xz_xy, g_z_0_xz_xz, g_z_0_xz_yy, g_z_0_xz_yz, g_z_0_xz_zz, g_z_0_z_xx, g_z_0_z_xxx, g_z_0_z_xxy, g_z_0_z_xxz, g_z_0_z_xy, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_xzz, g_z_0_z_yy, g_z_0_z_yz, g_z_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_xx[k] = -g_z_0_z_xx[k] * ab_x + g_z_0_z_xxx[k];

                g_z_0_xz_xy[k] = -g_z_0_z_xy[k] * ab_x + g_z_0_z_xxy[k];

                g_z_0_xz_xz[k] = -g_z_0_z_xz[k] * ab_x + g_z_0_z_xxz[k];

                g_z_0_xz_yy[k] = -g_z_0_z_yy[k] * ab_x + g_z_0_z_xyy[k];

                g_z_0_xz_yz[k] = -g_z_0_z_yz[k] * ab_x + g_z_0_z_xyz[k];

                g_z_0_xz_zz[k] = -g_z_0_z_zz[k] * ab_x + g_z_0_z_xzz[k];
            }

            /// Set up 90-96 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_xx = cbuffer.data(dd_geom_10_off + 90 * ccomps * dcomps);

            auto g_z_0_yy_xy = cbuffer.data(dd_geom_10_off + 91 * ccomps * dcomps);

            auto g_z_0_yy_xz = cbuffer.data(dd_geom_10_off + 92 * ccomps * dcomps);

            auto g_z_0_yy_yy = cbuffer.data(dd_geom_10_off + 93 * ccomps * dcomps);

            auto g_z_0_yy_yz = cbuffer.data(dd_geom_10_off + 94 * ccomps * dcomps);

            auto g_z_0_yy_zz = cbuffer.data(dd_geom_10_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xx, g_z_0_y_xxy, g_z_0_y_xy, g_z_0_y_xyy, g_z_0_y_xyz, g_z_0_y_xz, g_z_0_y_yy, g_z_0_y_yyy, g_z_0_y_yyz, g_z_0_y_yz, g_z_0_y_yzz, g_z_0_y_zz, g_z_0_yy_xx, g_z_0_yy_xy, g_z_0_yy_xz, g_z_0_yy_yy, g_z_0_yy_yz, g_z_0_yy_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_xx[k] = -g_z_0_y_xx[k] * ab_y + g_z_0_y_xxy[k];

                g_z_0_yy_xy[k] = -g_z_0_y_xy[k] * ab_y + g_z_0_y_xyy[k];

                g_z_0_yy_xz[k] = -g_z_0_y_xz[k] * ab_y + g_z_0_y_xyz[k];

                g_z_0_yy_yy[k] = -g_z_0_y_yy[k] * ab_y + g_z_0_y_yyy[k];

                g_z_0_yy_yz[k] = -g_z_0_y_yz[k] * ab_y + g_z_0_y_yyz[k];

                g_z_0_yy_zz[k] = -g_z_0_y_zz[k] * ab_y + g_z_0_y_yzz[k];
            }

            /// Set up 96-102 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_xx = cbuffer.data(dd_geom_10_off + 96 * ccomps * dcomps);

            auto g_z_0_yz_xy = cbuffer.data(dd_geom_10_off + 97 * ccomps * dcomps);

            auto g_z_0_yz_xz = cbuffer.data(dd_geom_10_off + 98 * ccomps * dcomps);

            auto g_z_0_yz_yy = cbuffer.data(dd_geom_10_off + 99 * ccomps * dcomps);

            auto g_z_0_yz_yz = cbuffer.data(dd_geom_10_off + 100 * ccomps * dcomps);

            auto g_z_0_yz_zz = cbuffer.data(dd_geom_10_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xx, g_z_0_yz_xy, g_z_0_yz_xz, g_z_0_yz_yy, g_z_0_yz_yz, g_z_0_yz_zz, g_z_0_z_xx, g_z_0_z_xxy, g_z_0_z_xy, g_z_0_z_xyy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_yy, g_z_0_z_yyy, g_z_0_z_yyz, g_z_0_z_yz, g_z_0_z_yzz, g_z_0_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_xx[k] = -g_z_0_z_xx[k] * ab_y + g_z_0_z_xxy[k];

                g_z_0_yz_xy[k] = -g_z_0_z_xy[k] * ab_y + g_z_0_z_xyy[k];

                g_z_0_yz_xz[k] = -g_z_0_z_xz[k] * ab_y + g_z_0_z_xyz[k];

                g_z_0_yz_yy[k] = -g_z_0_z_yy[k] * ab_y + g_z_0_z_yyy[k];

                g_z_0_yz_yz[k] = -g_z_0_z_yz[k] * ab_y + g_z_0_z_yyz[k];

                g_z_0_yz_zz[k] = -g_z_0_z_zz[k] * ab_y + g_z_0_z_yzz[k];
            }

            /// Set up 102-108 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_xx = cbuffer.data(dd_geom_10_off + 102 * ccomps * dcomps);

            auto g_z_0_zz_xy = cbuffer.data(dd_geom_10_off + 103 * ccomps * dcomps);

            auto g_z_0_zz_xz = cbuffer.data(dd_geom_10_off + 104 * ccomps * dcomps);

            auto g_z_0_zz_yy = cbuffer.data(dd_geom_10_off + 105 * ccomps * dcomps);

            auto g_z_0_zz_yz = cbuffer.data(dd_geom_10_off + 106 * ccomps * dcomps);

            auto g_z_0_zz_zz = cbuffer.data(dd_geom_10_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xx, g_z_0_z_xxz, g_z_0_z_xy, g_z_0_z_xyz, g_z_0_z_xz, g_z_0_z_xzz, g_z_0_z_yy, g_z_0_z_yyz, g_z_0_z_yz, g_z_0_z_yzz, g_z_0_z_zz, g_z_0_z_zzz, g_z_0_zz_xx, g_z_0_zz_xy, g_z_0_zz_xz, g_z_0_zz_yy, g_z_0_zz_yz, g_z_0_zz_zz, g_z_xx, g_z_xy, g_z_xz, g_z_yy, g_z_yz, g_z_zz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_xx[k] = -g_z_xx[k] - g_z_0_z_xx[k] * ab_z + g_z_0_z_xxz[k];

                g_z_0_zz_xy[k] = -g_z_xy[k] - g_z_0_z_xy[k] * ab_z + g_z_0_z_xyz[k];

                g_z_0_zz_xz[k] = -g_z_xz[k] - g_z_0_z_xz[k] * ab_z + g_z_0_z_xzz[k];

                g_z_0_zz_yy[k] = -g_z_yy[k] - g_z_0_z_yy[k] * ab_z + g_z_0_z_yyz[k];

                g_z_0_zz_yz[k] = -g_z_yz[k] - g_z_0_z_yz[k] * ab_z + g_z_0_z_yzz[k];

                g_z_0_zz_zz[k] = -g_z_zz[k] - g_z_0_z_zz[k] * ab_z + g_z_0_z_zzz[k];
            }
        }
    }
}

} // erirec namespace

