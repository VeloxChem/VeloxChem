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

#include "ElectronRepulsionGeom1000ContrRecFSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_fsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_fsxx,
                                            const size_t idx_dsxx,
                                            const size_t idx_geom_10_dsxx,
                                            const size_t idx_geom_10_dpxx,
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
            /// Set up components of auxilary buffer : DSSS

            const auto ds_off = idx_dsxx + i * dcomps + j;

            auto g_xx_0 = cbuffer.data(ds_off + 0 * ccomps * dcomps);

            auto g_yy_0 = cbuffer.data(ds_off + 3 * ccomps * dcomps);

            auto g_zz_0 = cbuffer.data(ds_off + 5 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DSSS

            const auto ds_geom_10_off = idx_geom_10_dsxx + i * dcomps + j;

            auto g_x_0_xx_0 = cbuffer.data(ds_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xy_0 = cbuffer.data(ds_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xz_0 = cbuffer.data(ds_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_yy_0 = cbuffer.data(ds_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_yz_0 = cbuffer.data(ds_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_zz_0 = cbuffer.data(ds_geom_10_off + 5 * ccomps * dcomps);

            auto g_y_0_xx_0 = cbuffer.data(ds_geom_10_off + 6 * ccomps * dcomps);

            auto g_y_0_xy_0 = cbuffer.data(ds_geom_10_off + 7 * ccomps * dcomps);

            auto g_y_0_xz_0 = cbuffer.data(ds_geom_10_off + 8 * ccomps * dcomps);

            auto g_y_0_yy_0 = cbuffer.data(ds_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_yz_0 = cbuffer.data(ds_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_zz_0 = cbuffer.data(ds_geom_10_off + 11 * ccomps * dcomps);

            auto g_z_0_xx_0 = cbuffer.data(ds_geom_10_off + 12 * ccomps * dcomps);

            auto g_z_0_xy_0 = cbuffer.data(ds_geom_10_off + 13 * ccomps * dcomps);

            auto g_z_0_xz_0 = cbuffer.data(ds_geom_10_off + 14 * ccomps * dcomps);

            auto g_z_0_yy_0 = cbuffer.data(ds_geom_10_off + 15 * ccomps * dcomps);

            auto g_z_0_yz_0 = cbuffer.data(ds_geom_10_off + 16 * ccomps * dcomps);

            auto g_z_0_zz_0 = cbuffer.data(ds_geom_10_off + 17 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DPSS

            const auto dp_geom_10_off = idx_geom_10_dpxx + i * dcomps + j;

            auto g_x_0_xx_x = cbuffer.data(dp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_y = cbuffer.data(dp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_z = cbuffer.data(dp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xy_y = cbuffer.data(dp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xz_y = cbuffer.data(dp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xz_z = cbuffer.data(dp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_yy_y = cbuffer.data(dp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_yz_y = cbuffer.data(dp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_zz_y = cbuffer.data(dp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_zz_z = cbuffer.data(dp_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_xx_x = cbuffer.data(dp_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_xy_x = cbuffer.data(dp_geom_10_off + 21 * ccomps * dcomps);

            auto g_y_0_xz_x = cbuffer.data(dp_geom_10_off + 24 * ccomps * dcomps);

            auto g_y_0_yy_x = cbuffer.data(dp_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_yy_y = cbuffer.data(dp_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_yy_z = cbuffer.data(dp_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_yz_x = cbuffer.data(dp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_yz_z = cbuffer.data(dp_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_zz_x = cbuffer.data(dp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_zz_z = cbuffer.data(dp_geom_10_off + 35 * ccomps * dcomps);

            auto g_z_0_xx_x = cbuffer.data(dp_geom_10_off + 36 * ccomps * dcomps);

            auto g_z_0_xy_x = cbuffer.data(dp_geom_10_off + 39 * ccomps * dcomps);

            auto g_z_0_xz_x = cbuffer.data(dp_geom_10_off + 42 * ccomps * dcomps);

            auto g_z_0_yy_x = cbuffer.data(dp_geom_10_off + 45 * ccomps * dcomps);

            auto g_z_0_yy_y = cbuffer.data(dp_geom_10_off + 46 * ccomps * dcomps);

            auto g_z_0_yz_x = cbuffer.data(dp_geom_10_off + 48 * ccomps * dcomps);

            auto g_z_0_yz_y = cbuffer.data(dp_geom_10_off + 49 * ccomps * dcomps);

            auto g_z_0_zz_x = cbuffer.data(dp_geom_10_off + 51 * ccomps * dcomps);

            auto g_z_0_zz_y = cbuffer.data(dp_geom_10_off + 52 * ccomps * dcomps);

            auto g_z_0_zz_z = cbuffer.data(dp_geom_10_off + 53 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fsxx

            const auto fs_geom_10_off = idx_geom_10_fsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_0, g_x_0_xx_x, g_x_0_xxx_0, g_xx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_0[k] = -g_xx_0[k] - g_x_0_xx_0[k] * ab_x + g_x_0_xx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_0, g_x_0_xx_y, g_x_0_xxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_0[k] = -g_x_0_xx_0[k] * ab_y + g_x_0_xx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_0, g_x_0_xx_z, g_x_0_xxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_0[k] = -g_x_0_xx_0[k] * ab_z + g_x_0_xx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xy_0, g_x_0_xy_y, g_x_0_xyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_0[k] = -g_x_0_xy_0[k] * ab_y + g_x_0_xy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyz_0, g_x_0_xz_0, g_x_0_xz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_0[k] = -g_x_0_xz_0[k] * ab_y + g_x_0_xz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xz_0, g_x_0_xz_z, g_x_0_xzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_0[k] = -g_x_0_xz_0[k] * ab_z + g_x_0_xz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_0, g_x_0_yy_y, g_x_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_0[k] = -g_x_0_yy_0[k] * ab_y + g_x_0_yy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyz_0, g_x_0_yz_0, g_x_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_0[k] = -g_x_0_yz_0[k] * ab_y + g_x_0_yz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzz_0, g_x_0_zz_0, g_x_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_0[k] = -g_x_0_zz_0[k] * ab_y + g_x_0_zz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_0, g_x_0_zz_z, g_x_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_0[k] = -g_x_0_zz_0[k] * ab_z + g_x_0_zz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_0, g_y_0_xx_x, g_y_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_0[k] = -g_y_0_xx_0[k] * ab_x + g_y_0_xx_x[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_0, g_y_0_xy_0, g_y_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_0[k] = -g_y_0_xy_0[k] * ab_x + g_y_0_xy_x[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxz_0, g_y_0_xz_0, g_y_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_0[k] = -g_y_0_xz_0[k] * ab_x + g_y_0_xz_x[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_0, g_y_0_yy_0, g_y_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_0[k] = -g_y_0_yy_0[k] * ab_x + g_y_0_yy_x[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyz_0, g_y_0_yz_0, g_y_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_0[k] = -g_y_0_yz_0[k] * ab_x + g_y_0_yz_x[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzz_0, g_y_0_zz_0, g_y_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_0[k] = -g_y_0_zz_0[k] * ab_x + g_y_0_zz_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_0, g_y_0_yy_y, g_y_0_yyy_0, g_yy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_0[k] = -g_yy_0[k] - g_y_0_yy_0[k] * ab_y + g_y_0_yy_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_0, g_y_0_yy_z, g_y_0_yyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_0[k] = -g_y_0_yy_0[k] * ab_z + g_y_0_yy_z[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yz_0, g_y_0_yz_z, g_y_0_yzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_0[k] = -g_y_0_yz_0[k] * ab_z + g_y_0_yz_z[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_0, g_y_0_zz_z, g_y_0_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_0[k] = -g_y_0_zz_0[k] * ab_z + g_y_0_zz_z[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_0, g_z_0_xx_x, g_z_0_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_0[k] = -g_z_0_xx_0[k] * ab_x + g_z_0_xx_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxy_0, g_z_0_xy_0, g_z_0_xy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_0[k] = -g_z_0_xy_0[k] * ab_x + g_z_0_xy_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_0, g_z_0_xz_0, g_z_0_xz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_0[k] = -g_z_0_xz_0[k] * ab_x + g_z_0_xz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyy_0, g_z_0_yy_0, g_z_0_yy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_0[k] = -g_z_0_yy_0[k] * ab_x + g_z_0_yy_x[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyz_0, g_z_0_yz_0, g_z_0_yz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_0[k] = -g_z_0_yz_0[k] * ab_x + g_z_0_yz_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_0, g_z_0_zz_0, g_z_0_zz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_0[k] = -g_z_0_zz_0[k] * ab_x + g_z_0_zz_x[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_0, g_z_0_yy_y, g_z_0_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_0[k] = -g_z_0_yy_0[k] * ab_y + g_z_0_yy_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_0, g_z_0_yz_0, g_z_0_yz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_0[k] = -g_z_0_yz_0[k] * ab_y + g_z_0_yz_y[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_0, g_z_0_zz_0, g_z_0_zz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_0[k] = -g_z_0_zz_0[k] * ab_y + g_z_0_zz_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_0, g_z_0_zz_z, g_z_0_zzz_0, g_zz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_0[k] = -g_zz_0[k] - g_z_0_zz_0[k] * ab_z + g_z_0_zz_z[k];
            }
        }
    }
}

} // erirec namespace

