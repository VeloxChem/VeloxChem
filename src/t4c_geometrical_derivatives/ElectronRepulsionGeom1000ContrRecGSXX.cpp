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

#include "ElectronRepulsionGeom1000ContrRecGSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_gsxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_gsxx,
                                            const size_t idx_fsxx,
                                            const size_t idx_geom_10_fsxx,
                                            const size_t idx_geom_10_fpxx,
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
            /// Set up components of auxilary buffer : FSSS

            const auto fs_off = idx_fsxx + i * dcomps + j;

            auto g_xxx_0 = cbuffer.data(fs_off + 0 * ccomps * dcomps);

            auto g_yyy_0 = cbuffer.data(fs_off + 6 * ccomps * dcomps);

            auto g_zzz_0 = cbuffer.data(fs_off + 9 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FSSS

            const auto fs_geom_10_off = idx_geom_10_fsxx + i * dcomps + j;

            auto g_x_0_xxx_0 = cbuffer.data(fs_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxy_0 = cbuffer.data(fs_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxz_0 = cbuffer.data(fs_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xyy_0 = cbuffer.data(fs_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xyz_0 = cbuffer.data(fs_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xzz_0 = cbuffer.data(fs_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_yyy_0 = cbuffer.data(fs_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_yyz_0 = cbuffer.data(fs_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_yzz_0 = cbuffer.data(fs_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_zzz_0 = cbuffer.data(fs_geom_10_off + 9 * ccomps * dcomps);

            auto g_y_0_xxx_0 = cbuffer.data(fs_geom_10_off + 10 * ccomps * dcomps);

            auto g_y_0_xxy_0 = cbuffer.data(fs_geom_10_off + 11 * ccomps * dcomps);

            auto g_y_0_xxz_0 = cbuffer.data(fs_geom_10_off + 12 * ccomps * dcomps);

            auto g_y_0_xyy_0 = cbuffer.data(fs_geom_10_off + 13 * ccomps * dcomps);

            auto g_y_0_xyz_0 = cbuffer.data(fs_geom_10_off + 14 * ccomps * dcomps);

            auto g_y_0_xzz_0 = cbuffer.data(fs_geom_10_off + 15 * ccomps * dcomps);

            auto g_y_0_yyy_0 = cbuffer.data(fs_geom_10_off + 16 * ccomps * dcomps);

            auto g_y_0_yyz_0 = cbuffer.data(fs_geom_10_off + 17 * ccomps * dcomps);

            auto g_y_0_yzz_0 = cbuffer.data(fs_geom_10_off + 18 * ccomps * dcomps);

            auto g_y_0_zzz_0 = cbuffer.data(fs_geom_10_off + 19 * ccomps * dcomps);

            auto g_z_0_xxx_0 = cbuffer.data(fs_geom_10_off + 20 * ccomps * dcomps);

            auto g_z_0_xxy_0 = cbuffer.data(fs_geom_10_off + 21 * ccomps * dcomps);

            auto g_z_0_xxz_0 = cbuffer.data(fs_geom_10_off + 22 * ccomps * dcomps);

            auto g_z_0_xyy_0 = cbuffer.data(fs_geom_10_off + 23 * ccomps * dcomps);

            auto g_z_0_xyz_0 = cbuffer.data(fs_geom_10_off + 24 * ccomps * dcomps);

            auto g_z_0_xzz_0 = cbuffer.data(fs_geom_10_off + 25 * ccomps * dcomps);

            auto g_z_0_yyy_0 = cbuffer.data(fs_geom_10_off + 26 * ccomps * dcomps);

            auto g_z_0_yyz_0 = cbuffer.data(fs_geom_10_off + 27 * ccomps * dcomps);

            auto g_z_0_yzz_0 = cbuffer.data(fs_geom_10_off + 28 * ccomps * dcomps);

            auto g_z_0_zzz_0 = cbuffer.data(fs_geom_10_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FPSS

            const auto fp_geom_10_off = idx_geom_10_fpxx + i * dcomps + j;

            auto g_x_0_xxx_x = cbuffer.data(fp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_y = cbuffer.data(fp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_z = cbuffer.data(fp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxy_y = cbuffer.data(fp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxz_y = cbuffer.data(fp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxz_z = cbuffer.data(fp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xyy_y = cbuffer.data(fp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xyz_y = cbuffer.data(fp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xzz_y = cbuffer.data(fp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xzz_z = cbuffer.data(fp_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_yyy_y = cbuffer.data(fp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_yyz_y = cbuffer.data(fp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_yzz_y = cbuffer.data(fp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_zzz_y = cbuffer.data(fp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_zzz_z = cbuffer.data(fp_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_xxx_x = cbuffer.data(fp_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_xxy_x = cbuffer.data(fp_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_xxz_x = cbuffer.data(fp_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_xyy_x = cbuffer.data(fp_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_xyz_x = cbuffer.data(fp_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_xzz_x = cbuffer.data(fp_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_yyy_x = cbuffer.data(fp_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_yyy_y = cbuffer.data(fp_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_yyy_z = cbuffer.data(fp_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_yyz_x = cbuffer.data(fp_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_yyz_z = cbuffer.data(fp_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_yzz_x = cbuffer.data(fp_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_yzz_z = cbuffer.data(fp_geom_10_off + 56 * ccomps * dcomps);

            auto g_y_0_zzz_x = cbuffer.data(fp_geom_10_off + 57 * ccomps * dcomps);

            auto g_y_0_zzz_z = cbuffer.data(fp_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_xxx_x = cbuffer.data(fp_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_xxy_x = cbuffer.data(fp_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_xxz_x = cbuffer.data(fp_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_xyy_x = cbuffer.data(fp_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_xyz_x = cbuffer.data(fp_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_xzz_x = cbuffer.data(fp_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_yyy_x = cbuffer.data(fp_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_yyy_y = cbuffer.data(fp_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_yyz_x = cbuffer.data(fp_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_yyz_y = cbuffer.data(fp_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_yzz_x = cbuffer.data(fp_geom_10_off + 84 * ccomps * dcomps);

            auto g_z_0_yzz_y = cbuffer.data(fp_geom_10_off + 85 * ccomps * dcomps);

            auto g_z_0_zzz_x = cbuffer.data(fp_geom_10_off + 87 * ccomps * dcomps);

            auto g_z_0_zzz_y = cbuffer.data(fp_geom_10_off + 88 * ccomps * dcomps);

            auto g_z_0_zzz_z = cbuffer.data(fp_geom_10_off + 89 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gsxx

            const auto gs_geom_10_off = idx_geom_10_gsxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_0, g_x_0_xxx_x, g_x_0_xxxx_0, g_xxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxx_0[k] = -g_xxx_0[k] - g_x_0_xxx_0[k] * ab_x + g_x_0_xxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_0, g_x_0_xxx_y, g_x_0_xxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxy_0[k] = -g_x_0_xxx_0[k] * ab_y + g_x_0_xxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxx_0, g_x_0_xxx_z, g_x_0_xxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxz_0[k] = -g_x_0_xxx_0[k] * ab_z + g_x_0_xxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxy_0, g_x_0_xxy_y, g_x_0_xxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyy_0[k] = -g_x_0_xxy_0[k] * ab_y + g_x_0_xxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyz_0, g_x_0_xxz_0, g_x_0_xxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyz_0[k] = -g_x_0_xxz_0[k] * ab_y + g_x_0_xxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxz_0, g_x_0_xxz_z, g_x_0_xxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzz_0[k] = -g_x_0_xxz_0[k] * ab_z + g_x_0_xxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyy_0, g_x_0_xyy_y, g_x_0_xyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyy_0[k] = -g_x_0_xyy_0[k] * ab_y + g_x_0_xyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyz_0, g_x_0_xyz_0, g_x_0_xyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyz_0[k] = -g_x_0_xyz_0[k] * ab_y + g_x_0_xyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzz_0, g_x_0_xzz_0, g_x_0_xzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzz_0[k] = -g_x_0_xzz_0[k] * ab_y + g_x_0_xzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzz_0, g_x_0_xzz_z, g_x_0_xzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzz_0[k] = -g_x_0_xzz_0[k] * ab_z + g_x_0_xzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyy_0, g_x_0_yyy_y, g_x_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyy_0[k] = -g_x_0_yyy_0[k] * ab_y + g_x_0_yyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyz_0, g_x_0_yyz_0, g_x_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyz_0[k] = -g_x_0_yyz_0[k] * ab_y + g_x_0_yyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzz_0, g_x_0_yzz_0, g_x_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzz_0[k] = -g_x_0_yzz_0[k] * ab_y + g_x_0_yzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzz_0, g_x_0_zzz_0, g_x_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzz_0[k] = -g_x_0_zzz_0[k] * ab_y + g_x_0_zzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzz_0, g_x_0_zzz_z, g_x_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzz_0[k] = -g_x_0_zzz_0[k] * ab_z + g_x_0_zzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxx_0, g_y_0_xxx_x, g_y_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxx_0[k] = -g_y_0_xxx_0[k] * ab_x + g_y_0_xxx_x[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxy_0, g_y_0_xxy_0, g_y_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxy_0[k] = -g_y_0_xxy_0[k] * ab_x + g_y_0_xxy_x[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxz_0, g_y_0_xxz_0, g_y_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxz_0[k] = -g_y_0_xxz_0[k] * ab_x + g_y_0_xxz_x[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyy_0, g_y_0_xyy_0, g_y_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyy_0[k] = -g_y_0_xyy_0[k] * ab_x + g_y_0_xyy_x[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyz_0, g_y_0_xyz_0, g_y_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyz_0[k] = -g_y_0_xyz_0[k] * ab_x + g_y_0_xyz_x[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzz_0, g_y_0_xzz_0, g_y_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzz_0[k] = -g_y_0_xzz_0[k] * ab_x + g_y_0_xzz_x[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyy_0, g_y_0_yyy_0, g_y_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyy_0[k] = -g_y_0_yyy_0[k] * ab_x + g_y_0_yyy_x[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyz_0, g_y_0_yyz_0, g_y_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyz_0[k] = -g_y_0_yyz_0[k] * ab_x + g_y_0_yyz_x[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzz_0, g_y_0_yzz_0, g_y_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzz_0[k] = -g_y_0_yzz_0[k] * ab_x + g_y_0_yzz_x[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzz_0, g_y_0_zzz_0, g_y_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzz_0[k] = -g_y_0_zzz_0[k] * ab_x + g_y_0_zzz_x[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_0, g_y_0_yyy_y, g_y_0_yyyy_0, g_yyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyy_0[k] = -g_yyy_0[k] - g_y_0_yyy_0[k] * ab_y + g_y_0_yyy_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyy_0, g_y_0_yyy_z, g_y_0_yyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyz_0[k] = -g_y_0_yyy_0[k] * ab_z + g_y_0_yyy_z[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyz_0, g_y_0_yyz_z, g_y_0_yyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzz_0[k] = -g_y_0_yyz_0[k] * ab_z + g_y_0_yyz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzz_0, g_y_0_yzz_z, g_y_0_yzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzz_0[k] = -g_y_0_yzz_0[k] * ab_z + g_y_0_yzz_z[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzz_0, g_y_0_zzz_z, g_y_0_zzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzz_0[k] = -g_y_0_zzz_0[k] * ab_z + g_y_0_zzz_z[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxx_0 = cbuffer.data(gs_geom_10_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxx_0, g_z_0_xxx_x, g_z_0_xxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxx_0[k] = -g_z_0_xxx_0[k] * ab_x + g_z_0_xxx_x[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxy_0 = cbuffer.data(gs_geom_10_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxy_0, g_z_0_xxy_0, g_z_0_xxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxy_0[k] = -g_z_0_xxy_0[k] * ab_x + g_z_0_xxy_x[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxz_0 = cbuffer.data(gs_geom_10_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxz_0, g_z_0_xxz_0, g_z_0_xxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxz_0[k] = -g_z_0_xxz_0[k] * ab_x + g_z_0_xxz_x[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyy_0 = cbuffer.data(gs_geom_10_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyy_0, g_z_0_xyy_0, g_z_0_xyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyy_0[k] = -g_z_0_xyy_0[k] * ab_x + g_z_0_xyy_x[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyz_0 = cbuffer.data(gs_geom_10_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyz_0, g_z_0_xyz_0, g_z_0_xyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyz_0[k] = -g_z_0_xyz_0[k] * ab_x + g_z_0_xyz_x[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzz_0 = cbuffer.data(gs_geom_10_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzz_0, g_z_0_xzz_0, g_z_0_xzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzz_0[k] = -g_z_0_xzz_0[k] * ab_x + g_z_0_xzz_x[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyy_0 = cbuffer.data(gs_geom_10_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyy_0, g_z_0_yyy_0, g_z_0_yyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyy_0[k] = -g_z_0_yyy_0[k] * ab_x + g_z_0_yyy_x[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyz_0 = cbuffer.data(gs_geom_10_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyz_0, g_z_0_yyz_0, g_z_0_yyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyz_0[k] = -g_z_0_yyz_0[k] * ab_x + g_z_0_yyz_x[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzz_0 = cbuffer.data(gs_geom_10_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzz_0, g_z_0_yzz_0, g_z_0_yzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzz_0[k] = -g_z_0_yzz_0[k] * ab_x + g_z_0_yzz_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzz_0 = cbuffer.data(gs_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzz_0, g_z_0_zzz_0, g_z_0_zzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzz_0[k] = -g_z_0_zzz_0[k] * ab_x + g_z_0_zzz_x[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyy_0 = cbuffer.data(gs_geom_10_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyy_0, g_z_0_yyy_y, g_z_0_yyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyy_0[k] = -g_z_0_yyy_0[k] * ab_y + g_z_0_yyy_y[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyz_0 = cbuffer.data(gs_geom_10_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyz_0, g_z_0_yyz_0, g_z_0_yyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyz_0[k] = -g_z_0_yyz_0[k] * ab_y + g_z_0_yyz_y[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzz_0 = cbuffer.data(gs_geom_10_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzz_0, g_z_0_yzz_0, g_z_0_yzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzz_0[k] = -g_z_0_yzz_0[k] * ab_y + g_z_0_yzz_y[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzz_0 = cbuffer.data(gs_geom_10_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzz_0, g_z_0_zzz_0, g_z_0_zzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzz_0[k] = -g_z_0_zzz_0[k] * ab_y + g_z_0_zzz_y[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzz_0 = cbuffer.data(gs_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzz_0, g_z_0_zzz_z, g_z_0_zzzz_0, g_zzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzz_0[k] = -g_zzz_0[k] - g_z_0_zzz_0[k] * ab_z + g_z_0_zzz_z[k];
            }
        }
    }
}

} // erirec namespace

