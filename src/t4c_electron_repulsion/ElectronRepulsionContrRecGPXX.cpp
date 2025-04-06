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

#include "ElectronRepulsionContrRecGPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_hrr_electron_repulsion_gpxx(CSimdArray<double>& cbuffer,
                                     const size_t idx_gpxx,
                                     const size_t idx_fpxx,
                                     const size_t idx_fdxx,
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
            /// Set up components of auxilary buffer : FPSS

            const auto fp_off = idx_fpxx + i * dcomps + j;

            auto g_xxx_x = cbuffer.data(fp_off + 0 * ccomps * dcomps);

            auto g_xxx_y = cbuffer.data(fp_off + 1 * ccomps * dcomps);

            auto g_xxx_z = cbuffer.data(fp_off + 2 * ccomps * dcomps);

            auto g_xxy_x = cbuffer.data(fp_off + 3 * ccomps * dcomps);

            auto g_xxy_y = cbuffer.data(fp_off + 4 * ccomps * dcomps);

            auto g_xxy_z = cbuffer.data(fp_off + 5 * ccomps * dcomps);

            auto g_xxz_x = cbuffer.data(fp_off + 6 * ccomps * dcomps);

            auto g_xxz_y = cbuffer.data(fp_off + 7 * ccomps * dcomps);

            auto g_xxz_z = cbuffer.data(fp_off + 8 * ccomps * dcomps);

            auto g_xyy_x = cbuffer.data(fp_off + 9 * ccomps * dcomps);

            auto g_xyy_y = cbuffer.data(fp_off + 10 * ccomps * dcomps);

            auto g_xyy_z = cbuffer.data(fp_off + 11 * ccomps * dcomps);

            auto g_xyz_x = cbuffer.data(fp_off + 12 * ccomps * dcomps);

            auto g_xyz_y = cbuffer.data(fp_off + 13 * ccomps * dcomps);

            auto g_xyz_z = cbuffer.data(fp_off + 14 * ccomps * dcomps);

            auto g_xzz_x = cbuffer.data(fp_off + 15 * ccomps * dcomps);

            auto g_xzz_y = cbuffer.data(fp_off + 16 * ccomps * dcomps);

            auto g_xzz_z = cbuffer.data(fp_off + 17 * ccomps * dcomps);

            auto g_yyy_x = cbuffer.data(fp_off + 18 * ccomps * dcomps);

            auto g_yyy_y = cbuffer.data(fp_off + 19 * ccomps * dcomps);

            auto g_yyy_z = cbuffer.data(fp_off + 20 * ccomps * dcomps);

            auto g_yyz_x = cbuffer.data(fp_off + 21 * ccomps * dcomps);

            auto g_yyz_y = cbuffer.data(fp_off + 22 * ccomps * dcomps);

            auto g_yyz_z = cbuffer.data(fp_off + 23 * ccomps * dcomps);

            auto g_yzz_x = cbuffer.data(fp_off + 24 * ccomps * dcomps);

            auto g_yzz_y = cbuffer.data(fp_off + 25 * ccomps * dcomps);

            auto g_yzz_z = cbuffer.data(fp_off + 26 * ccomps * dcomps);

            auto g_zzz_x = cbuffer.data(fp_off + 27 * ccomps * dcomps);

            auto g_zzz_y = cbuffer.data(fp_off + 28 * ccomps * dcomps);

            auto g_zzz_z = cbuffer.data(fp_off + 29 * ccomps * dcomps);

            /// Set up components of auxilary buffer : FDSS

            const auto fd_off = idx_fdxx + i * dcomps + j;

            auto g_xxx_xx = cbuffer.data(fd_off + 0 * ccomps * dcomps);

            auto g_xxx_xy = cbuffer.data(fd_off + 1 * ccomps * dcomps);

            auto g_xxx_xz = cbuffer.data(fd_off + 2 * ccomps * dcomps);

            auto g_xxy_xx = cbuffer.data(fd_off + 6 * ccomps * dcomps);

            auto g_xxy_xy = cbuffer.data(fd_off + 7 * ccomps * dcomps);

            auto g_xxy_xz = cbuffer.data(fd_off + 8 * ccomps * dcomps);

            auto g_xxz_xx = cbuffer.data(fd_off + 12 * ccomps * dcomps);

            auto g_xxz_xy = cbuffer.data(fd_off + 13 * ccomps * dcomps);

            auto g_xxz_xz = cbuffer.data(fd_off + 14 * ccomps * dcomps);

            auto g_xyy_xx = cbuffer.data(fd_off + 18 * ccomps * dcomps);

            auto g_xyy_xy = cbuffer.data(fd_off + 19 * ccomps * dcomps);

            auto g_xyy_xz = cbuffer.data(fd_off + 20 * ccomps * dcomps);

            auto g_xyz_xx = cbuffer.data(fd_off + 24 * ccomps * dcomps);

            auto g_xyz_xy = cbuffer.data(fd_off + 25 * ccomps * dcomps);

            auto g_xyz_xz = cbuffer.data(fd_off + 26 * ccomps * dcomps);

            auto g_xzz_xx = cbuffer.data(fd_off + 30 * ccomps * dcomps);

            auto g_xzz_xy = cbuffer.data(fd_off + 31 * ccomps * dcomps);

            auto g_xzz_xz = cbuffer.data(fd_off + 32 * ccomps * dcomps);

            auto g_yyy_xx = cbuffer.data(fd_off + 36 * ccomps * dcomps);

            auto g_yyy_xy = cbuffer.data(fd_off + 37 * ccomps * dcomps);

            auto g_yyy_xz = cbuffer.data(fd_off + 38 * ccomps * dcomps);

            auto g_yyy_yy = cbuffer.data(fd_off + 39 * ccomps * dcomps);

            auto g_yyy_yz = cbuffer.data(fd_off + 40 * ccomps * dcomps);

            auto g_yyz_xx = cbuffer.data(fd_off + 42 * ccomps * dcomps);

            auto g_yyz_xy = cbuffer.data(fd_off + 43 * ccomps * dcomps);

            auto g_yyz_xz = cbuffer.data(fd_off + 44 * ccomps * dcomps);

            auto g_yyz_yy = cbuffer.data(fd_off + 45 * ccomps * dcomps);

            auto g_yyz_yz = cbuffer.data(fd_off + 46 * ccomps * dcomps);

            auto g_yzz_xx = cbuffer.data(fd_off + 48 * ccomps * dcomps);

            auto g_yzz_xy = cbuffer.data(fd_off + 49 * ccomps * dcomps);

            auto g_yzz_xz = cbuffer.data(fd_off + 50 * ccomps * dcomps);

            auto g_yzz_yy = cbuffer.data(fd_off + 51 * ccomps * dcomps);

            auto g_yzz_yz = cbuffer.data(fd_off + 52 * ccomps * dcomps);

            auto g_zzz_xx = cbuffer.data(fd_off + 54 * ccomps * dcomps);

            auto g_zzz_xy = cbuffer.data(fd_off + 55 * ccomps * dcomps);

            auto g_zzz_xz = cbuffer.data(fd_off + 56 * ccomps * dcomps);

            auto g_zzz_yy = cbuffer.data(fd_off + 57 * ccomps * dcomps);

            auto g_zzz_yz = cbuffer.data(fd_off + 58 * ccomps * dcomps);

            auto g_zzz_zz = cbuffer.data(fd_off + 59 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_gpxx

            const auto gp_off = idx_gpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xxxx_x = cbuffer.data(gp_off + 0 * ccomps * dcomps);

            auto g_xxxx_y = cbuffer.data(gp_off + 1 * ccomps * dcomps);

            auto g_xxxx_z = cbuffer.data(gp_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxx_x, g_xxx_xx, g_xxx_xy, g_xxx_xz, g_xxx_y, g_xxx_z, g_xxxx_x, g_xxxx_y, g_xxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxx_x[k] = -g_xxx_x[k] * ab_x + g_xxx_xx[k];

                g_xxxx_y[k] = -g_xxx_y[k] * ab_x + g_xxx_xy[k];

                g_xxxx_z[k] = -g_xxx_z[k] * ab_x + g_xxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xxxy_x = cbuffer.data(gp_off + 3 * ccomps * dcomps);

            auto g_xxxy_y = cbuffer.data(gp_off + 4 * ccomps * dcomps);

            auto g_xxxy_z = cbuffer.data(gp_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxy_x, g_xxxy_y, g_xxxy_z, g_xxy_x, g_xxy_xx, g_xxy_xy, g_xxy_xz, g_xxy_y, g_xxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxy_x[k] = -g_xxy_x[k] * ab_x + g_xxy_xx[k];

                g_xxxy_y[k] = -g_xxy_y[k] * ab_x + g_xxy_xy[k];

                g_xxxy_z[k] = -g_xxy_z[k] * ab_x + g_xxy_xz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xxxz_x = cbuffer.data(gp_off + 6 * ccomps * dcomps);

            auto g_xxxz_y = cbuffer.data(gp_off + 7 * ccomps * dcomps);

            auto g_xxxz_z = cbuffer.data(gp_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxxz_x, g_xxxz_y, g_xxxz_z, g_xxz_x, g_xxz_xx, g_xxz_xy, g_xxz_xz, g_xxz_y, g_xxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxxz_x[k] = -g_xxz_x[k] * ab_x + g_xxz_xx[k];

                g_xxxz_y[k] = -g_xxz_y[k] * ab_x + g_xxz_xy[k];

                g_xxxz_z[k] = -g_xxz_z[k] * ab_x + g_xxz_xz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_xxyy_x = cbuffer.data(gp_off + 9 * ccomps * dcomps);

            auto g_xxyy_y = cbuffer.data(gp_off + 10 * ccomps * dcomps);

            auto g_xxyy_z = cbuffer.data(gp_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxyy_x, g_xxyy_y, g_xxyy_z, g_xyy_x, g_xyy_xx, g_xyy_xy, g_xyy_xz, g_xyy_y, g_xyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxyy_x[k] = -g_xyy_x[k] * ab_x + g_xyy_xx[k];

                g_xxyy_y[k] = -g_xyy_y[k] * ab_x + g_xyy_xy[k];

                g_xxyy_z[k] = -g_xyy_z[k] * ab_x + g_xyy_xz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_xxyz_x = cbuffer.data(gp_off + 12 * ccomps * dcomps);

            auto g_xxyz_y = cbuffer.data(gp_off + 13 * ccomps * dcomps);

            auto g_xxyz_z = cbuffer.data(gp_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxyz_x, g_xxyz_y, g_xxyz_z, g_xyz_x, g_xyz_xx, g_xyz_xy, g_xyz_xz, g_xyz_y, g_xyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxyz_x[k] = -g_xyz_x[k] * ab_x + g_xyz_xx[k];

                g_xxyz_y[k] = -g_xyz_y[k] * ab_x + g_xyz_xy[k];

                g_xxyz_z[k] = -g_xyz_z[k] * ab_x + g_xyz_xz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_xxzz_x = cbuffer.data(gp_off + 15 * ccomps * dcomps);

            auto g_xxzz_y = cbuffer.data(gp_off + 16 * ccomps * dcomps);

            auto g_xxzz_z = cbuffer.data(gp_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xxzz_x, g_xxzz_y, g_xxzz_z, g_xzz_x, g_xzz_xx, g_xzz_xy, g_xzz_xz, g_xzz_y, g_xzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xxzz_x[k] = -g_xzz_x[k] * ab_x + g_xzz_xx[k];

                g_xxzz_y[k] = -g_xzz_y[k] * ab_x + g_xzz_xy[k];

                g_xxzz_z[k] = -g_xzz_z[k] * ab_x + g_xzz_xz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_xyyy_x = cbuffer.data(gp_off + 18 * ccomps * dcomps);

            auto g_xyyy_y = cbuffer.data(gp_off + 19 * ccomps * dcomps);

            auto g_xyyy_z = cbuffer.data(gp_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyyy_x, g_xyyy_y, g_xyyy_z, g_yyy_x, g_yyy_xx, g_yyy_xy, g_yyy_xz, g_yyy_y, g_yyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyyy_x[k] = -g_yyy_x[k] * ab_x + g_yyy_xx[k];

                g_xyyy_y[k] = -g_yyy_y[k] * ab_x + g_yyy_xy[k];

                g_xyyy_z[k] = -g_yyy_z[k] * ab_x + g_yyy_xz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_xyyz_x = cbuffer.data(gp_off + 21 * ccomps * dcomps);

            auto g_xyyz_y = cbuffer.data(gp_off + 22 * ccomps * dcomps);

            auto g_xyyz_z = cbuffer.data(gp_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyyz_x, g_xyyz_y, g_xyyz_z, g_yyz_x, g_yyz_xx, g_yyz_xy, g_yyz_xz, g_yyz_y, g_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyyz_x[k] = -g_yyz_x[k] * ab_x + g_yyz_xx[k];

                g_xyyz_y[k] = -g_yyz_y[k] * ab_x + g_yyz_xy[k];

                g_xyyz_z[k] = -g_yyz_z[k] * ab_x + g_yyz_xz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_xyzz_x = cbuffer.data(gp_off + 24 * ccomps * dcomps);

            auto g_xyzz_y = cbuffer.data(gp_off + 25 * ccomps * dcomps);

            auto g_xyzz_z = cbuffer.data(gp_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xyzz_x, g_xyzz_y, g_xyzz_z, g_yzz_x, g_yzz_xx, g_yzz_xy, g_yzz_xz, g_yzz_y, g_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xyzz_x[k] = -g_yzz_x[k] * ab_x + g_yzz_xx[k];

                g_xyzz_y[k] = -g_yzz_y[k] * ab_x + g_yzz_xy[k];

                g_xyzz_z[k] = -g_yzz_z[k] * ab_x + g_yzz_xz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_xzzz_x = cbuffer.data(gp_off + 27 * ccomps * dcomps);

            auto g_xzzz_y = cbuffer.data(gp_off + 28 * ccomps * dcomps);

            auto g_xzzz_z = cbuffer.data(gp_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xzzz_x, g_xzzz_y, g_xzzz_z, g_zzz_x, g_zzz_xx, g_zzz_xy, g_zzz_xz, g_zzz_y, g_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xzzz_x[k] = -g_zzz_x[k] * ab_x + g_zzz_xx[k];

                g_xzzz_y[k] = -g_zzz_y[k] * ab_x + g_zzz_xy[k];

                g_xzzz_z[k] = -g_zzz_z[k] * ab_x + g_zzz_xz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_yyyy_x = cbuffer.data(gp_off + 30 * ccomps * dcomps);

            auto g_yyyy_y = cbuffer.data(gp_off + 31 * ccomps * dcomps);

            auto g_yyyy_z = cbuffer.data(gp_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyy_x, g_yyy_xy, g_yyy_y, g_yyy_yy, g_yyy_yz, g_yyy_z, g_yyyy_x, g_yyyy_y, g_yyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyyy_x[k] = -g_yyy_x[k] * ab_y + g_yyy_xy[k];

                g_yyyy_y[k] = -g_yyy_y[k] * ab_y + g_yyy_yy[k];

                g_yyyy_z[k] = -g_yyy_z[k] * ab_y + g_yyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_yyyz_x = cbuffer.data(gp_off + 33 * ccomps * dcomps);

            auto g_yyyz_y = cbuffer.data(gp_off + 34 * ccomps * dcomps);

            auto g_yyyz_z = cbuffer.data(gp_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyyz_x, g_yyyz_y, g_yyyz_z, g_yyz_x, g_yyz_xy, g_yyz_y, g_yyz_yy, g_yyz_yz, g_yyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyyz_x[k] = -g_yyz_x[k] * ab_y + g_yyz_xy[k];

                g_yyyz_y[k] = -g_yyz_y[k] * ab_y + g_yyz_yy[k];

                g_yyyz_z[k] = -g_yyz_z[k] * ab_y + g_yyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_yyzz_x = cbuffer.data(gp_off + 36 * ccomps * dcomps);

            auto g_yyzz_y = cbuffer.data(gp_off + 37 * ccomps * dcomps);

            auto g_yyzz_z = cbuffer.data(gp_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_yyzz_x, g_yyzz_y, g_yyzz_z, g_yzz_x, g_yzz_xy, g_yzz_y, g_yzz_yy, g_yzz_yz, g_yzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yyzz_x[k] = -g_yzz_x[k] * ab_y + g_yzz_xy[k];

                g_yyzz_y[k] = -g_yzz_y[k] * ab_y + g_yzz_yy[k];

                g_yyzz_z[k] = -g_yzz_z[k] * ab_y + g_yzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_yzzz_x = cbuffer.data(gp_off + 39 * ccomps * dcomps);

            auto g_yzzz_y = cbuffer.data(gp_off + 40 * ccomps * dcomps);

            auto g_yzzz_z = cbuffer.data(gp_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_yzzz_x, g_yzzz_y, g_yzzz_z, g_zzz_x, g_zzz_xy, g_zzz_y, g_zzz_yy, g_zzz_yz, g_zzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yzzz_x[k] = -g_zzz_x[k] * ab_y + g_zzz_xy[k];

                g_yzzz_y[k] = -g_zzz_y[k] * ab_y + g_zzz_yy[k];

                g_yzzz_z[k] = -g_zzz_z[k] * ab_y + g_zzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_zzzz_x = cbuffer.data(gp_off + 42 * ccomps * dcomps);

            auto g_zzzz_y = cbuffer.data(gp_off + 43 * ccomps * dcomps);

            auto g_zzzz_z = cbuffer.data(gp_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_zzz_x, g_zzz_xz, g_zzz_y, g_zzz_yz, g_zzz_z, g_zzz_zz, g_zzzz_x, g_zzzz_y, g_zzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zzzz_x[k] = -g_zzz_x[k] * ab_z + g_zzz_xz[k];

                g_zzzz_y[k] = -g_zzz_y[k] * ab_z + g_zzz_yz[k];

                g_zzzz_z[k] = -g_zzz_z[k] * ab_z + g_zzz_zz[k];
            }
        }
    }
}

} // erirec namespace

